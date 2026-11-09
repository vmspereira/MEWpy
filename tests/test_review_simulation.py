"""Regression tests for the simulation/solvers findings of the 2026-11-09 review.

Each test reproduces one finding; see the docstrings for the finding id.
"""

import copy
import math
import multiprocessing
import pickle
import subprocess
import sys
import time
import types
import unittest
import warnings
from pathlib import Path
from unittest import mock

from mewpy.model.kinetic import ODEModel

DATA = Path(__file__).parent.joinpath("data")
EC_CORE_GZ = str(DATA.joinpath("e_coli_core.xml.gz"))
EC_CORE_XML = str(DATA.joinpath("e_coli_core.xml"))
R_BIOMASS = "R_BIOMASS_Ecoli_core_w_GAM"
C_BIOMASS = "BIOMASS_Ecoli_core_w_GAM"
G_BIOMASS = "Biomass_Ecoli_core"
WT_GROWTH = 0.8739


def reframed_sim(**kwargs):
    from reframed.io.sbml import load_cbmodel

    from mewpy.simulation import get_simulator

    return get_simulator(load_cbmodel(EC_CORE_GZ), **kwargs)


def cobra_sim(**kwargs):
    from cobra.io import read_sbml_model

    from mewpy.simulation import get_simulator

    return get_simulator(read_sbml_model(EC_CORE_GZ), **kwargs)


def germ_sim(**kwargs):
    from mewpy.io import read_sbml
    from mewpy.simulation import get_simulator

    return get_simulator(read_sbml(EC_CORE_XML, regulatory=False, warnings=False), **kwargs)


class SlowODE(ODEModel):
    """An ODE model whose right-hand side takes longer than the simulation timeout."""

    def get_ode(self, r_dict=None, params=None, factors=None):
        time.sleep(30)
        return super().get_ode(r_dict=r_dict, params=params, factors=factors)


class FakeKineticModel:
    metabolites = {"A": None}

    def __init__(self, rhs=None):
        self.rhs = rhs or (lambda t, y: y)

    def get_ode(self, r_dict=None, params=None, factors=None):
        return self.rhs


class TestCobraReferenceMethods(unittest.TestCase):
    """F1: COBRA lMOMA/ROOM used the knocked-out model as its own reference."""

    def test_lmoma_uses_wild_type_reference(self):
        from cobra.flux_analysis import moma, pfba

        sim = cobra_sim()
        ko = {"PGI": 0}
        fba_ko = sim.simulate(constraints=ko).fluxes[C_BIOMASS]
        lmoma = sim.simulate(method="lMOMA", constraints=ko).fluxes[C_BIOMASS]

        model = sim.model
        wt = pfba(model)
        with model:
            model.reactions.PGI.bounds = (0, 0)
            expected = moma(model, solution=wt, linear=True).fluxes[C_BIOMASS]
        self.assertGreater(abs(lmoma - fba_ko), 1e-3)
        self.assertAlmostEqual(lmoma, expected, places=3)

    def test_lmoma_honours_reference_argument(self):
        sim = cobra_sim()
        ref = sim.simulate(method="pFBA").fluxes
        zero_ref = {k: 0.0 for k in ref}
        r1 = sim.simulate(method="lMOMA", constraints={"PGI": 0}, reference=ref)
        r2 = sim.simulate(method="lMOMA", constraints={"PGI": 0}, reference=zero_ref)
        self.assertNotAlmostEqual(r1.objective_value, r2.objective_value, places=3)

    def test_room_uses_wild_type_reference(self):
        # A real ROOM run is a MILP that GLPK takes minutes to solve, so check the reference handed to cobra.
        import mewpy.simulation.cobra as csim

        sim = cobra_sim()
        wt_pgi = sim.simulate(method="pFBA").fluxes["PGI"]
        seen = {}

        def spy(model, solution=None, **kwargs):
            seen["PGI"] = None if solution is None else solution.fluxes["PGI"]
            return model.optimize()

        with mock.patch.object(csim, "room", spy):
            sim.simulate(method="ROOM", constraints={"PGI": 0})
        self.assertGreater(abs(wt_pgi), 1.0)
        self.assertIsNotNone(seen["PGI"], "ROOM must receive a reference solution")
        self.assertAlmostEqual(seen["PGI"], wt_pgi, places=4)


class TestReframedScaleFactor(unittest.TestCase):
    """F2: scalefactor was applied after the solver was built."""

    def test_scalefactor_preserves_objective(self):
        for reset in (True, False):
            sim = reframed_sim(reset_solver=reset)
            base = sim.simulate().objective_value
            scaled = sim.simulate(scalefactor=10).objective_value
            self.assertAlmostEqual(scaled, base, places=5, msg=f"reset_solver={reset}")
            self.assertEqual(sim.get_reaction_bounds("R_EX_glc__D_e"), (-10.0, math.inf))

    def test_scalefactor_with_slim(self):
        sim = reframed_sim()
        self.assertAlmostEqual(sim.simulate(scalefactor=10, slim=True), WT_GROWTH, places=4)


class TestReframedCachedSolver(unittest.TestCase):
    """F3: reset_solver=False broke copies, pickling, ROOM and lMOMA."""

    def test_copy_and_pickle_with_cached_solver(self):
        from mewpy.simulation import get_simulator

        sim = reframed_sim(reset_solver=False)
        sim.simulate()
        clone = get_simulator(sim)
        self.assertAlmostEqual(clone.simulate().objective_value, WT_GROWTH, places=4)
        restored = pickle.loads(pickle.dumps(sim))
        self.assertAlmostEqual(restored.simulate().objective_value, WT_GROWTH, places=4)
        self.assertIsNotNone(sim.solver, "the original keeps its cached solver")

    def test_wyield_same_with_and_without_cached_solver(self):
        from reframed.io.sbml import load_cbmodel

        from mewpy.optimization.evaluation import WYIELD
        from mewpy.problems import RKOProblem

        values = []
        for reset in (True, False):
            problem = RKOProblem(load_cbmodel(EC_CORE_GZ), [WYIELD(R_BIOMASS, "R_EX_succ_e")], reset_solver=reset)
            values.append(problem.evaluate_solution({"R_PFL": 0, "R_ACKr": 0, "R_LDH_D": 0}, decode=False)[0])
        self.assertGreater(values[0], 0.0)
        self.assertAlmostEqual(values[0], values[1], places=6)

    def test_room_repeatable_with_cached_solver(self):
        sim = reframed_sim(reset_solver=False)
        wt = reframed_sim().simulate(method="pFBA").fluxes
        reference = {r: wt[r] for r in ("R_PGI", "R_G6PDH2r", "R_PFK", "R_PYK", "R_CS")}
        first = sim.simulate(method="ROOM", constraints={"R_PGI": 0}, reference=reference)
        second = sim.simulate(method="ROOM", constraints={"R_PGI": 0}, reference=reference)
        self.assertGreater(first.objective_value, 0.5)
        self.assertAlmostEqual(first.objective_value, second.objective_value, places=4)

    def test_fba_fluxes_are_reactions_after_pfba_with_cached_solver(self):
        sim = reframed_sim(reset_solver=False)
        sim.simulate(method="pFBA")
        self.assertEqual(set(sim.simulate().fluxes), set(sim.reactions))

    def test_lmoma_follows_new_reference_with_cached_solver(self):
        wt = reframed_sim().simulate(method="pFBA").fluxes
        zero_ref = {k: 0.0 for k in wt}
        cached = reframed_sim(reset_solver=False)
        cached.simulate(method="lMOMA", constraints={"R_PGI": 0}, reference=wt)
        second = cached.simulate(method="lMOMA", constraints={"R_PGI": 0}, reference=zero_ref)
        fresh = reframed_sim().simulate(method="lMOMA", constraints={"R_PGI": 0}, reference=zero_ref)
        self.assertAlmostEqual(second.objective_value, fresh.objective_value, places=4)


class TestSCIPSolverStages(unittest.TestCase):
    """F4: SCIP wrapper modifications after solve, empty constraints, pending cache."""

    def setUp(self):
        from mewpy.solvers.pyscipopt_solver import PySCIPOptSolver

        self.solver = PySCIPOptSolver(reframed_sim())

    def test_add_and_remove_constraint_after_solve(self):
        s = self.solver
        self.assertAlmostEqual(s.solve({R_BIOMASS: 1}, minimize=False).fobj, WT_GROWTH, places=4)
        s.add_constraint("cap", {R_BIOMASS: 1}, "<", 0.5)
        self.assertAlmostEqual(s.solve({R_BIOMASS: 1}, minimize=False).fobj, 0.5, places=5)
        s.remove_constraint("cap")
        s.add_variable("extra", 0, 1)
        self.assertAlmostEqual(s.solve({R_BIOMASS: 1}, minimize=False).fobj, WT_GROWTH, places=4)

    def test_empty_constraint_is_added(self):
        from mewpy.solvers import Status

        s = self.solver
        s.add_constraint("empty", {}, "<", -1)
        self.assertEqual(s.solve({R_BIOMASS: 1}, minimize=False).status, Status.INFEASIBLE)

    def test_solution_pool_returns_list(self):
        from mewpy.solvers.pyscipopt_solver import PySCIPOptSolver
        from mewpy.solvers.solver import VarType

        s = PySCIPOptSolver()
        for v in ("x", "y", "z"):
            s.add_variable(v, 0, 1, vartype=VarType.BINARY)
        s.add_constraint("pick_one", {"x": 1, "y": 1, "z": 1}, ">", 1)
        pool = s.solve({"x": 1, "y": 1, "z": 1}, minimize=True, get_values=["x", "y"], pool_size=5, pool_gap=0.5)
        self.assertIsInstance(pool, list)
        self.assertGreaterEqual(len(pool), 1)
        self.assertTrue(all(abs(sol.fobj - 1.0) < 1e-6 for sol in pool))
        self.assertEqual(set(pool[0].values), {"x", "y"})

    def test_solve_flushes_pending_rows(self):
        s = self.solver
        s.add_variable("aux", 0, 0.25, update=False)
        s.add_constraint("link", {R_BIOMASS: 1, "aux": -1}, "=", 0, update=False)
        self.assertAlmostEqual(s.solve({R_BIOMASS: 1}, minimize=False).fobj, 0.25, places=5)


class TestOptLangSolver(unittest.TestCase):
    """F5 (infinite bounds abort GLPK), set_variable_bounds truthiness and set_parameter no-op."""

    def test_infinite_bounds_do_not_abort(self):
        code = (
            "import warnings; warnings.filterwarnings('ignore')\n"
            "from reframed.io.sbml import load_cbmodel\n"
            "from mewpy.simulation import get_simulator\n"
            "from mewpy.solvers.optlang_solver import OptLangSolver\n"
            f"s = OptLangSolver(get_simulator(load_cbmodel({EC_CORE_GZ!r})))\n"
            f"obj = {{{R_BIOMASS!r}: 1}}\n"
            "print(s.solve(obj, minimize=False).fobj)\n"
            "print(s.solve(obj, minimize=False, constraints={'R_EX_glc__D_e': (-5, float('inf'))}).fobj)\n"
            "print(s.solve(obj, minimize=False).fobj)\n"
        )
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stderr[-1000:])
        wt, limited, restored = (float(x) for x in result.stdout.split())
        self.assertAlmostEqual(wt, WT_GROWTH, places=4)
        self.assertLess(limited, wt - 0.1)
        self.assertAlmostEqual(restored, WT_GROWTH, places=4)

    def test_set_variable_bounds_to_zero(self):
        from mewpy.solvers.optlang_solver import OptLangSolver

        s = OptLangSolver(cobra_sim())
        s.set_variable_bounds(C_BIOMASS, 0, 0)
        self.assertAlmostEqual(s.solve({C_BIOMASS: 1}, minimize=False).fobj, 0.0, places=6)

    def test_set_parameter_time_limit(self):
        from mewpy.solvers.optlang_solver import OptLangSolver
        from mewpy.solvers.solver import Parameter

        s = OptLangSolver()
        s.set_parameter(Parameter.TIME_LIMIT, 5)
        self.assertEqual(s.problem.configuration.timeout, 5)


class TestGermSimulator(unittest.TestCase):
    """F6 (constraint precedence) and F7 (essential genes as objects)."""

    def test_simulation_constraints_override_persistent(self):
        sim = germ_sim(constraints={"PGI": (-1000, 1000)})
        res = sim.simulate(constraints={"PGI": 0})
        self.assertAlmostEqual(res.fluxes["PGI"], 0.0, places=6)
        self.assertLess(res.objective_value, WT_GROWTH - 0.005)

    def test_essential_genes_are_ids(self):
        sim = germ_sim()
        essential = sim.essential_genes()
        self.assertGreater(len(essential), 0)
        self.assertTrue(all(isinstance(g, str) for g in essential))
        self.assertTrue(set(essential) <= set(sim.genes))


class TestCobraObjective(unittest.TestCase):
    """F8: COBRA dict objectives kept only the first key; slim ignored the sense."""

    def setUp(self):
        self.sim = cobra_sim()

    def test_multi_term_objective(self):
        res = self.sim.simulate(objective={C_BIOMASS: 1, "EX_ac_e": 1})
        self.assertAlmostEqual(res.objective_value, 20.0, places=4)

    def test_negative_coefficient(self):
        res = self.sim.simulate(objective={"EX_glc__D_e": -1})
        self.assertAlmostEqual(res.objective_value, 10.0, places=4)

    def test_slim_minimize(self):
        self.assertAlmostEqual(self.sim.simulate(objective={"EX_ac_e": 1}, maximize=False, slim=True), 0.0, places=6)
        self.assertAlmostEqual(self.sim.simulate(objective={"EX_ac_e": 1}, slim=True), 20.0, places=4)

    def test_objective_restored(self):
        self.sim.simulate(objective={"EX_ac_e": 1}, maximize=False, slim=True)
        self.assertEqual(self.sim.objective, {C_BIOMASS: 1.0})
        self.assertAlmostEqual(self.sim.simulate().objective_value, WT_GROWTH, places=4)


class TestInfeasibleHandling(unittest.TestCase):
    """F10: infeasible problems raised on COBRA pFBA/lMOMA and GERM pFBA."""

    def test_cobra_infeasible_returns_status(self):
        from mewpy.simulation import SStatus

        sim = cobra_sim()
        for method in ("FBA", "pFBA", "lMOMA"):
            res = sim.simulate(method=method, constraints={"ATPM": (500, 500)})
            self.assertEqual(res.status, SStatus.INFEASIBLE, method)
            self.assertIsNone(res.objective_value, method)
        self.assertAlmostEqual(sim.simulate().objective_value, WT_GROWTH, places=4)

    def test_germ_pfba_infeasible_returns_status(self):
        from mewpy.simulation import SStatus

        res = germ_sim().simulate(method="pFBA", constraints={"ATPM": (500, 500)})
        self.assertEqual(res.status, SStatus.INFEASIBLE)


class TestKinetic(unittest.TestCase):
    """F11: kinetic failure paths."""

    def test_negative_concentration_returns_error(self):
        from mewpy.simulation import kinetic
        from mewpy.solvers import ODEStatus

        fake = types.SimpleNamespace(solve=lambda y0, ts: ([-1.0], [0.0], [[-1.0]]))
        with mock.patch.object(kinetic, "ode_solver_instance", return_value=fake):
            result = kinetic.kinetic_solve(FakeKineticModel(), [1.0], [0, 1])
        self.assertEqual(len(result), 5)
        self.assertEqual(result[0], ODEStatus.ERROR)

    def test_failed_integration_is_an_error(self):
        from mewpy.simulation import kinetic
        from mewpy.solvers import ODEStatus
        from mewpy.solvers.ode import ODEMethod
        from mewpy.solvers.scipy_solver import ScipySolver

        model = FakeKineticModel(rhs=lambda t, y: [y[0] ** 2])  # blows up at t=1
        with mock.patch.object(kinetic, "ode_solver_instance", lambda f, m: ScipySolver(f, ODEMethod.RK45)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = kinetic.kinetic_solve(model, [1.0], [0, 10])
        self.assertEqual(result[0], ODEStatus.ERROR)

    def test_timeout_terminates_worker(self):
        from mewpy.simulation.kinetic import KineticSimulation

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = KineticSimulation(SlowODE("slow"), timeout=1).simulate()
            alive = [p for p in multiprocessing.active_children() if type(p).__name__ == "KineticThread"]
            self.assertEqual(alive, [])
            self.assertNotEqual(str(result.status), "ODEStatus.OPTIMAL")
        finally:
            for p in multiprocessing.active_children():
                if type(p).__name__ == "KineticThread":
                    p.terminate()


class TestResultConstraintsAreCopies(unittest.TestCase):
    """F-OPT-5 (simulation side): results must not alias the caller's constraints dict."""

    def check(self, sim, ko):
        constraints = dict(ko)
        res = sim.simulate(constraints=constraints)
        res.simulation_constraints["__marker__"] = 1
        res.envcond["__marker__"] = 1
        res.model_constraints["__marker__"] = 1
        self.assertNotIn("__marker__", constraints)
        self.assertNotIn("__marker__", sim.environmental_conditions)
        self.assertNotIn("__marker__", sim.simulate().model_constraints)

    def test_reframed(self):
        self.check(reframed_sim(envcond={"R_EX_glc__D_e": (-10, 1000)}), {"R_PGI": 0})

    def test_cobra(self):
        self.check(cobra_sim(envcond={"EX_glc__D_e": (-10, 1000)}), {"PGI": 0})

    def test_germ(self):
        self.check(germ_sim(envcond={"EX_glc__D_e": (-10, 1000)}), {"PGI": 0})


class TestSmallItems(unittest.TestCase):
    """Items left out of the 12-finding cap."""

    def test_reframed_moma_warns_when_unavailable(self):
        import mewpy.simulation.reframed as rsim

        if rsim.MOMA is not rsim.lMOMA:
            self.skipTest("reframed provides MOMA")
        with self.assertWarns(UserWarning):
            reframed_sim().simulate(method="MOMA", constraints={"R_PGI": 0})

    def test_set_default_solver_glpk(self):
        import mewpy.simulation as msim
        import mewpy.solvers as msolvers

        previous_sim, previous_lp = msim.get_default_solver(), msolvers.get_default_solver()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                msim.set_default_solver("glpk")
            self.assertEqual(msim.get_default_solver(), "glpk")
            self.assertEqual(msolvers.get_default_solver(), previous_lp)
        finally:
            msim.set_default_solver(previous_sim)
            msolvers.set_default_solver(previous_lp)

    def test_get_container_unknown_model(self):
        from mewpy.simulation import get_container

        with self.assertRaises(ValueError):
            get_container(object())

    def test_hybrid_mapping_errors(self):
        from mewpy.simulation.hybrid import HybridSimulation

        empty = types.SimpleNamespace(get_mapping=lambda: {})
        with self.assertRaises(ValueError):
            HybridSimulation.mapping_conversion(empty, {"a": 1.0})
        with self.assertRaises(ValueError):
            HybridSimulation.mapping_bounds(empty, {"a": 1.0}, {"a": 2.0})

    def test_hybrid_nsamples_failures_and_nans(self):
        from mewpy.simulation import hybrid

        outcomes = iter([RuntimeError("boom"), {"a": float("nan")}, {"a": 1.0}])

        class FakeKSim:
            def __init__(self, *args, **kwargs):
                pass

            def simulate(self, parameters=None):
                outcome = next(outcomes)
                if isinstance(outcome, Exception):
                    raise outcome
                return types.SimpleNamespace(fluxes=outcome)

        holder = types.SimpleNamespace(get_kinetic_model=lambda: None, t_points=[0, 1], timeout=0)
        with mock.patch.object(hybrid, "KineticSimulation", FakeKSim):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                df = hybrid.HybridSimulation.nsamples(holder, {"v": 1.0}, n=3)
        self.assertEqual(len(df), 1)

    def test_simulator_deepcopy_keeps_results(self):
        sim = reframed_sim()
        clone = copy.deepcopy(sim)
        self.assertAlmostEqual(clone.simulate().objective_value, WT_GROWTH, places=4)


if __name__ == "__main__":
    unittest.main()
