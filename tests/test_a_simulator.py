import re
import unittest
import warnings
from functools import lru_cache
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"
EC_CORE_MODEL2 = Path(__file__).parent.joinpath("data", "e_coli_core.xml")
MIN_GROWTH = 0.1
# Orth et al. (2010): e_coli_core grows at 0.8739 /h on 10 mmol/gDW/h glucose, aerobic
AEROBIC_GROWTH = 0.8739
# ROOM keeps a flux unchanged within ref +- (delta * |ref| + epsilon) (reframed and COBRApy defaults)
ROOM_DELTA, ROOM_EPSILON = 0.03, 0.001
# GECKO single-pool proteins: the P33312 knockout stops growth (0.0 with GLPK through OptLangSolver, and with
# SCIP); Q12184 carries no flux in the wild type
GECKO_ESSENTIAL, GECKO_UNUSED = "P33312", "Q12184"


def _strip(identifier):
    """The identifier without its SBML prefix (REFRAMED R_/M_/G_), as COBRApy names it."""
    return re.sub(r"^[RMG]_", "", identifier)


@lru_cache(maxsize=None)
def _cobra_model(path):
    """COBRApy reads the same SBML file: the external oracle. Changes are made in `with model:` blocks."""
    from cobra.io import read_sbml_model

    return read_sbml_model(str(path))


@lru_cache(maxsize=None)
def _cobra_reference(path):
    """Wild-type pFBA fluxes (COBRApy): the reference of MOMA, lMOMA and ROOM."""
    from cobra.flux_analysis import pfba

    return pfba(_cobra_model(path)).fluxes.to_dict()


@lru_cache(maxsize=None)
def _cobra_lmoma(path, knockout):
    """COBRApy linear MOMA fluxes of a knockout, against the wild-type pFBA reference."""
    import pandas as pd
    from cobra.core import Solution
    from cobra.flux_analysis import moma

    model = _cobra_model(path)
    reference = Solution(None, "optimal", fluxes=pd.Series(_cobra_reference(path)))
    with model:
        model.reactions.get_by_id(knockout).knock_out()
        return moma(model, solution=reference, linear=True).fluxes.to_dict()


@lru_cache(maxsize=None)
def _cobra_room_changes(path, knockout):
    """Number of fluxes that COBRApy's ROOM (MILP) moves out of the reference window after a knockout."""
    import pandas as pd
    from cobra.core import Solution
    from cobra.flux_analysis import room

    model = _cobra_model(path)
    ref = _cobra_reference(path)
    with model:
        model.reactions.get_by_id(knockout).knock_out()
        fluxes = room(model, solution=Solution(None, "optimal", fluxes=pd.Series(ref)), linear=False).fluxes
    return _room_changes(ref, {r: fluxes[r] for r in ref})


def _room_changes(reference, fluxes):
    return sorted(r for r, v in reference.items() if abs(fluxes[r] - v) > ROOM_DELTA * abs(v) + ROOM_EPSILON + 1e-6)


def _cobra_envelope(path, r_x, r_y, steps):
    """Range of r_y at `steps` values of r_x, from its minimum to its maximum (COBRApy FVA)."""
    from cobra.flux_analysis import flux_variability_analysis

    model = _cobra_model(path)
    x_range = flux_variability_analysis(model, reaction_list=[r_x], fraction_of_optimum=0)
    xs = np.linspace(x_range.loc[r_x, "minimum"], x_range.loc[r_x, "maximum"], steps)
    ymins, ymaxs = [], []
    for x in xs:
        with model:
            model.reactions.get_by_id(r_x).bounds = (x, x)
            y = flux_variability_analysis(model, reaction_list=[r_y], fraction_of_optimum=0)
        ymins.append(y.loc[r_y, "minimum"])
        ymaxs.append(y.loc[r_y, "maximum"])
    return xs, np.array(ymins), np.array(ymaxs)


class TestReframedSimul(unittest.TestCase):
    """Tests the REFRAMED Simulator"""

    MODEL_FILE = EC_CORE_MODEL

    def setUp(self):
        """Set up
        Loads a model
        """
        from reframed.io.sbml import load_cbmodel

        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.simulation import get_simulator

        self.simul = get_simulator(model)
        self.BIOMASS_ID = model.biomass_reaction
        self.SUCC = "R_EX_succ_e"

    @property
    def native(self):
        """COBRApy reaction id -> the simulator's reaction id"""
        return {_strip(r): r for r in self.simul.reactions}

    @property
    def reference(self):
        """The COBRApy wild-type pFBA fluxes, keyed by the simulator's reaction ids"""
        native = self.native
        return {native[r]: v for r, v in _cobra_reference(self.MODEL_FILE).items()}

    def distance(self, fluxes, power=1):
        """Distance of a flux distribution (simulator ids) to the reference"""
        native = self.native
        return sum(abs(fluxes[native[r]] - v) ** power for r, v in _cobra_reference(self.MODEL_FILE).items())

    def test_essential_reactions(self):
        """Essential reactions are those COBRApy finds (knockout growth below 1% of the wild type)"""
        from cobra.flux_analysis import find_essential_reactions

        essential = {_strip(r) for r in self.simul.essential_reactions()}
        expected = {r.id for r in find_essential_reactions(_cobra_model(self.MODEL_FILE), processes=1)}
        self.assertEqual(essential, expected)

    def test_essential_genes(self):
        """Essential genes are those COBRApy finds; b1779 is the only gene of the essential GAPD"""
        from cobra.flux_analysis import find_essential_genes

        essential = {_strip(g) for g in self.simul.essential_genes()}
        expected = {g.id for g in find_essential_genes(_cobra_model(self.MODEL_FILE), processes=1)}
        self.assertEqual(essential, expected)
        self.assertIn("b1779", essential)

    def test_uptake_reactions(self):
        """Uptake reactions are the exchanges open for uptake: the medium"""
        uptake = {_strip(r) for r in self.simul.get_uptake_reactions()}
        medium = {r.id for r in _cobra_model(self.MODEL_FILE).exchanges if r.lower_bound < 0}
        self.assertEqual(uptake, medium)

    def test_transport_reactions(self):
        """Transport reactions move metabolites between compartments; exchanges are not transport reactions"""
        model = _cobra_model(self.MODEL_FILE)
        transport = {_strip(r) for r in self.simul.get_transport_reactions()}
        for r_id in sorted(transport):
            rxn = model.reactions.get_by_id(r_id)
            self.assertFalse(rxn.boundary, f"{r_id} is an exchange reaction")
            self.assertGreater(len(rxn.compartments), 1, f"{r_id} has a single compartment")
        # one metabolite crossing the membrane (e.g. O2t: o2_e <=> o2_c) is a transport under any definition
        uniporters = {
            r.id
            for r in model.reactions
            if len(r.reactants) == 1
            and len(r.products) == 1
            and r.reactants[0].compartment != r.products[0].compartment
        }
        self.assertIn("O2t", uniporters)
        self.assertLessEqual(uniporters, transport)

    def test_fba(self):
        """FBA: the textbook aerobic growth on glucose, as computed by COBRApy"""
        from mewpy.simulation import SStatus

        res = self.simul.simulate()
        self.assertEqual(res.status, SStatus.OPTIMAL)
        self.assertAlmostEqual(res.objective_value, _cobra_model(self.MODEL_FILE).slim_optimize(), places=6)
        self.assertAlmostEqual(res.objective_value, AEROBIC_GROWTH, places=4)

    def test_pfba(self):
        """pFBA: optimal growth (REFRAMED keeps 99.9%) with the minimal total flux, as COBRApy's pfba"""
        from cobra.flux_analysis import pfba

        res = self.simul.simulate(method="pFBA")
        wild_type = _cobra_model(self.MODEL_FILE).slim_optimize()
        growth = res.fluxes[self.BIOMASS_ID]
        self.assertGreaterEqual(growth, 0.999 * wild_type - 1e-6)
        self.assertLessEqual(growth, wild_type + 1e-6)
        expected = pfba(_cobra_model(self.MODEL_FILE), fraction_of_optimum=min(1.0, growth / wild_type)).objective_value
        self.assertAlmostEqual(sum(abs(v) for v in res.fluxes.values()), expected, delta=1e-6 * expected)

    def test_moma(self):
        """MOMA of a PGI knockout: no feasible flux distribution is closer (L2) to the wild type
        Note: MOMA requires a QP-capable solver (CPLEX, Gurobi)
        """
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                res = self.simul.simulate(method="MOMA", constraints={self.native["PGI"]: 0}, reference=self.reference)
        except Exception as e:
            if "QP-capable solver" in str(e) or "SolverNotFound" in str(type(e).__name__):
                self.skipTest("MOMA requires a QP-capable solver (CPLEX, Gurobi)")
            else:
                raise
        self.assertAlmostEqual(res.fluxes[self.native["PGI"]], 0.0, places=6)
        lmoma = _cobra_lmoma(self.MODEL_FILE, "PGI")
        lmoma = {self.native[r]: v for r, v in lmoma.items()}
        if any("running lMOMA instead" in str(w.message) for w in caught):
            # without a quadratic MOMA, reframed runs lMOMA and says so
            self.assertAlmostEqual(self.distance(res.fluxes), self.distance(lmoma), delta=1e-6)
        else:
            self.assertLessEqual(self.distance(res.fluxes, power=2), self.distance(lmoma, power=2) + 1e-6)

    def test_lmoma(self):
        """lMOMA: the L1 distance to the wild type after a PGI knockout equals COBRApy's linear MOMA"""
        native = self.native
        res = self.simul.simulate(method="lMOMA", constraints={native["PGI"]: 0}, reference=self.reference)
        self.assertAlmostEqual(res.fluxes[native["PGI"]], 0.0, places=6)
        expected = self.distance({native[r]: v for r, v in _cobra_lmoma(self.MODEL_FILE, "PGI").items()})
        self.assertGreater(expected, 1.0)
        self.assertAlmostEqual(self.distance(res.fluxes), expected, delta=1e-6 * expected)
        self.assertLessEqual(res.fluxes[self.BIOMASS_ID], _cobra_model(self.MODEL_FILE).slim_optimize() + 1e-6)
        # without changes lMOMA reproduces the reference
        res = self.simul.simulate(method="lMOMA", reference=self.reference)
        self.assertAlmostEqual(self.distance(res.fluxes), 0.0, places=6)

    def check_room(self):
        """ROOM of a PYK knockout moves as few fluxes out of the reference window as COBRApy's ROOM"""
        from mewpy.simulation import SStatus

        native = self.native
        res = self.simul.simulate(method="ROOM", constraints={native["PYK"]: 0}, reference=self.reference)
        self.assertEqual(res.status, SStatus.OPTIMAL)
        self.assertAlmostEqual(res.fluxes[native["PYK"]], 0.0, places=6)
        reference = _cobra_reference(self.MODEL_FILE)
        changed = _room_changes(reference, {r: res.fluxes[native[r]] for r in reference})
        expected = _cobra_room_changes(self.MODEL_FILE, "PYK")
        self.assertIn("PYK", expected)
        self.assertEqual(len(changed), len(expected), f"changed {changed}, COBRApy changes {expected}")

    @pytest.mark.xfail(
        strict=True,
        reason="reframed ROOM uses big-M 1e6: binaries of ~1e-7 (SCIP integrality tolerance) let fluxes leave "
        "the ROOM window uncounted (objective 4, 9 fluxes changed; COBRApy's MILP optimum is 5)",
    )
    def test_room(self):
        """Tests ROOM"""
        self.check_room()

    def test_FVA(self):
        """FVA at 90% of the optimum (the default) matches COBRApy and brackets the FBA flux"""
        from cobra.flux_analysis import flux_variability_analysis

        res = self.simul.FVA(reactions=[self.SUCC])
        succ = _strip(self.SUCC)
        expected = flux_variability_analysis(
            _cobra_model(self.MODEL_FILE), reaction_list=[succ], fraction_of_optimum=0.9
        )
        low, high = res[self.SUCC]
        self.assertAlmostEqual(low, expected.loc[succ, "minimum"], places=6)
        self.assertAlmostEqual(high, expected.loc[succ, "maximum"], places=6)
        self.assertGreater(high, low + 1.0)
        flux = self.simul.simulate().fluxes[self.SUCC]
        self.assertLessEqual(low - 1e-6, flux)
        self.assertLessEqual(flux, high + 1e-6)

    def test_envelope(self):
        """The plotted envelope is the succinate range at each growth rate (COBRApy FVA)"""
        import matplotlib.pyplot as plt

        from mewpy.visualization.envelope import plot_flux_envelope

        ax = plot_flux_envelope(self.simul, self.BIOMASS_ID, self.SUCC)
        try:
            xs, ymins = ax.lines[0].get_data()
            _, ymaxs = ax.lines[1].get_data()
        finally:
            plt.close(ax.figure)
        exp_x, exp_min, exp_max = _cobra_envelope(self.MODEL_FILE, _strip(self.BIOMASS_ID), _strip(self.SUCC), len(xs))
        np.testing.assert_allclose(xs, exp_x, atol=1e-6)
        np.testing.assert_allclose(ymins, exp_min, atol=1e-6)
        np.testing.assert_allclose(ymaxs, exp_max, atol=1e-6)
        self.assertAlmostEqual(exp_x[-1], AEROBIC_GROWTH, places=4)

    def test_solver(self):
        """A solver built from the simulator holds the model: steady state, bounds and the FBA optimum"""
        from cobra.util.array import create_stoichiometric_matrix

        from mewpy.solvers import solver_instance
        from mewpy.solvers.solution import Status

        model = _cobra_model(self.MODEL_FILE)
        stoichiometry = create_stoichiometric_matrix(model)
        native = self.native
        for objective in (None, {self.BIOMASS_ID: 1}):
            solver = solver_instance(self.simul)
            solution = solver.solve(objective, minimize=False) if objective else solver.solve()
            self.assertEqual(solution.status, Status.OPTIMAL)
            v = np.array([solution.values[native[r.id]] for r in model.reactions])
            np.testing.assert_allclose(stoichiometry @ v, 0, atol=1e-6)
            for r, flux in zip(model.reactions, v):
                self.assertTrue(r.lower_bound - 1e-6 <= flux <= r.upper_bound + 1e-6, r.id)
        self.assertAlmostEqual(solution.fobj, model.slim_optimize(), places=6)


class TestCobra(TestReframedSimul):
    """Tests COBRApy Simulator"""

    MODEL_FILE = EC_CORE_MODEL

    def setUp(self):
        """Set up
        Loads a model
        """
        from cobra.io.sbml import read_sbml_model

        model = read_sbml_model(EC_CORE_MODEL)
        from mewpy.simulation import get_simulator

        self.simul = get_simulator(model)
        k = list(self.simul.objective.keys())
        self.BIOMASS_ID = k[0]
        self.SUCC = "EX_succ_e"

    def test_transport_reactions(self):
        """Tests transport reactions"""
        super().test_transport_reactions()

    def test_room(self):
        """Tests ROOM"""
        self.check_room()


class TestGERM(TestReframedSimul):
    """Tests GERM Simulator"""

    MODEL_FILE = EC_CORE_MODEL2

    def setUp(self):
        """Set up
        Loads a model
        """
        from mewpy.io import read_sbml

        model = read_sbml(EC_CORE_MODEL2, regulatory=False)
        from mewpy.simulation import get_simulator

        self.simul = get_simulator(model)
        k = list(self.simul.objective.keys())
        self.BIOMASS_ID = k[0]
        self.SUCC = "EX_succ_e"

    def test_moma(self):
        """MOMA is not available for GERM models: it raises instead of returning another method's result"""
        with self.assertRaises(NotImplementedError):
            self.simul.simulate(method="MOMA")

    def test_lmoma(self):
        """lMOMA is not available for GERM models"""
        with self.assertRaises(NotImplementedError):
            self.simul.simulate(method="lMOMA")

    def test_room(self):
        """ROOM is not available for GERM models"""
        with self.assertRaises(NotImplementedError):
            self.simul.simulate(method="ROOM")

    def test_FVA(self):
        """FVA of the first two reactions matches COBRApy, in the documented {reaction: [min, max]} format"""
        from cobra.flux_analysis import flux_variability_analysis

        reactions = self.simul.reactions[0:2]
        res = self.simul.FVA(reactions=reactions)
        expected = flux_variability_analysis(
            _cobra_model(self.MODEL_FILE), reaction_list=reactions, fraction_of_optimum=0.9
        )
        self.assertEqual(set(res), set(reactions))
        for r_id in reactions:
            low, high = res[r_id]
            self.assertAlmostEqual(low, expected.loc[r_id, "minimum"], places=6)
            self.assertAlmostEqual(high, expected.loc[r_id, "maximum"], places=6)

    def test_envelope(self):
        """Tests the flux envelope"""
        super().test_envelope()


class TestGeckoLoad(unittest.TestCase):
    """Tests GECKO simulator"""

    def test_gecko(self):
        """A single-pool GECKO model: every enzyme is drawn from the common protein pool"""
        from mewpy.model.gecko import GeckoModel

        model = GeckoModel("single-pool")
        self.assertEqual(model.individual_proteins, frozenset())
        self.assertEqual(model.pool_proteins, model.proteins)
        self.assertEqual(len(model.pool_protein_exchanges), len(model.proteins))
        pool = model.protein_pool_exchange
        ((pool_metabolite, coefficient),) = pool.stoichiometry.items()
        self.assertEqual(coefficient, 1)
        for r_id in model.pool_protein_exchanges:
            stoichiometry = model.reactions[r_id].stoichiometry
            self.assertLess(stoichiometry[pool_metabolite], 0, r_id)
            self.assertEqual(sorted(v for m, v in stoichiometry.items() if m != pool_metabolite), [1], r_id)

    def test_simulator(self):
        """The GECKO simulator's FBA is REFRAMED's FBA and keeps protein usage within the pool"""
        from reframed import FBA

        from mewpy.model.gecko import GeckoModel

        model = GeckoModel("single-pool")
        from mewpy.simulation import SStatus, get_simulator

        simulator = get_simulator(model)
        res = simulator.simulate()
        self.assertEqual(res.status, SStatus.OPTIMAL)
        # the same solver set-up (with SCIP, a zero threshold below the model's ~3e-10 coefficients)
        self.assertAlmostEqual(res.objective_value, FBA(model, solver=simulator._solver_instance()).fobj, places=6)
        pool = model.protein_pool_exchange
        self.assertLessEqual(res.fluxes[pool.id], pool.ub + 1e-9)


class TestGeckoSimul(unittest.TestCase):

    def setUp(self):
        from mewpy.model.gecko import GeckoModel

        model = GeckoModel("single-pool")
        from mewpy.simulation import get_simulator

        self.simul = get_simulator(model)

    def test_essential_proteins(self):
        """
        P33312 is essential; Q12184, without wild-type flux, cannot be (the wild-type optimum stays feasible).
        Two proteins only: all 816 take about 15 minutes.
        """
        from mewpy.model.gecko import GeckoModel

        model = self.simul.model
        draw = {
            p: next(r for r in model.pool_protein_exchanges if r.endswith(p)) for p in (GECKO_ESSENTIAL, GECKO_UNUSED)
        }
        self.assertEqual(self.simul.simulate().fluxes[draw[GECKO_UNUSED]], 0)
        proteins = frozenset(draw)
        with mock.patch.object(GeckoModel, "proteins", new_callable=mock.PropertyMock, return_value=proteins):
            essential = self.simul.essential_proteins()
        self.assertEqual(set(essential), {GECKO_ESSENTIAL})  # protein ids, as on COBRApy


if __name__ == "__main__":
    unittest.main()
