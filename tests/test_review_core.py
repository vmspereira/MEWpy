"""Regression tests for blocking defects found in the 2026-11-09 review.

F-PKG-1: mewpy.model imported the undeclared 'six' package.
F-PKG-2: jmetalpy >= 1.10 calls problem.create_solution(rng).
F-SIM-1: the COBRA simulator assigned mewpy's solver name to an optlang model.
F-SIM-2: get_simulator hid simulator construction errors.
"""

import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock

EC_CORE_MODEL = str(Path(__file__).parent.joinpath("data", "e_coli_core.xml.gz"))
BIOMASS_ID = "R_BIOMASS_Ecoli_core_w_GAM"
SUCC_ID = "R_EX_succ_e"


class TestReviewCoreRegressions(unittest.TestCase):
    def test_package_does_not_import_six(self):
        # A scan of the source, because blocking six at run time also blocks dependencies that
        # import it (python-dateutil < 4, installed on Python <= 3.10).
        import ast

        import mewpy

        found = []
        for path in sorted(Path(mewpy.__file__).parent.rglob("*.py")):
            for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"), str(path))):
                if isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.level == 0:
                    names = [node.module]
                else:
                    continue
                if any(name == "six" or name.startswith("six.") for name in names):
                    found.append("%s:%d" % (path, node.lineno))
        self.assertEqual(found, [])

    def test_cobra_simulator_with_solver_unknown_to_optlang(self):
        from cobra.io import read_sbml_model

        import mewpy.solvers as solvers
        from mewpy.simulation import get_simulator

        previous = solvers.get_default_solver()
        # "optlang" is a mewpy solver name but not an optlang interface name
        solvers.set_default_solver("optlang")
        try:
            simulator = get_simulator(read_sbml_model(EC_CORE_MODEL))
            self.assertGreater(simulator.simulate().objective_value, 0.8)
        finally:
            solvers.set_default_solver(previous)

    def test_get_simulator_reports_construction_error(self):
        from cobra.io import read_sbml_model

        from mewpy.simulation import get_simulator

        model = read_sbml_model(EC_CORE_MODEL)
        target = "mewpy.simulation.cobra.Simulation.__init__"
        with mock.patch(target, side_effect=RuntimeError("boom")):
            with self.assertRaises(Exception) as ctx:
                get_simulator(model)
        self.assertIn("boom", f"{ctx.exception} {ctx.exception.__cause__}")

    def test_environment_apply_inplace_on_simulator(self):
        """Environment.apply(simulator, inplace=True) must change that simulator, not a copy."""
        from reframed.io.sbml import load_cbmodel

        from mewpy.simulation import get_simulator
        from mewpy.simulation.environment import Environment

        simulator = get_simulator(load_cbmodel(EC_CORE_MODEL))
        environment = Environment()
        environment["R_EX_glc__D_e"] = (-1.0, 1000.0)
        environment.apply(simulator, exclusive=False, inplace=True)
        self.assertEqual(simulator.get_reaction_bounds("R_EX_glc__D_e")[0], -1.0)

    def test_kinetic_worker_uses_fork_when_available(self):
        """mewpy no longer forces a global 'fork'; the kinetic worker must still use fork."""
        import multiprocessing

        from mewpy.simulation.kinetic import KineticThread

        if "fork" not in multiprocessing.get_all_start_methods():
            self.skipTest("fork start method not available")
        self.assertEqual(KineticThread._start_method, "fork")

    def test_cobra_fva_in_unguarded_spawn_script(self):
        """COBRApy FVA must not open a process pool by default: an unguarded script hangs under spawn."""
        import tempfile

        script = (
            "import multiprocessing; multiprocessing.set_start_method('spawn', force=True)\n"
            "from cobra.io import read_sbml_model\n"
            "from mewpy.simulation import get_simulator\n"
            f"sim = get_simulator(read_sbml_model({EC_CORE_MODEL!r}))\n"
            # >= 2 reactions: cobra uses min(processes, n_reactions) workers, so one reaction never opens a pool
            "print(sim.FVA(reactions=['PGI', 'ACKr'], obj_frac=0.9))\n"
        )
        # a real file: spawn re-imports __main__ from its file (python -c would not reproduce the hang)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "unguarded_fva.py"
            path.write_text(script)
            try:
                result = subprocess.run([sys.executable, str(path)], capture_output=True, text=True, timeout=120)
            except subprocess.TimeoutExpired:
                self.fail("FVA hung in an unguarded script using the spawn start method")
        self.assertEqual(result.returncode, 0, result.stderr[-2000:])
        self.assertIn("PGI", result.stdout)

    def test_jmetal_engine_runs(self):
        from reframed.io.sbml import load_cbmodel

        from mewpy.optimization import EA, get_default_engine, set_default_engine
        from mewpy.optimization.evaluation import BPCY
        from mewpy.optimization.settings import set_default_population_size
        from mewpy.problems import RKOProblem

        previous = get_default_engine()
        set_default_engine("jmetal")
        set_default_population_size(4)
        try:
            model = load_cbmodel(EC_CORE_MODEL)
            problem = RKOProblem(model, [BPCY(BIOMASS_ID, SUCC_ID)], candidate_max_size=3)
            ea = EA(problem, max_generations=1, mp=False, algorithm="GA")
            self.assertGreater(len(ea.run()), 0)
        finally:
            set_default_engine(previous)


if __name__ == "__main__":
    unittest.main()
