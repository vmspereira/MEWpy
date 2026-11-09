import unittest
from pathlib import Path

# the GERM problems run the checks of the REFRAMED ones (test_b_problem), with GERM models and COBRApy oracles.
# A module import: pytest collects every TestCase class in this module's namespace, so the classes are not imported.
import test_b_problem as b

EC_CORE_MODEL = Path(__file__).parent.joinpath("data", "e_coli_core.xml")

MIN_GROWTH = 0.1

MAX_TRIES = 3


class TestRKOP(b.TestRKOP):

    MODEL_FILE = EC_CORE_MODEL

    def setUp(self):
        from mewpy.io import read_sbml

        model = read_sbml(EC_CORE_MODEL, regulatory=False, warnings=False)
        from mewpy.problems import RKOProblem

        self.problem = RKOProblem(model, [])


class TestROUP(b.TestROUP):

    MODEL_FILE = EC_CORE_MODEL

    def setUp(self):
        from mewpy.io import read_sbml

        model = read_sbml(EC_CORE_MODEL, regulatory=False, warnings=False)
        from mewpy.problems import ROUProblem

        self.problem = ROUProblem(model, [])


class TestGKOP(b.TestGKOP):

    MODEL_FILE = EC_CORE_MODEL

    def setUp(self):
        from mewpy.io import read_sbml

        model = read_sbml(EC_CORE_MODEL, regulatory=False, warnings=False)
        from mewpy.problems import GKOProblem

        self.problem = GKOProblem(model, [])

    def test_targets(self):
        """Tests the target list"""
        super().test_targets()


class TestGOUP(b.TestGOUP):

    MODEL_FILE = EC_CORE_MODEL

    def setUp(self):
        from mewpy.io import read_sbml

        model = read_sbml(EC_CORE_MODEL, regulatory=False, warnings=False)
        from mewpy.problems import GOUProblem

        self.problem = GOUProblem(model, [])

    def test_targets(self):
        """Tests the target list"""
        super().test_targets()


if __name__ == "__main__":
    unittest.main()
