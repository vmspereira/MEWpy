import unittest

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"


class TestCommReframed(unittest.TestCase):

    def setUp(self):
        """Set up - Uses COBRApy models which work properly with community model construction"""
        from cobra.io.sbml import read_sbml_model

        from mewpy.model import CommunityModel

        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "m1"

        model2 = model1.copy()
        model2.id = "m2"

        model3 = model1.copy()
        model3.id = "m3"
        self.models = [model1, model2, model3]
        self.comm = CommunityModel(self.models)

    def test_FBA(self):
        sim = self.comm.get_community_model()
        res = sim.simulate()
        self.assertGreater(res.objective_value, 0)

    def test_SteadyCom(self):
        """
        SteadyCom requires change_coefficients support.
        Currently supported by: CPLEX, Gurobi, PySCIPOpt
        Not supported by: OptLang (GLPK)
        """
        from mewpy.com.steadycom import SteadyCom
        from mewpy.solvers import get_default_solver

        solver_name = get_default_solver()
        if solver_name == "optlang":
            self.skipTest("SteadyCom requires change_coefficients support (not available in OptLang/GLPK)")

        result = SteadyCom(self.comm)
        self.assertIsNotNone(result)
        self.assertGreater(result.growth, 0)

    def test_SteadyComVA(self):
        """
        SteadyComVA requires change_coefficients support.
        Currently supported by: CPLEX, Gurobi, PySCIPOpt
        Not supported by: OptLang (GLPK)
        """
        from mewpy.com.steadycom import SteadyComVA
        from mewpy.solvers import get_default_solver

        solver_name = get_default_solver()
        if solver_name == "optlang":
            self.skipTest("SteadyComVA requires change_coefficients support (not available in OptLang/GLPK)")

        result = SteadyComVA(self.comm)
        self.assertIsNotNone(result)
        self.assertGreater(len(result), 0)


class TestCommCobra(TestCommReframed):

    def setUp(self):
        """Set up"""
        from cobra.io.sbml import read_sbml_model

        from mewpy.model import CommunityModel

        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "m1"

        model2 = model1.copy()
        model2.id = "m2"

        model3 = model1.copy()
        model3.id = "m3"
        self.models = [model1, model2, model3]
        self.comm = CommunityModel(self.models)
