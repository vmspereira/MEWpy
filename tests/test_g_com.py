import math
import unittest
from functools import lru_cache

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"
# SteadyCom finds the growth rate by binary search, to this relative tolerance
STEADYCOM_REL_TOL = 1e-4


@lru_cache(maxsize=None)
def open_medium_growth():
    """Growth of one member alone in the community medium (COBRApy): ATPM off and every exchange unbounded,
    as the community model leaves its exchanges open. Identical members share it: the community grows as one."""
    from cobra.io.sbml import read_sbml_model

    model = read_sbml_model(EC_CORE_MODEL)
    model.reactions.get_by_id("ATPM").bounds = (0, 0)
    for rxn in model.exchanges:
        rxn.bounds = (-math.inf, math.inf)
    return model.slim_optimize()


def check_community_fba(test, community, res):
    """Community FBA of identical members: the growth of one member alone, every member growing at the community
    rate, and each shared metabolite's exchange equal to the sum of the members' exchanges (mass balance)"""
    from mewpy.simulation import SStatus

    test.assertEqual(res.status, SStatus.OPTIMAL)
    test.assertAlmostEqual(res.objective_value, open_medium_growth(), delta=1e-6 * open_medium_growth())
    for org_id, biomass in community.organisms_biomass.items():
        test.assertAlmostEqual(res.fluxes[biomass], res.objective_value, delta=1e-6 * res.objective_value, msg=org_id)
    sim = community.get_community_model()
    balanced = 0
    for r_id in sim.get_exchange_reactions():
        members = [f"{r_id}_{org_id}" for org_id in community.organisms]
        if all(m in res.fluxes for m in members):
            test.assertAlmostEqual(res.fluxes[r_id], sum(res.fluxes[m] for m in members), delta=1e-6, msg=r_id)
            balanced += 1
    test.assertEqual(balanced, 20)


def check_steadycom(test, community, result):
    """SteadyCom of identical members: the growth of one member alone (to the binary search tolerance) and
    abundances that sum to 1"""
    growth = open_medium_growth()
    test.assertLessEqual(result.growth, growth + 1e-6)
    test.assertGreaterEqual(result.growth, growth * (1 - STEADYCOM_REL_TOL))
    test.assertEqual(set(result.abundance), set(community.organisms))
    test.assertAlmostEqual(sum(result.abundance.values()), 1.0, places=6)
    for org_id, x in result.abundance.items():
        test.assertGreaterEqual(x, -1e-9, org_id)
        test.assertLessEqual(x, 1 + 1e-9, org_id)


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
        """Community FBA of three identical members"""
        sim = self.comm.get_community_model()
        res = sim.simulate()
        check_community_fba(self, self.comm, res)

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
        check_steadycom(self, self.comm, result)

    def test_SteadyComVA(self):
        """
        SteadyComVA requires change_coefficients support.
        Currently supported by: CPLEX, Gurobi, PySCIPOpt
        Not supported by: OptLang (GLPK)

        Identical members: at the maximal growth any split is optimal, so each abundance ranges over [0, 1].
        """
        from mewpy.com.steadycom import SteadyComVA
        from mewpy.solvers import get_default_solver

        solver_name = get_default_solver()
        if solver_name == "optlang":
            self.skipTest("SteadyComVA requires change_coefficients support (not available in OptLang/GLPK)")

        result = SteadyComVA(self.comm)
        self.assertEqual(set(result), set(self.comm.organisms))
        for org_id, (low, high) in result.items():
            self.assertAlmostEqual(low, 0.0, places=6, msg=org_id)
            self.assertAlmostEqual(high, 1.0, places=6, msg=org_id)


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
