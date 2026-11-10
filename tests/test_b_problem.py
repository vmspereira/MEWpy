import ast
import math
import random
import re
import unittest
from functools import lru_cache
from pathlib import Path

import pytest

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"

OPTRAM_MODEL = MODELS_PATH + "yeast_7.6-optram.xml"
OPTRAM_GENES = MODELS_PATH + "mgene.csv"
OPTRAM_TFS = MODELS_PATH + "TFnames.csv"
OPTRAM_REGNET = MODELS_PATH + "regnet.csv"

EC_CORE_MODEL2 = Path(__file__).parent.joinpath("data", "e_coli_core.xml")
EC_CORE_REG_MODEL = Path(__file__).parent.joinpath("data", "e_coli_core_trn.csv")

MIN_GROWTH = 0.1

MAX_TRIES = 3


def _strip(identifier):
    """The identifier without its SBML prefix (REFRAMED R_/M_/G_), as COBRApy names it."""
    return re.sub(r"^[RMG]_", "", identifier)


@lru_cache(maxsize=None)
def _cobra_model(path):
    """COBRApy reads the same SBML file: the external oracle. Changes are made in `with model:` blocks."""
    from cobra.io import read_sbml_model

    return read_sbml_model(str(path))


@lru_cache(maxsize=None)
def _cobra_essential(path):
    """Essential reactions and genes (COBRApy: knockout growth below 1% of the wild type)."""
    from cobra.flux_analysis import find_essential_genes, find_essential_reactions

    model = _cobra_model(path)
    reactions = {r.id for r in find_essential_reactions(model, processes=1)}
    genes = {g.id for g in find_essential_genes(model, processes=1)}
    return reactions, genes


def _gpr_level(node, levels):
    """Level of a COBRApy GPR tree: min over AND, max over OR, 1 for genes without a level."""
    if isinstance(node, ast.Name):
        return levels.get(node.id, 1)
    values = [_gpr_level(v, levels) for v in node.values]
    return min(values) if isinstance(node.op, ast.And) else max(values)


def _ou_bounds(level, wild_type):
    """Over/under expression bounds: the bound on the side of the wild-type flux becomes level * flux."""
    from mewpy.util.constants import ModelConstants

    if level == 0:
        return (0, 0)
    if level > 1:
        return (
            (level * wild_type, ModelConstants.REACTION_UPPER_BOUND)
            if wild_type > 0
            else (
                -ModelConstants.REACTION_UPPER_BOUND,
                level * wild_type,
            )
        )
    return (0, level * wild_type) if wild_type > 0 else (level * wild_type, 0)


class TestRKOP(unittest.TestCase):

    MODEL_FILE = EC_CORE_MODEL

    def setUp(self):
        from reframed.io.sbml import load_cbmodel

        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.problems import RKOProblem

        self.problem = RKOProblem(model, [])

    @property
    def is_ou(self):
        from mewpy.problems.problem import AbstractOUProblem

        return isinstance(self.problem, AbstractOUProblem)

    def candidates(self, seed, n):
        """n decoded candidates from a seeded generator"""
        rng = random.Random(seed)
        return [self.problem.decode(self.problem.generator(rng)) for _ in range(n)]

    def test_targets(self):
        """Targets exclude essential genes/reactions (COBRApy) and exchanges, and include every other
        enzyme or reaction of a single compartment (transport reactions may be excluded)"""
        from mewpy.problems import GKOProblem, RKOProblem, ROUProblem

        model = _cobra_model(self.MODEL_FILE)
        target = {_strip(t) for t in self.problem.target_list}
        self.assertEqual(len(target), len(self.problem.target_list))
        essential_reactions, essential_genes = _cobra_essential(self.MODEL_FILE)
        internal = [r for r in model.reactions if not r.boundary and len(r.compartments) == 1]
        if isinstance(self.problem, RKOProblem):
            self.assertLessEqual(target, {r.id for r in model.reactions if not r.boundary} - essential_reactions)
            self.assertLessEqual({r.id for r in internal} - essential_reactions, target)
        elif isinstance(self.problem, ROUProblem):
            self.assertEqual(target, {r.id for r in model.reactions})
        else:
            excluded = essential_genes if isinstance(self.problem, GKOProblem) else set()
            self.assertLessEqual(target, {g.id for g in model.genes} - excluded)
            self.assertLessEqual({g.id for r in internal for g in r.genes} - excluded, target)

    def test_generator(self):
        """A candidate has between the minimum and maximum size of distinct targets (and valid levels)"""
        rng = random.Random(0)
        n_targets = len(self.problem.target_list)
        max_size = min(self.problem.candidate_max_size, n_targets)
        for _ in range(20):
            candidate = self.problem.generator(rng)
            self.assertGreaterEqual(len(candidate), self.problem.candidate_min_size)
            self.assertLessEqual(len(candidate), max_size)
            if self.is_ou:
                indexes = [idx for idx, _ in candidate]
                self.assertLessEqual({lv for _, lv in candidate}, set(range(len(self.problem.levels))))
            else:
                indexes = list(candidate)
            self.assertEqual(len(set(indexes)), len(indexes))
            self.assertLessEqual(set(indexes), set(range(n_targets)))

    def test_decode(self):
        import random

        candidate = self.problem.generator(random)
        solution = self.problem.decode(candidate)
        n_candidate = self.problem.encode(solution)
        self.assertEqual(candidate, n_candidate)

    def expected_constraints(self, solution):
        """The constraints a decoded candidate must give, computed independently with COBRApy GPRs"""
        from mewpy.problems import GKOProblem, GOUProblem, RKOProblem, ROUProblem

        model = _cobra_model(self.MODEL_FILE)
        native = {_strip(r): r for r in self.problem.simulator.reactions}
        if isinstance(self.problem, RKOProblem):
            return {r: 0 for r in solution}
        if isinstance(self.problem, GKOProblem):
            knockouts = {_strip(g) for g in solution}
            return {native[r.id]: 0 for r in model.reactions if not r.gpr.eval(knockouts=knockouts)}
        if isinstance(self.problem, ROUProblem):
            levels = {_strip(r): lv for r, lv in solution.items()}
        elif isinstance(self.problem, GOUProblem):
            genes = {_strip(g): lv for g, lv in solution.items()}
            levels = {r.id: _gpr_level(r.gpr.body, genes) for r in model.reactions if r.gpr.body is not None}
        else:
            raise NotImplementedError(type(self.problem).__name__)
        reference = self.problem.reference
        expected = {}
        for r_id, lv in levels.items():
            wild_type = reference[native[r_id]]
            if lv == 0 or (lv != 1 and wild_type != 0):
                expected[native[r_id]] = _ou_bounds(lv, wild_type)
        return expected

    def test_to_constraints(self):
        """Constraints: knocked-out reactions (or those whose GPR becomes false) are closed; over/under
        expression sets the bound on the side of the wild-type flux to level * flux"""
        non_empty = 0
        for solution in self.candidates(seed=1, n=5):
            constraints = self.problem.solution_to_constraints(solution)
            expected = self.expected_constraints(solution)
            self.assertEqual(set(constraints), set(expected), solution)
            for r_id, bounds in expected.items():
                value = constraints[r_id]
                value = value if isinstance(value, tuple) else (value, value)
                bounds = bounds if isinstance(bounds, tuple) else (bounds, bounds)
                self.assertAlmostEqual(value[0], bounds[0], places=9, msg=r_id)
                self.assertAlmostEqual(value[1], bounds[1], places=9, msg=r_id)
            non_empty += bool(constraints)
        self.assertGreater(non_empty, 0)

    def test_simul_constraints(self):
        """Simulating a candidate gives COBRApy's FBA optimum under the same bounds, and respects them"""
        from mewpy.simulation import SStatus

        model = _cobra_model(self.MODEL_FILE)
        for solution in self.candidates(seed=2, n=3):
            constraints = self.problem.solution_to_constraints(solution)
            res = self.problem.simulator.simulate(constraints=constraints)
            with model:
                for r_id, bounds in constraints.items():
                    model.reactions.get_by_id(_strip(r_id)).bounds = (
                        bounds if isinstance(bounds, tuple) else (bounds,) * 2
                    )
                expected = model.slim_optimize(error_value=float("nan"))
            if math.isnan(expected):
                self.assertNotEqual(res.status, SStatus.OPTIMAL, solution)
                continue
            self.assertEqual(res.status, SStatus.OPTIMAL, solution)
            self.assertAlmostEqual(res.objective_value, expected, places=6, msg=solution)
            for r_id, bounds in constraints.items():
                lb, ub = bounds if isinstance(bounds, tuple) else (bounds, bounds)
                self.assertTrue(lb - 1e-6 <= res.fluxes[r_id] <= ub + 1e-6, r_id)


class TestROUP(TestRKOP):

    def setUp(self):
        from reframed.io.sbml import load_cbmodel

        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.problems import ROUProblem

        self.problem = ROUProblem(model, [])


class TestGKOP(TestRKOP):

    def setUp(self):
        from reframed.io.sbml import load_cbmodel

        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.problems import GKOProblem

        self.problem = GKOProblem(model, [])


class TestGOUP(TestRKOP):

    def setUp(self):
        from reframed.io.sbml import load_cbmodel

        model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.problems import GOUProblem

        self.problem = GOUProblem(model, [])


class TestOptRAM(TestRKOP):

    def setUp(self):
        from mewpy.problems import OptRamProblem, load_optram

        regnet = load_optram(OPTRAM_GENES, OPTRAM_TFS, OPTRAM_REGNET, gene_prefix="G_")
        from reframed.io.sbml import load_cbmodel

        model = load_cbmodel(OPTRAM_MODEL)
        self.regnet = regnet
        self.problem = OptRamProblem(model, [], regnet)

    def test_targets(self):
        """OptRAM targets are the regulatory network's metabolic genes, which are model genes, and its TFs"""
        target = self.problem.target_list
        self.assertEqual(len(set(target)), len(target))
        self.assertEqual(set(target), set(self.regnet.genes) | set(self.regnet.tfs))
        self.assertLessEqual(set(self.regnet.genes), set(self.problem.simulator.genes))

    @pytest.mark.xfail
    def test_decode(self):
        import random

        candidate = self.problem.generator(random)
        solution = self.problem.decode(candidate)
        n_candidate = self.problem.encode(solution)
        self.assertEqual(candidate, n_candidate)

    def test_simul_constraints(self):
        """Simulating a candidate gives REFRAMED's FBA optimum under the same bounds, and respects them"""
        from reframed import FBA

        from mewpy.simulation import SStatus

        for solution in self.candidates(seed=2, n=2):
            constraints = self.problem.solution_to_constraints(solution)
            res = self.problem.simulator.simulate(constraints=constraints)
            expected = FBA(self.problem.simulator.model, constraints=constraints)
            self.assertEqual(res.status, SStatus.OPTIMAL, solution)
            self.assertAlmostEqual(res.objective_value, expected.fobj, places=6, msg=solution)
            for r_id, (lb, ub) in constraints.items():
                self.assertTrue(lb - 1e-6 <= res.fluxes[r_id] <= ub + 1e-6, r_id)

    def test_to_constraints(self):
        """OptRAM constraints are bounds (lb <= ub) of model reactions catalysed by the candidate's genes"""
        reactions = set(self.problem.simulator.reactions)
        non_empty = 0
        for solution in self.candidates(seed=1, n=3):
            constraints = self.problem.solution_to_constraints(solution)
            self.assertLessEqual(set(constraints), reactions)
            for r_id, (lb, ub) in constraints.items():
                self.assertLessEqual(lb, ub, r_id)
            non_empty += bool(constraints)
        self.assertGreater(non_empty, 0)


class TestOptORF(unittest.TestCase):

    def setUp(self):

        _BIOMASS_ID = "Biomass_Ecoli_core"
        _O2 = "EX_o2_e"  # noqa: F841
        _GLC = "EX_glc__D_e"  # noqa: F841
        _FUM = "EX_fum_e"  # noqa: F841
        _AC = "EX_ac_e"  # noqa: F841
        _GLU = "EX_glu__L_e"  # noqa: F841
        _LAC = "EX_lac__D_e"  # noqa: F841
        _SUC = "EX_succ_e"  # noqa: F841

        from mewpy.io import Engines, Reader, read_model

        metabolic_reader = Reader(Engines.MetabolicSBML, EC_CORE_MODEL2)
        regulatory_reader = Reader(
            Engines.BooleanRegulatoryCSV, EC_CORE_REG_MODEL, sep=",", id_col=0, rule_col=2, aliases_cols=[1], header=0
        )

        model = read_model(metabolic_reader, regulatory_reader)

        envcond = {"EX_glc__D_e": (-10, 100000.0)}
        from mewpy.simulation import get_simulator

        sim = get_simulator(model, envcond=envcond)
        sim.objective = _BIOMASS_ID

        _PRODUCT_ID = "EX_succ_e"
        from mewpy.optimization.evaluation import BPCY, WYIELD

        evaluator_1 = BPCY(_BIOMASS_ID, _PRODUCT_ID, method="pFBA")
        evaluator_2 = WYIELD(_BIOMASS_ID, _PRODUCT_ID)
        from mewpy.problems import OptORFProblem

        self.model = model
        self.problem = OptORFProblem(model, [evaluator_1, evaluator_2], candidate_max_size=6)

    def test_targets(self):
        """OptORF targets are the regulators, targets and genes that are neither reactions nor metabolites"""
        model = self.model
        expected = {
            v.id
            for v in (*model.yield_regulators(), *model.yield_targets(), *model.yield_genes())
            if not v.is_reaction() and not v.is_metabolite()
        }
        target = self.problem.target_list
        self.assertEqual(len(set(target)), len(target))
        self.assertEqual(set(target), expected)
        self.assertLessEqual(set(model.genes), set(target))

    def test_generator(self):
        """A candidate has 1 to 6 (candidate_max_size) distinct target indexes"""
        rng = random.Random(0)
        n_targets = len(self.problem.target_list)
        for _ in range(20):
            candidate = self.problem.generator(rng)
            self.assertGreaterEqual(len(candidate), self.problem.candidate_min_size)
            self.assertLessEqual(len(candidate), 6)
            self.assertLessEqual(set(candidate), set(range(n_targets)))


if __name__ == "__main__":
    unittest.main()
