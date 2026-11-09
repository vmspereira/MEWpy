import itertools
import math
import unittest
from functools import lru_cache

import pytest

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"
BIOMASS = "BIOMASS_Ecoli_core_w_GAM"
# eight genes, each of one or more e_coli_core reactions, with distinct expression levels
GENES8 = ["b1779", "b4025", "b2926", "b1702", "b0720", "b2914", "b1478", "b3403"]
LEVELS8 = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]


@lru_cache(maxsize=None)
def _cobra_model():
    """COBRApy model of the same SBML file: the external oracle. Changes are made in `with model:` blocks."""
    from cobra.io.sbml import read_sbml_model

    return read_sbml_model(EC_CORE_MODEL)


def _gimme_optimum(coefficients, growth_frac=0.9):
    """GIMME LP in COBRApy: minimise sum(c_r * |v_r|) with growth >= growth_frac * the wild-type growth"""
    model = _cobra_model()
    with model:
        model.reactions.get_by_id(BIOMASS).lower_bound = growth_frac * model.slim_optimize()
        expression = sum(
            c * (model.reactions.get_by_id(r).forward_variable + model.reactions.get_by_id(r).reverse_variable)
            for r, c in coefficients.items()
        )
        model.objective = model.problem.Objective(expression, direction="min")
        return model.slim_optimize()


def _imat_optimum(low, high, epsilon=1):
    """Most agreements (iMAT): highly expressed reactions with |v| >= epsilon plus lowly expressed ones with
    |v| <= epsilon, over every feasible choice (brute force, COBRApy)"""
    model = _cobra_model()
    choices = [(None, "fwd", "rev")] * len(high) + [(None, "off")] * len(low)
    best = 0
    for combo in itertools.product(*choices):
        count = sum(c is not None for c in combo)
        if count <= best:
            continue
        with model:
            for r_id, choice in zip(list(high) + list(low), combo):
                rxn = model.reactions.get_by_id(r_id)
                if choice == "fwd":
                    rxn.lower_bound = max(rxn.lower_bound, epsilon)
                elif choice == "rev":
                    rxn.upper_bound = min(rxn.upper_bound, -epsilon)
                elif choice == "off":
                    rxn.bounds = (max(rxn.lower_bound, -epsilon), min(rxn.upper_bound, epsilon))
            if not math.isnan(model.slim_optimize(error_value=float("nan"))):
                best = count
    return best


class TestExpressionSet(unittest.TestCase):

    def setUp(self):
        import numpy as np

        self.genes = ["b0356", "b1478", "b3734", "b3731"]
        self.conditions = ["Exp#1", "Exp#2", "Exp#3"]
        self.expression = np.array(([0.17, 0.20, 0.93], [0.36, 0.83, 0.77], [0.87, 0.65, 0.07], [0.55, 0.49, 0.52]))
        from mewpy.omics import ExpressionSet

        self.expr = ExpressionSet(self.genes, self.conditions, self.expression)
        self.expr8 = ExpressionSet(GENES8, ["c1"], np.array([[v] for v in LEVELS8]))
        from cobra.io.sbml import read_sbml_model

        model = read_sbml_model(EC_CORE_MODEL)
        from mewpy.simulation import get_simulator

        self.sim = get_simulator(model)

    def reaction_expression(self, expr):
        from mewpy.omics import Preprocessing

        return Preprocessing(self.sim, expr).reactions_expression(0)

    def test_GIMME(self):
        """GIMME minimises the use of lowly expressed reactions keeping 90% growth: its objective is the optimum
        of that LP (COBRApy) and the weighted use of the returned fluxes"""
        from mewpy.omics import GIMME, Preprocessing
        from mewpy.simulation import SStatus

        wild_type = _cobra_model().slim_optimize()
        for expr in (self.expr, self.expr8):
            solution = GIMME(self.sim, expr, cutoff=100)
            coefficients, _ = Preprocessing(self.sim, expr).percentile(0, cutoff=100)
            self.assertEqual(solution.status, SStatus.OPTIMAL)
            self.assertGreaterEqual(solution.fluxes[BIOMASS], 0.9 * wild_type - 1e-6)
            expected = _gimme_optimum(coefficients)
            self.assertAlmostEqual(solution.objective_value, expected, delta=1e-6)
            used = sum(c * abs(solution.fluxes[r]) for r, c in coefficients.items())
            self.assertAlmostEqual(used, solution.objective_value, delta=1e-6)
        # GENES8: penalised reactions must carry flux at 90% growth, so the optimum is not trivial
        self.assertGreater(expected, 1.0)

    def test_GIMME_build(self):
        """The GIMME context-specific model keeps the flux of the GIMME solution, so it grows at least 90%"""
        from mewpy.omics import GIMME

        wild_type = _cobra_model().slim_optimize()
        solution, sim = GIMME(self.sim, self.expr, build_model=True)
        growth = sim.simulate().objective_value
        self.assertGreaterEqual(solution.fluxes[BIOMASS], 0.9 * wild_type - 1e-6)
        self.assertGreaterEqual(growth, 0.9 * wild_type - 1e-6)
        self.assertLessEqual(growth, wild_type + 1e-6)

    def test_eFlux(self):
        """eFlux bounds every reaction by its expression relative to the maximum (1 without data): the
        objective is FBA with those bounds (COBRApy), and the built model reproduces it"""
        from mewpy.omics import eFlux

        solution = eFlux(self.sim, self.expr)
        rxn_exp = self.reaction_expression(self.expr)
        max_exp = max(rxn_exp.values())
        model = _cobra_model()
        with model:
            for rxn in model.reactions:
                level = rxn_exp.get(rxn.id, max_exp) / max_exp
                rxn.bounds = (-level if rxn.lower_bound < 0 else 0, level if rxn.upper_bound > 0 else 0)
            expected = model.slim_optimize()
        self.assertGreater(expected, 0.01)
        self.assertAlmostEqual(solution.objective_value, expected, delta=1e-6)
        built_solution, built = eFlux(self.sim, self.expr, build_model=True)
        self.assertAlmostEqual(built.simulate().objective_value, built_solution.objective_value, delta=1e-6)

    def test_iMAT(self):
        """iMAT maximises agreements: reactions above the 75th expression percentile carrying |v| >= 1 plus
        reactions below the 25th carrying |v| <= 1; its objective is the most agreements (brute force)"""
        import numpy as np

        from mewpy.omics import iMAT
        from mewpy.simulation import SStatus

        for expr in (self.expr, self.expr8):
            solution = iMAT(self.sim, expr)
            rxn_exp = self.reaction_expression(expr)
            low_cut, high_cut = np.percentile(list(rxn_exp.values()), [25, 75])
            low = [r for r, v in rxn_exp.items() if v < low_cut]
            high = [r for r, v in rxn_exp.items() if v > high_cut]
            self.assertEqual(solution.status, SStatus.OPTIMAL)
            self.assertLessEqual(solution.objective_value, len(low) + len(high) + 1e-6)
            self.assertAlmostEqual(solution.objective_value, _imat_optimum(low, high), delta=1e-6)


if __name__ == "__main__":
    test = TestExpressionSet()
    test.setUp()
    test.test_GIMME()
