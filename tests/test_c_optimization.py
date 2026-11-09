import importlib.util
import math
import re
import unittest
from functools import lru_cache
from pathlib import Path

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"
EC_CORE_MODEL2 = Path(__file__).parent.joinpath("data", "e_coli_core.xml")
BIOMASS_ID = "R_BIOMASS_Ecoli_core_w_GAM"
SUCC_ID = "R_EX_succ_e"
# e_coli_core.xml (GERM tests) names its biomass reaction differently
GERM_BIOMASS_ID = "Biomass_Ecoli_core"
GERM_SUCC_ID = "EX_succ_e"
MIN_GROWTH = 0.1
POPULATION_SIZE = 10
CANDIDATE_MAX_SIZE = 6


def _strip(identifier):
    """The identifier without its SBML prefix (REFRAMED R_/M_/G_), as COBRApy names it."""
    return re.sub(r"^[RMG]_", "", identifier)


@lru_cache(maxsize=None)
def _cobra_model(path):
    """COBRApy reads the same SBML file: the external oracle. Changes are made in `with model:` blocks."""
    from cobra.io import read_sbml_model

    return read_sbml_model(str(path))


def _native(simulator, identifier):
    """The simulator's id of a reaction given with or without the SBML prefix"""
    return identifier if identifier in simulator.reactions else _strip(identifier)


def _cobra_bounds(model, constraints):
    for r_id, bounds in constraints.items():
        model.reactions.get_by_id(_strip(r_id)).bounds = bounds if isinstance(bounds, tuple) else (bounds,) * 2


def _cobra_growth(path, constraints, biomass):
    """Optimal growth under the candidate's constraints (COBRApy), nan if infeasible"""
    model = _cobra_model(path)
    with model:
        _cobra_bounds(model, constraints)
        model.objective = _strip(biomass)
        return model.slim_optimize(error_value=float("nan"))


def _cobra_range(path, constraints, biomass, min_growth, product):
    """Range of the product flux under the candidate's constraints and growth >= min_growth (COBRApy)"""
    model = _cobra_model(path)
    with model:
        _cobra_bounds(model, constraints)
        model.reactions.get_by_id(_strip(biomass)).lower_bound = min_growth
        model.objective = _strip(product)
        high = model.slim_optimize()
        model.objective_direction = "min"
        low = model.slim_optimize()
    return low, high


class TestOptInspyred(unittest.TestCase):
    """Unittests of Inspyred based optimizations."""

    MODEL_FILE = EC_CORE_MODEL
    BIOMASS = BIOMASS_ID

    def setUp(self):
        """Sets up the the model"""
        from reframed.io.sbml import load_cbmodel

        self.model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.optimization.settings import set_default_population_size

        set_default_population_size(10)
        from mewpy.optimization import get_available_engines, set_default_engine

        if len(get_available_engines()):
            set_default_engine("inspyred")

    def check_population(self, ea, problem):
        """The final population: at most the population size of distinct candidates within the size limit,
        whose constraints decode their modifications and whose fitness is recomputed independently.
        WYIELD and the BPCY_FVA range come from COBRApy FBA/FVA; fluxes of other methods from a fresh simulation,
        except FBA fluxes (alternative optima), bounded by their range at the optimal growth; ModificationType
        from the modifications."""
        from mewpy.optimization.evaluation import BPCY, BPCY_FVA, WYIELD, ModificationType, TargetFlux
        from mewpy.simulation import SimulationMethod, SStatus, get_simulator

        population = ea.final_population
        self.assertGreaterEqual(len(population), 1)
        self.assertLessEqual(len(population), POPULATION_SIZE)
        self.assertEqual(len({frozenset(s.values.items()) for s in population}), len(population))
        simulator = get_simulator(problem.simulator.model)
        biomass = _native(simulator, self.BIOMASS)
        wild_type = _cobra_growth(self.MODEL_FILE, {}, biomass)
        for solution in population:
            constraints = solution.constraints
            # simplification may remove every modification (the wild type has the same fitness)
            self.assertLessEqual(len(solution.values), CANDIDATE_MAX_SIZE)
            self.assertEqual(constraints, problem.solution_to_constraints(solution.values))
            self.assertEqual(len(solution.fitness), len(problem.fevaluation))
            if not solution.values and any(isinstance(f, ModificationType) for f in problem.fevaluation):
                # ModificationType divides by the number of modifications: the evaluation fails
                self.assertEqual(solution.fitness, [f.no_solution for f in problem.fevaluation])
                continue
            for f, fitness in zip(problem.fevaluation, solution.fitness):
                msg = f"{type(f).__name__} of {solution.values}"
                if isinstance(f, ModificationType):
                    penalties = [5 if v == 0 else 2 if v < 1 else 0 for v in solution.values.values()]
                    self.assertAlmostEqual(fitness, sum(penalties) / len(penalties), places=9, msg=msg)
                    continue
                product = _native(simulator, f.reaction if isinstance(f, TargetFlux) else f.productId)
                if isinstance(f, WYIELD):
                    growth = _cobra_growth(self.MODEL_FILE, constraints, biomass)
                    if math.isnan(growth) or f.obj_frac * growth <= f.min_biomass_per * wild_type:
                        self.assertEqual(fitness, f.no_solution, msg)
                    else:
                        low, high = _cobra_range(self.MODEL_FILE, constraints, biomass, f.obj_frac * growth, product)
                        self.assertAlmostEqual(fitness, f.alpha * high + (1 - f.alpha) * low, delta=1e-6, msg=msg)
                    continue
                res = simulator.simulate(constraints=constraints, method=f.method)
                if res.status != SStatus.OPTIMAL:
                    self.assertEqual(fitness, f.no_solution, msg)
                    continue
                growth, flux = res.fluxes[biomass], res.fluxes[product]
                if f.method == SimulationMethod.FBA:
                    # any flux of its range at the optimal growth is an FBA flux
                    low, high = _cobra_range(self.MODEL_FILE, constraints, biomass, growth - 1e-9, product)
                    scale = growth if isinstance(f, BPCY) else 1.0
                    self.assertGreaterEqual(fitness, scale * low - 1e-6, msg)
                    self.assertLessEqual(fitness, scale * high + 1e-6, msg)
                elif isinstance(f, BPCY):
                    self.assertAlmostEqual(fitness, growth * flux, delta=1e-6, msg=msg)
                elif isinstance(f, BPCY_FVA):
                    low, high = _cobra_range(self.MODEL_FILE, constraints, biomass, f.obj_frac * growth, product)
                    if abs(high) < 1e-9:
                        expected = f.worst_fitness
                    elif abs(abs(high) - abs(low)) < 1e-9:
                        expected = growth * flux
                    else:
                        expected = growth * flux * (1 - math.log(abs((high - low) / (high + low))))
                    self.assertAlmostEqual(fitness, expected, delta=1e-6 * max(1.0, abs(expected)), msg=msg)
                else:
                    self.assertAlmostEqual(fitness, flux, delta=1e-6, msg=msg)

    def test_engine(self):
        """The available engines are the installed ones: inspyred and jmetalpy"""
        from mewpy.optimization import get_available_engines

        installed = {
            name
            for name, package in (("inspyred", "inspyred"), ("jmetal", "jmetal"))
            if importlib.util.find_spec(package) is not None
        }
        self.assertEqual(set(get_available_engines()), installed)
        self.assertIn("inspyred", installed)

    def test_KOProblem(self):
        """Tests KO problems"""
        from mewpy.optimization.evaluation import BPCY, WYIELD

        f1 = BPCY(BIOMASS_ID, SUCC_ID, method="lMOMA")
        f2 = WYIELD(BIOMASS_ID, SUCC_ID)
        from mewpy.problems import RKOProblem

        problem = RKOProblem(self.model, [f1, f2], candidate_max_size=6)
        from mewpy.optimization import EA

        ea = EA(problem, max_generations=2)
        ea.run()
        self.assertEqual(ea.get_population_size(), 10)
        self.check_population(ea, problem)

    def test_OUProblem(self):
        """Tests OU problems"""
        from mewpy.optimization.evaluation import BPCY_FVA, ModificationType, TargetFlux

        f1 = BPCY_FVA(BIOMASS_ID, SUCC_ID, method="lMOMA")
        f2 = TargetFlux(SUCC_ID)
        f3 = ModificationType()
        from mewpy.problems import ROUProblem

        problem = ROUProblem(self.model, [f1, f2, f3], candidate_max_size=6)
        from mewpy.optimization import EA

        ea = EA(problem, max_generations=1)
        ea.run()
        self.assertEqual(ea.get_population_size(), 10)
        self.check_population(ea, problem)


class TestOptJMetal(TestOptInspyred):
    """Unittests for JMetalPy based optimizations."""

    def setUp(self):
        """Sets up the the model"""
        from mewpy.optimization import get_available_engines

        available = get_available_engines()
        if "jmetal" not in available:
            raise unittest.SkipTest("JMetal optimization engine not available")

        from reframed.io.sbml import load_cbmodel

        self.model = load_cbmodel(EC_CORE_MODEL)
        from mewpy.optimization.settings import set_default_population_size

        set_default_population_size(10)
        from mewpy.optimization import set_default_engine

        set_default_engine("jmetal")


class TestGERMOptInspyred(TestOptInspyred):
    """Unittests for Inspyred based optimizations using germ models."""

    MODEL_FILE = EC_CORE_MODEL2
    BIOMASS = GERM_BIOMASS_ID

    def setUp(self):
        """Sets up the the model"""
        from mewpy.io import read_sbml

        self.model = read_sbml(EC_CORE_MODEL2, regulatory=False, warnings=False)
        from mewpy.optimization.settings import set_default_population_size

        set_default_population_size(10)
        from mewpy.optimization import get_available_engines, set_default_engine

        if len(get_available_engines()):
            set_default_engine("inspyred")

    def test_KOProblem(self):
        """Tests KO problems"""
        from mewpy.optimization.evaluation import BPCY, WYIELD

        # "FBA": the SimulationMethod value ("fba" is no method, so every candidate got the worst fitness)
        f1 = BPCY(GERM_BIOMASS_ID, GERM_SUCC_ID, method="FBA")
        f2 = WYIELD(GERM_BIOMASS_ID, GERM_SUCC_ID, method="FBA")
        from mewpy.problems import RKOProblem

        problem = RKOProblem(self.model, [f1, f2], candidate_max_size=6)
        from mewpy.optimization import EA

        ea = EA(problem, max_generations=2)
        ea.run()
        df = ea.dataframe()
        self.assertEqual(len(df), len(ea.final_population))
        self.assertEqual(list(df.columns), ["Modification", "Size", "BPCY", "WYIELD"])
        self.assertEqual(ea.get_population_size(), 10)
        self.check_population(ea, problem)

    def test_OUProblem(self):
        """Tests OU problems"""
        from mewpy.optimization.evaluation import BPCY_FVA, ModificationType, TargetFlux

        f1 = BPCY_FVA(GERM_BIOMASS_ID, GERM_SUCC_ID)
        f2 = TargetFlux(GERM_SUCC_ID, method="FBA")
        f3 = ModificationType()
        from mewpy.problems import ROUProblem

        problem = ROUProblem(self.model, [f1, f2, f3], candidate_max_size=6)
        from mewpy.optimization import EA

        ea = EA(problem, max_generations=1)
        ea.run()
        self.assertEqual(ea.get_population_size(), 10)
        self.check_population(ea, problem)


class TestGERMOptJMetal(TestGERMOptInspyred):
    """Unittests for JMetalPy based optimizations using germ models."""

    def setUp(self):
        """Sets up the the model"""
        from mewpy.optimization import get_available_engines

        available = get_available_engines()
        if "jmetal" not in available:
            raise unittest.SkipTest("JMetal optimization engine not available")

        from mewpy.io import read_sbml

        self.model = read_sbml(EC_CORE_MODEL2, regulatory=False, warnings=False)
        from mewpy.optimization.settings import set_default_population_size

        set_default_population_size(10)
        from mewpy.optimization import set_default_engine

        set_default_engine("jmetal")


if __name__ == "__main__":
    unittest.main()
