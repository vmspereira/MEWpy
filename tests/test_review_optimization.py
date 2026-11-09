"""Regression tests for optimization/problems defects found in the 2026-11-09 review.

Each test names the review finding it guards (F1..F12, L1..L4).
"""

import contextlib
import logging
import pickle
import random
import signal
import subprocess
import sys
import unittest
import warnings
from pathlib import Path
from types import SimpleNamespace as NS
from unittest import mock

DATA = Path(__file__).parent.joinpath("data")
EC_CORE_MODEL = str(DATA.joinpath("e_coli_core.xml.gz"))
BIOMASS_ID = "R_BIOMASS_Ecoli_core_w_GAM"
SUCC_ID = "R_EX_succ_e"
SMALL_TARGETS = ["R_PFL", "R_ACKr", "R_LDH_D"]


@contextlib.contextmanager
def time_limit(seconds):
    """Raises TimeoutError if the block runs longer than `seconds` (POSIX only)."""
    if not hasattr(signal, "SIGALRM"):
        yield
        return

    def _handler(signum, frame):
        raise TimeoutError(f"still running after {seconds} s")

    previous = signal.signal(signal.SIGALRM, _handler)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def load_ecoli():
    from reframed.io.sbml import load_cbmodel

    return load_cbmodel(EC_CORE_MODEL)


class Doubler:
    """Picklable evaluable used by the multiprocessing evaluator tests."""

    def evaluate(self, candidate):
        return candidate * 2

    def evaluator(self, candidates, *args):
        return [self.evaluate(c) for c in candidates]


# Exercises the 'mp' and 'nodaemon' pool evaluators in a fresh interpreter
POOL_SCRIPT = """
import multiprocessing, time
from mewpy.util import process


class Doubler:
    def evaluate(self, candidate):
        return candidate * 2


class Failing:
    def evaluate(self, candidate):
        raise ValueError("evaluation failed")


if __name__ == "__main__":
    for kind in ("mp", "nodaemon"):
        print(kind, process.get_evaluator(Doubler(), n_mp=2, evaluator=kind)([1, 2, 3], None))
        try:
            process.get_evaluator(Failing(), n_mp=2, evaluator=kind)([1, 2, 3], None)
        except ValueError as e:
            print(kind, "raised", e)
        deadline = time.time() + 5
        while multiprocessing.active_children() and time.time() < deadline:
            time.sleep(0.1)
        print(kind, "children", len(multiprocessing.active_children()))
"""


class TestEngines(unittest.TestCase):
    def test_f1_jmetal_mo_with_single_objective_algorithm_name(self):
        from mewpy.optimization.evaluation import BPCY, CandidateSize
        from mewpy.optimization.jmetal.ea import EA as JMetalEA
        from mewpy.problems import RKOProblem

        fevaluation = [BPCY(BIOMASS_ID, SUCC_ID), CandidateSize()]
        problem = RKOProblem(load_ecoli(), fevaluation, target=SMALL_TARGETS)
        ea = JMetalEA(problem, max_generations=1, mp=False, algorithm="GA", population_size=4)
        self.assertGreater(len(ea.run(simplify=False)), 0)

    def _run_twice(self, engine, seed):
        from mewpy.optimization import EA, get_default_engine, set_default_engine
        from mewpy.optimization.evaluation import BPCY
        from mewpy.problems import RKOProblem

        previous = get_default_engine()
        set_default_engine(engine)
        try:
            problem = RKOProblem(load_ecoli(), [BPCY(BIOMASS_ID, SUCC_ID)], candidate_max_size=3)
            runs = []
            kwargs = dict(max_generations=2, mp=False, algorithm="GA", population_size=6, seed=seed)
            for _ in range(2):
                ea = EA(problem, **kwargs)
                runs.append(sorted(str(s) for s in ea.run(simplify=False)))
            return runs
        finally:
            set_default_engine(previous)

    def test_f11_seed_makes_inspyred_reproducible(self):
        first, second = self._run_twice("inspyred", seed=7)
        self.assertEqual(first, second)

    def test_f11_seed_makes_jmetal_reproducible(self):
        first, second = self._run_twice("jmetal", seed=7)
        self.assertEqual(first, second)


class TestEvaluationFunctions(unittest.TestCase):
    def test_f2_minimised_objective_penalises_lethal_candidate(self):
        from mewpy.optimization.evaluation import CandidateSize, TargetFlux
        from mewpy.problems import RKOProblem

        f = TargetFlux("R_EX_co2_e", biomass=BIOMASS_ID, maximize=False, min_biomass_value=0.1)
        problem = RKOProblem(load_ecoli(), [f], target=["R_PFL", "R_GAPD"])
        problem.pre_process()
        viable = problem.evaluate_solution({0})[0]  # R_PFL KO grows
        lethal = problem.evaluate_solution({1})[0]  # R_GAPD KO does not grow
        self.assertGreater(lethal, viable)
        self.assertGreater(CandidateSize(maximize=False).no_solution, 1000)
        self.assertEqual(CandidateSize(maximize=True).no_solution, 0.0)

    def test_f5_bpcy_fva_does_not_mutate_simulation_constraints(self):
        from mewpy.optimization.evaluation import BPCY_FVA, CandidateSize
        from mewpy.problems import RKOProblem

        fevaluation = [BPCY_FVA(BIOMASS_ID, SUCC_ID), CandidateSize()]
        problem = RKOProblem(load_ecoli(), fevaluation, target=SMALL_TARGETS)
        problem.pre_process()
        self.assertEqual(problem.evaluate_solution({0})[1], 1)

        constraints = {"R_PFL": 0}
        result = problem.simulator.simulate(constraints=constraints, method="pFBA")
        result.simulation_constraints = constraints
        BPCY_FVA(BIOMASS_ID, SUCC_ID).get_fitness({"pFBA": result}, {"R_PFL": 0})
        self.assertEqual(constraints, {"R_PFL": 0})

    def _fva_problem(self):
        from mewpy.optimization.evaluation import BPCY_FVA, WYIELD
        from mewpy.problems import RKOProblem

        fevaluation = [WYIELD(BIOMASS_ID, SUCC_ID), BPCY_FVA(BIOMASS_ID, SUCC_ID)]
        envcond = {"R_EX_o2_e": (0, 0)}
        problem = RKOProblem(load_ecoli(), fevaluation, target=SMALL_TARGETS, envcond=envcond)
        problem.pre_process()
        return problem

    def test_f_sim_12_fva_functions_reuse_one_simulator_copy(self):
        import mewpy.optimization.evaluation.phenotype as phenotype

        problem = self._fva_problem()
        with mock.patch.object(phenotype, "get_simulator", wraps=phenotype.get_simulator) as spy:
            for candidate in ({0}, {1}, {0, 2}):
                problem.evaluate_solution(candidate)
        self.assertEqual(spy.call_count, 2)  # one copy per evaluation function
        restored = pickle.loads(pickle.dumps(problem.fevaluation))
        self.assertTrue(all("_sim_cache" not in f.__dict__ for f in restored))

    def test_f_sim_12_reused_simulator_copy_gives_fresh_values(self):
        reused = self._fva_problem()
        values = [reused.evaluate_solution(c) for c in ({0}, {1}, {0, 2}, {1})]
        fresh = [self._fva_problem().evaluate_solution(c) for c in ({0}, {1}, {0, 2}, {1})]
        self.assertEqual(values, fresh)

    def test_f_sim_problem_simulate_does_not_mutate_caller_constraints(self):
        from mewpy.problems import RKOProblem

        problem = RKOProblem(load_ecoli(), [], target=SMALL_TARGETS)
        user_constraints = {"R_EX_glc__D_e": (-5, 0)}
        problem.simulate(solution={"R_PFL": 0}, constraints=user_constraints)
        problem.FVA(reactions=["R_PFL"], solution={"R_ACKr": 0}, constraints=user_constraints)
        self.assertEqual(user_constraints, {"R_EX_glc__D_e": (-5, 0)})

    def test_f6_target_flux_unknown_reaction_raises(self):
        from mewpy.optimization.evaluation import TargetFlux
        from mewpy.simulation import get_simulator

        result = get_simulator(load_ecoli()).simulate(method="pFBA")
        with self.assertRaises(ValueError):
            TargetFlux("R_EX_succ_e_TYPO").get_fitness({"pFBA": result}, {})

    def test_f6_target_flux_without_reaction_uses_model_objective(self):
        from mewpy.optimization.evaluation import TargetFlux
        from mewpy.simulation import get_simulator

        result = get_simulator(load_ecoli()).simulate(method="pFBA")
        value = TargetFlux(None).get_fitness({"pFBA": result}, {})
        self.assertAlmostEqual(value, result.fluxes[BIOMASS_ID], places=6)

    def test_f7_evaluation_failure_is_logged_once(self):
        from mewpy.optimization.evaluation import BPCY
        from mewpy.problems import RKOProblem

        problem = RKOProblem(load_ecoli(), [BPCY("R_BIOMASS_TYPO", SUCC_ID)], target=SMALL_TARGETS)
        problem.pre_process()
        with self.assertLogs("mewpy.problems.problem", level="WARNING") as logs:
            self.assertEqual(problem.evaluate_solution({0}), [0.0])
        self.assertIn("R_BIOMASS_TYPO", "\n".join(logs.output))
        with self.assertNoLogs("mewpy.problems.problem", level="WARNING"):
            problem.evaluate_solution({1})

    def test_l2_min_cand_size_constructs(self):
        from mewpy.optimization.evaluation import MinCandSize

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            f = MinCandSize()
        self.assertFalse(f.maximize)
        self.assertTrue(any(issubclass(w.category, DeprecationWarning) for w in caught))


class TestCandidateSize(unittest.TestCase):
    def test_f3_ko_generator_with_target_list_smaller_than_max_size(self):
        from mewpy.problems import RKOProblem

        # default candidate_max_size is 10
        problem = RKOProblem(load_ecoli(), [], target=SMALL_TARGETS)
        for i in range(30):
            candidate = problem.generator(random.Random(i))
            self.assertTrue(1 <= len(candidate) <= len(SMALL_TARGETS))

    def test_f3_community_ko_smoke(self):
        from reframed.io.sbml import load_cbmodel

        from mewpy.optimization.evaluation import TargetFlux
        from mewpy.problems.com import CommunityKOProblem

        m1, m2 = load_cbmodel(EC_CORE_MODEL), load_cbmodel(EC_CORE_MODEL)
        m1.id, m2.id = "ec1", "ec2"
        problem = CommunityKOProblem([m1, m2], [TargetFlux("R_EX_ac_e")])
        problem.pre_process()
        for i in range(10):
            candidate = problem.generator(random.Random(i))
            self.assertTrue(1 <= len(candidate) <= 2)
        fitness = problem.evaluate_solution({0})
        self.assertEqual(len(fitness), 1)

    def test_f3_grow_mutations_terminate_on_full_candidate(self):
        from mewpy.optimization.inspyred.operators import grow_mutation_KO, grow_mutation_OU
        from mewpy.optimization.jmetal.operators import GrowMutationKO, GrowMutationOU
        from mewpy.optimization.jmetal.problem import KOSolution, OUSolution
        from mewpy.problems.problem import KOBounder, OUBounder

        ko_args = {"_ec": NS(bounder=KOBounder(0, 4)), "gs_mutation_rate": 1.0}
        ko_args.update(candidate_min_size=1, candidate_max_size=10)
        ou_args = dict(ko_args, _ec=NS(bounder=OUBounder([0, 0], [4, 3])))
        full_ou = {(i, 1) for i in range(5)}
        with time_limit(3):
            full_ko = set(range(5))
            self.assertEqual(grow_mutation_KO(random.Random(0), [set(full_ko)], ko_args), [full_ko])
            self.assertEqual(grow_mutation_OU(random.Random(0), [set(full_ou)], ou_args), [full_ou])
            s = KOSolution(0, 4, 5, 1)
            s.variables = [0, 1, 2, 3, 4]
            self.assertEqual(GrowMutationKO(1.0, max_size=10).execute(s).variables, [0, 1, 2, 3, 4])
            s = OUSolution([0, 0], [4, 3], 5, 1)
            s.variables = sorted(full_ou)
            self.assertEqual(GrowMutationOU(1.0, max_size=10).execute(s).variables, sorted(full_ou))

    def test_f3_ou_ea_with_small_target_list_terminates(self):
        from mewpy.optimization.evaluation import TargetFlux
        from mewpy.optimization.inspyred.ea import EA as InspyredEA
        from mewpy.problems import ROUProblem

        problem = ROUProblem(load_ecoli(), [TargetFlux(SUCC_ID)], target=SMALL_TARGETS)
        problem.operators_param = {"gs_mutation_rate": 1.0}  # always try to grow
        with time_limit(30):
            ea = InspyredEA(problem, max_generations=3, mp=False, population_size=6, seed=1)
            self.assertGreater(len(ea.run(simplify=False)), 0)

    def test_l4_max_candidate_size_alias_is_honoured(self):
        from mewpy.problems import GKOProblem, ROUProblem, RKOProblem

        model = load_ecoli()
        for cls in (RKOProblem, ROUProblem, GKOProblem):
            with self.assertWarns(DeprecationWarning):
                problem = cls(model, [], max_candidate_size=4, min_candidate_size=2)
            sizes = (problem.candidate_min_size, problem.candidate_max_size)
            self.assertEqual(sizes, (2, 4), cls.__name__)
        with self.assertRaises(ValueError):
            RKOProblem(model, [], max_candidate_size=4, candidate_max_size=5)


class TestEncodeDecode(unittest.TestCase):
    def _optram(self, fevaluation):
        from reframed.io.sbml import load_cbmodel

        from mewpy.problems import OptRamProblem, load_optram

        files = [str(DATA / f) for f in ("mgene.csv", "TFnames.csv", "regnet.csv")]
        regnet = load_optram(*files, gene_prefix="G_")
        model = load_cbmodel(str(DATA / "yeast_7.6-optram.xml"))
        problem = OptRamProblem(model, fevaluation, regnet)
        return problem, len(regnet.genes)

    def test_f4_optram_simplify_after_conversion(self):
        from mewpy.optimization.evaluation import CandidateSize
        from mewpy.optimization.inspyred.ea import EA as InspyredEA

        problem, n_genes = self._optram([CandidateSize()])
        ea = InspyredEA(problem, mp=False)
        individual = NS(candidate={(n_genes, 6)}, fitness=1.0)  # a TF over-expressed 2x
        population = ea._convertPopulation([individual])
        simplified = problem.simplify_population(population)
        self.assertEqual(len(simplified), 1)
        self.assertEqual(simplified[0].values, population[0].values)

    def test_f4_optram_simplify_without_encoded_candidate_keeps_solution(self):
        from mewpy.optimization.ea import Solution
        from mewpy.optimization.evaluation import CandidateSize

        problem, n_genes = self._optram([CandidateSize()])
        values = problem.decode({(n_genes, 6)})
        solution = Solution(values, [len(values)])
        with self.assertLogs("mewpy.problems.problem", level="WARNING"):
            self.assertEqual(problem.simplify_population([solution]), [solution])

    def test_f4_optorf_decode_identifies_knockouts(self):
        from mewpy.io import Engines, Reader, read_model
        from mewpy.optimization.evaluation import CandidateSize
        from mewpy.problems import OptORFProblem

        metabolic = Reader(Engines.MetabolicSBML, str(DATA / "e_coli_core.xml"))
        regulatory = Reader(
            Engines.BooleanRegulatoryCSV,
            str(DATA / "e_coli_core_trn.csv"),
            sep=",",
            id_col=0,
            rule_col=2,
            aliases_cols=[1],
            header=0,
        )
        problem = OptORFProblem(read_model(metabolic, regulatory), [CandidateSize()])
        problem.pre_process()
        candidate = {problem.target_list.index("surplusFDP")}
        decoded = problem.decode(candidate)
        self.assertEqual(decoded, {"surplusFDP": 0})
        self.assertEqual(problem.encode(decoded), candidate)
        self.assertEqual(problem.evaluate_solution(candidate), [1])
        constraints = problem.solution_to_constraints(decoded)
        self.assertIn("FORt2", constraints)
        self.assertTrue(all(v == (0.0, 0.0) for v in constraints.values()))

    def test_f4_medium_problem_smoke(self):
        from mewpy.optimization.evaluation import TargetFlux
        from mewpy.problems.reactions import MediumProblem

        problem = MediumProblem(load_ecoli(), [TargetFlux(BIOMASS_ID)])
        problem.pre_process()
        candidate = problem.generator(random.Random(0))
        self.assertEqual(problem.encode(problem.decode(candidate)), candidate)
        self.assertEqual(len(problem.evaluate_solution(candidate)), 1)

    def test_f9_hybrid_gecko_applies_all_enzyme_constraints(self):
        from mewpy.problems import GeckoOUProblem
        from mewpy.problems.hybrid import HybridGeckoOUProblem

        p = HybridGeckoOUProblem.__new__(HybridGeckoOUProblem)
        p.vmaxs, p.prot_prefix, p.initcond = [], "draw_prot_", None
        p.apply_lb, p.lb_tolerance, p.gDW = False, 0.05, 564.0
        p.ksim = NS(simulate=lambda factors, initcon: NS(fluxes={"K1": 1.0, "K2": 2.0}))
        p.kmodel = NS(merge_constants=lambda: {"vmax1": 10.0, "vmax2": 20.0})
        p.enzyme_mapping = {
            "K1": NS(sense=1, forward={"PA": 100.0}, backward={}, vmax_id="vmax1"),
            "K2": NS(sense=1, forward={"PB": 100.0}, backward={}, vmax_id="vmax2"),
        }
        with mock.patch.object(GeckoOUProblem, "solution_to_constraints", lambda self, c: {}):
            constraints = p.solution_to_constraints({})
        self.assertEqual(set(constraints), {"draw_prot_PA", "draw_prot_PB"})

    def test_f9_hybrid_gecko_target_list_with_non_target(self):
        from mewpy.problems.hybrid import HybridGeckoOUProblem

        p = HybridGeckoOUProblem.__new__(HybridGeckoOUProblem)
        p.kmodel = NS(get_parameters=lambda exclude_compartments: ["vmax_A", "vmax_B", "km_A"])
        p.enzyme_mapping = NS(proteins=["PA"])
        p._simul = NS(proteins=["PA", "PB", "PC"])
        p.non_target, p._partial_solution = ["vmax_B", "PC"], {}
        p._build_target_list()
        self.assertEqual(sorted(p._trg_list), ["PB", "vmax_A"])


class TestModelSafety(unittest.TestCase):
    def test_f10_kcatopt_does_not_mutate_model(self):
        from mewpy.model.gecko import GeckoModel
        from mewpy.optimization.evaluation import TargetFlux
        from mewpy.problems.gecko import KcatOptProblem

        model = GeckoModel("single-pool")
        prot, rxn = "P38790", "r_2148No1"
        met = [m for m in model.reactions[rxn].stoichiometry if prot in m][0]
        before = model.reactions[rxn].stoichiometry[met]
        problem = KcatOptProblem(model, [prot], [TargetFlux(rxn)])
        for _ in range(2):  # kcat edits must not build up across candidates
            sim = problem.solution_to_constraints({(prot, rxn): 2.0})
            self.assertAlmostEqual(sim.model.reactions[rxn].stoichiometry[met], before / 2)
            self.assertEqual(model.reactions[rxn].stoichiometry[met], before)
            model._needs_update = True  # as after any model edit: refreshes the kcat lookup


class TestOperators(unittest.TestCase):
    def test_f8_jmetal_ou_crossover_keeps_unique_targets(self):
        from mewpy.optimization.jmetal.operators import UniformCrossoverOU
        from mewpy.optimization.jmetal.problem import OUSolution

        random.seed(0)
        for _ in range(500):
            mom = OUSolution([0, 0], [9, 10], 3, 1)
            mom.variables = [(1, 1), (2, 1), (3, 1)]
            dad = OUSolution([0, 0], [9, 10], 2, 1)
            dad.variables = [(1, 2), (4, 2)]
            for child in UniformCrossoverOU(1.0, 3).execute([mom, dad]):
                idx = [i for i, _ in child.variables]
                self.assertEqual(len(idx), len(set(idx)), child.variables)
                self.assertLessEqual(len(idx), 3)

    def test_f11_ou_mutators_use_the_given_prng(self):
        from mewpy.optimization.inspyred.operators import (
            single_mutation_OU,
            single_mutation_OU_level,
        )
        from mewpy.problems.problem import OUBounder

        args = {"_ec": NS(bounder=OUBounder([0, 0], [9, 10])), "mutation_rate": 1.0}
        for op in (single_mutation_OU, single_mutation_OU_level):
            results = []
            for global_seed in (1, 2, 3):
                random.seed(global_seed)
                results.append(op(random.Random(5), [{(1, 2), (3, 4)}], args))
            self.assertEqual(results[0], results[1], op.__name__)
            self.assertEqual(results[1], results[2], op.__name__)

    def test_f11_target_lists_are_sorted(self):
        from mewpy.problems import GKOProblem, ROUProblem, RKOProblem

        model = load_ecoli()
        for cls in (RKOProblem, ROUProblem, GKOProblem):
            targets = cls(model, []).target_list
            self.assertEqual(targets, sorted(targets), cls.__name__)

    def test_l1_jmetal_solution_equality(self):
        from mewpy.optimization.jmetal.problem import KOSolution, OUSolution

        a, b, c = KOSolution(0, 10, 2, 1), KOSolution(0, 10, 2, 1), KOSolution(0, 10, 2, 1)
        a.variables, b.variables, c.variables = [3, 1], [7, 9], [1, 3]
        self.assertNotEqual(a, b)
        self.assertEqual(a, c)
        self.assertEqual(a.variables, [3, 1])
        x, y = OUSolution([0, 0], [9, 9], 1, 1), OUSolution([0, 0], [9, 9], 1, 1)
        x.variables, y.variables = [(1, 2)], [(1, 3)]
        self.assertNotEqual(x, y)

    def test_l3_dominance_direction(self):
        from mewpy.optimization.ea import Solution, cmetric, dominance_test
        from mewpy.optimization.jmetal.problem import KOSolution

        self.assertEqual(dominance_test([2.0], [1.0], maximize=True), 1)
        self.assertEqual(dominance_test([2.0], [1.0], maximize=False), -1)
        a, b = Solution({"R1": 0}, [2.0, 2.0]), Solution({"R2": 0}, [1.0, 1.0])
        self.assertTrue(a > b)
        self.assertTrue(b < a)
        r1, r2, _, _ = cmetric([a], [b], maximize=True)
        self.assertEqual((r1, r2), (1.0, 0.0))
        # jMetal minimises: objectives [1] dominate [2]
        s1, s2 = KOSolution(0, 9, 1, 1), KOSolution(0, 9, 1, 1)
        s1.objectives, s2.objectives = [1.0], [2.0]
        self.assertTrue(s1 > s2)

    def test_l3_non_dominated_population_unchanged(self):
        from mewpy.optimization.ea import Solution, non_dominated_population

        a, b = Solution({"R1": 0}, [2.0, 2.0]), Solution({"R2": 0}, [1.0, 1.0])
        c = Solution({"R3": 0}, [3.0, 0.5])
        self.assertEqual(non_dominated_population([a, b, c], maximize=True), [a, c])
        self.assertEqual(non_dominated_population([a, b, c], maximize=False), [b, c])


class TestProcess(unittest.TestCase):
    def test_f12_get_evaluator_reads_mp_evaluator_at_call_time(self):
        from mewpy.util import process
        from mewpy.util.constants import ModelConstants

        previous = ModelConstants.MP_EVALUATOR
        ModelConstants.MP_EVALUATOR = "mp"
        try:
            with mock.patch.object(process, "RayEvaluator", lambda *a, **k: "ray", create=True):
                evaluator = process.get_evaluator(Doubler(), n_mp=2)
        finally:
            ModelConstants.MP_EVALUATOR = previous
        self.assertIsInstance(evaluator, process.MultiProcessorEvaluator)

    def test_f12_pool_evaluators(self):
        # run in a fresh interpreter: forking this process, which may run Ray threads, is unsafe
        cmd = [sys.executable, "-c", POOL_SCRIPT]
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        self.assertEqual(out.returncode, 0, out.stderr[-2000:])
        lines = out.stdout.splitlines()
        for kind in ("mp", "nodaemon"):
            self.assertIn(f"{kind} [2, 4, 6]", lines)
            self.assertIn(f"{kind} raised evaluation failed", lines)
            self.assertIn(f"{kind} children 0", lines)

    def test_f12_cpu_count_at_least_one(self):
        from mewpy.util import process
        from mewpy.util.constants import EAConstants

        with mock.patch.object(EAConstants, "NUM_CPUS", -1), mock.patch.object(
            process.multiprocessing, "cpu_count", return_value=1
        ):
            self.assertEqual(process.cpu_count(), 1)

    def test_f12_import_does_not_set_start_method(self):
        code = "import multiprocessing as m, mewpy.util.process"
        code += "; print(m.get_start_method(allow_none=True))"
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stderr[-2000:])
        self.assertEqual(out.stdout.strip().splitlines()[-1], "None")


if __name__ == "__main__":
    logging.basicConfig()
    unittest.main()
