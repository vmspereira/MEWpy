"""Problems and evaluation functions accept ids with or without the SBML prefix (R_/M_/G_).

REFRAMED keeps the prefixes (R_PFL), COBRApy drops them (PFL). The same problem, defined with
either form, must give the same targets and the same fitness on both backends.
"""

import unittest
from pathlib import Path

EC_CORE_MODEL = str(Path(__file__).parent.joinpath("data", "e_coli_core.xml.gz"))

# unprefixed ids: native for COBRApy, need the prefix for REFRAMED
BIOMASS = "BIOMASS_Ecoli_core_w_GAM"
PRODUCT = "EX_succ_e"
ANAEROBIC = {"EX_o2_e": (0, 0)}
REACTIONS = ["PFL", "ACKr", "LDH_D", "ALCD2x", "PTAr"]
GENES = ["b0902", "b3951", "b1241", "b2296", "b1380"]


def load(backend):
    if backend == "reframed":
        from reframed.io.sbml import load_cbmodel

        return load_cbmodel(EC_CORE_MODEL)
    from cobra.io import read_sbml_model

    return read_sbml_model(EC_CORE_MODEL)


def prefixed(ids, prefix):
    return [prefix + i for i in ids]


def fevaluation(biomass=BIOMASS, product=PRODUCT):
    from mewpy.optimization.evaluation import BPCY, WYIELD, TargetFlux

    # FBA, not pFBA: the REFRAMED simulator runs pFBA at 99.9% of the optimum and COBRApy at 100%
    return [
        BPCY(biomass, product, method="FBA"),
        WYIELD(biomass, product),
        TargetFlux(product, biomass=biomass, method="FBA"),
    ]


class TestPrefixAgnosticProblems(unittest.TestCase):
    def _fitness(self, cls, backend, targets, candidate, **kwargs):
        problem = cls(load(backend), fevaluation(**kwargs), target=targets, envcond=ANAEROBIC)
        problem.pre_process()
        return problem.evaluate_solution(problem.encode(candidate))

    def assertFitnessEqual(self, first, second):
        self.assertEqual(len(first), len(second))
        for a, b in zip(first, second):
            self.assertAlmostEqual(a, b, places=6)

    def test_rko_unprefixed_ids_same_fitness_on_both_backends(self):
        from mewpy.problems import RKOProblem

        candidate = ["PFL", "LDH_D"]
        reframed = self._fitness(RKOProblem, "reframed", REACTIONS, candidate)
        cobra = self._fitness(RKOProblem, "cobra", REACTIONS, candidate)
        self.assertGreater(reframed[0], 0.0)  # BPCY is positive under anaerobiosis
        self.assertFitnessEqual(reframed, cobra)

    def test_rko_prefixed_ids_on_cobra(self):
        from mewpy.problems import RKOProblem

        candidate = prefixed(["PFL", "LDH_D"], "R_")
        ids = dict(biomass="R_" + BIOMASS, product="R_" + PRODUCT)
        cobra = self._fitness(RKOProblem, "cobra", prefixed(REACTIONS, "R_"), candidate, **ids)
        reframed = self._fitness(RKOProblem, "reframed", REACTIONS, ["PFL", "LDH_D"])
        self.assertGreater(cobra[0], 0.0)
        self.assertFitnessEqual(cobra, reframed)

    def test_gko_unprefixed_ids_same_fitness_on_both_backends(self):
        from mewpy.problems import GKOProblem

        candidate = ["b0902", "b3951"]  # knocks out PFL
        reframed = self._fitness(GKOProblem, "reframed", GENES, candidate)
        cobra = self._fitness(GKOProblem, "cobra", GENES, candidate)
        self.assertGreater(reframed[0], 0.0)  # BPCY is 0 for the wild type
        self.assertFitnessEqual(reframed, cobra)

    def test_non_target_in_the_other_form_is_excluded(self):
        from mewpy.problems import GKOProblem, RKOProblem

        targets = RKOProblem(load("reframed"), [], non_target=["PFL", "LDH_D"]).target_list
        self.assertTrue({"R_PFL", "R_LDH_D"}.isdisjoint(targets))
        self.assertIn("R_ACKr", targets)
        targets = RKOProblem(load("cobra"), [], non_target=["R_PFL", "R_LDH_D"]).target_list
        self.assertTrue({"PFL", "LDH_D"}.isdisjoint(targets))
        self.assertIn("ACKr", targets)
        targets = GKOProblem(load("reframed"), [], non_target=["b1380"]).target_list
        self.assertNotIn("G_b1380", targets)
        self.assertIn("G_b2133", targets)

    def test_rou_reference_and_partial_solution_in_the_other_form(self):
        from mewpy.problems import ROUProblem
        from mewpy.simulation import get_simulator

        def fitness(prefix):
            model = load("reframed")
            reference = get_simulator(model, envcond=ANAEROBIC).simulate(method="pFBA").fluxes
            reference = {prefix + k[2:]: v for k, v in reference.items()}  # keyed in the given form
            problem = ROUProblem(
                model,
                fevaluation(prefix + BIOMASS, prefix + PRODUCT),
                target=prefixed(REACTIONS, prefix),
                envcond=ANAEROBIC,
                reference=reference,
                partial_solution={prefix + "PFL": 0},
            )
            problem.pre_process()
            return problem.evaluate_solution(problem.encode({prefix + "LDH_D": 0.5}))

        native = fitness("R_")
        other_form = fitness("")
        self.assertGreater(native[0], 0.0)
        self.assertFitnessEqual(native, other_form)

    def test_hybrid_gou_constraints_use_native_ids(self):
        from mewpy.problems.hybrid import HybridGOUProblem

        problem = HybridGOUProblem(load("reframed"), {"PFL": (0.0, 1.0)}, [], target=GENES)
        constraints = problem.solution_to_constraints({})
        self.assertEqual(constraints.get("R_PFL"), (0.0, 1.0))
        self.assertNotIn("PFL", constraints)

    def test_optram_regulatory_genes_without_the_model_gene_prefix(self):
        from reframed.io.sbml import load_cbmodel

        from mewpy.problems import OptRamProblem, load_optram

        data = Path(__file__).parent.joinpath("data")
        files = [str(data / f) for f in ("mgene.csv", "TFnames.csv", "regnet.csv")]
        model = load_cbmodel(str(data / "yeast_7.6-optram.xml"))
        prefixed_problem = OptRamProblem(model, [], load_optram(*files, gene_prefix="G_"))
        plain_problem = OptRamProblem(model, [], load_optram(*files, gene_prefix=""))
        for idx in range(len(prefixed_problem.regmodel.genes)):
            expected = prefixed_problem.solution_to_constraints(prefixed_problem.decode({(idx, 0)}))
            if expected:  # a gene whose under-expression constrains some reaction
                break
        self.assertTrue(expected)
        plain = plain_problem.solution_to_constraints(plain_problem.decode({(idx, 0)}))
        self.assertEqual(plain, expected)

    def test_seed_given_in_the_other_form_is_encoded(self):
        from mewpy.problems import RKOProblem

        problem = RKOProblem(load("reframed"), [], target=prefixed(REACTIONS, "R_"))
        self.assertEqual(problem.encode(["PFL"]), problem.encode(["R_PFL"]))


if __name__ == "__main__":
    unittest.main()
