"""Batch titer, yield and productivity objectives (dynamic FBA of each candidate)"""

import pickle
from pathlib import Path

import pytest

EC_CORE = str(Path(__file__).parent.joinpath("data", "e_coli_core.xml"))
ANAEROBIC = {"EX_o2_e": (0, 1000)}
SUCC, GLC = "EX_succ_e", "EX_glc__D_e"
# no ethanol or lactate: anaerobic e_coli_core secretes succinate
SUCCINATE_KO = {"ACALD": 0, "LDH_D": 0}


def _model(backend):
    if backend == "cobra":
        from cobra.io import read_sbml_model

        return read_sbml_model(EC_CORE)
    from reframed import load_cbmodel

    return load_cbmodel(EC_CORE)


def _culture():
    """DySScO's batch: 0.01 gDW/L, 20 mM glucose, vmax 20 mmol/gDW/h, Km 1 mM, anaerobic"""
    from mewpy.optimization.evaluation import BatchCulture
    from mewpy.simulation import MichaelisMenten

    return BatchCulture(0.01, {GLC: 20.0, SUCC: 0.0}, uptake={GLC: MichaelisMenten(20, 1)}, t_end=24)


def _problem(backend, culture):
    from mewpy.optimization.evaluation import BatchProductivity, BatchTiter, BatchYield
    from mewpy.problems import RKOProblem

    objectives = [BatchTiter(culture, SUCC), BatchYield(culture, SUCC, GLC), BatchProductivity(culture, SUCC)]
    return RKOProblem(_model(backend), objectives, envcond=ANAEROBIC, candidate_max_size=3)


@pytest.mark.parametrize("backend", ["cobra", "reframed"])
def test_objectives_from_the_candidate_culture(backend):
    """Titer, yield and productivity are those of the DynamicFBA of the candidate (its knockouts as constraints):
    final concentration, product per substrate consumed, product per batch time"""
    from mewpy.simulation import get_simulator

    culture = _culture()
    problem = _problem(backend, culture)
    titer, yld, productivity = problem.evaluate_solution(SUCCINATE_KO, decode=False)

    sim = get_simulator(_model(backend), envcond=ANAEROBIC)
    batch = culture.simulate(sim, SUCCINATE_KO)
    succ, glc = batch.concentrations[sim.resolve(SUCC)], batch.concentrations[sim.resolve(GLC)]
    assert titer == pytest.approx(succ[-1], rel=1e-9)
    assert yld == pytest.approx(succ[-1] / (glc[0] - glc[-1]), rel=1e-9)
    assert productivity == pytest.approx(succ[-1] / 24, rel=1e-9)
    assert titer > 10  # the wild type secretes no succinate (next test)
    assert 0 < yld < 1.5  # at most 6 glucose carbons -> 4 succinate carbons


def test_the_wild_type_makes_no_succinate():
    titer, yld, productivity = _problem("cobra", _culture()).evaluate_solution({}, decode=False)
    assert titer == pytest.approx(0, abs=1e-9) and yld == pytest.approx(0, abs=1e-9)


def test_objectives_that_share_a_culture_simulate_each_candidate_once(monkeypatch):
    from mewpy.simulation.dfba import DynamicFBA

    calls = []
    simulate = DynamicFBA.simulate
    monkeypatch.setattr(DynamicFBA, "simulate", lambda self, *a, **k: calls.append(1) or simulate(self, *a, **k))
    problem = _problem("cobra", _culture())
    first = problem.evaluate_solution(SUCCINATE_KO, decode=False)
    assert len(calls) == 1
    wild_type = problem.evaluate_solution({}, decode=False)
    assert len(calls) == 2 and wild_type != first
    assert problem.evaluate_solution(SUCCINATE_KO, decode=False) == first


def test_an_infeasible_candidate_gets_the_worst_fitness():
    """Without phosphoglycerate kinase there is no anaerobic growth: no culture is simulated"""
    problem = _problem("cobra", _culture())
    worst = [f.worst_fitness for f in problem.fevaluation]
    assert problem.evaluate_solution({"PGK": 0}, decode=False) == worst


def test_untracked_compounds_and_no_consumption():
    from mewpy.optimization.evaluation import BatchTiter, BatchYield
    from mewpy.simulation import SimulationMethod, get_simulator

    culture = _culture()
    sim = get_simulator(_model("cobra"), envcond=ANAEROBIC)
    results = {SimulationMethod.pFBA: sim.simulate(method=SimulationMethod.pFBA)}
    with pytest.raises(ValueError, match="not a tracked compound"):
        BatchTiter(culture, "EX_ac_e")(results, None)
    # succinate is not consumed: no yield on it
    no_yield = BatchYield(culture, GLC, SUCC)
    assert no_yield(results, None) == no_yield.no_solution


def test_a_pickled_culture_drops_its_caches():
    culture = _culture()
    problem = _problem("cobra", culture)
    expected = problem.evaluate_solution(SUCCINATE_KO, decode=False)
    assert "_last" in culture.__dict__ and "_sim_cache" in culture.__dict__
    copy = pickle.loads(pickle.dumps(problem))
    shared = {id(f.culture) for f in copy.fevaluation}
    assert len(shared) == 1  # the objectives still share one culture
    assert "_last" not in copy.fevaluation[0].culture.__dict__
    assert copy.evaluate_solution(SUCCINATE_KO, decode=False) == pytest.approx(expected, rel=1e-9)


def test_an_ea_runs_on_the_batch_objectives():
    from mewpy.optimization import EA

    problem = _problem("cobra", _culture())
    ea = EA(problem, max_generations=1, population_size=4, mp=False)
    solutions = ea.run(simplify=False)
    assert solutions and all(len(s.fitness) == 3 for s in solutions)
