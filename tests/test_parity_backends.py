"""Backend parity: the same SBML and user code give the same results on COBRApy and REFRAMED.

e_coli_core.xml is loaded with cobra.io.read_sbml_model and reframed.io.sbml.load_cbmodel. Ids
the user passes are tried in both forms ('PGI' and 'R_PGI'); results are compared on canonical
(unprefixed) ids. Where COBRApy has its own implementation (cobra.flux_analysis), it is the
oracle; otherwise the other backend or a documented constant is. Only well-defined quantities
are compared: objective values, FVA bounds, essential sets, knockout growth, and fluxes that FVA
shows are fixed at the optimum.

The helpers derive native and canonical ids themselves, not with Simulator.resolve, so the file
also runs against versions without the id API. ROOM is not covered: after a PGI knockout on
e_coli_core it took 17 s on REFRAMED/SCIP, and on COBRApy/GLPK it had not finished after 12 min.
"""

import csv
import math
import re
import statistics
from copy import deepcopy
from functools import lru_cache
from pathlib import Path

import pytest

from mewpy.simulation import SStatus, get_simulator

DATA = Path(__file__).parent.joinpath("data")
EC_CORE = str(DATA.joinpath("e_coli_core.xml"))
EC_CORE_TRN = str(DATA.joinpath("e_coli_core_trn.csv"))
TRN_KWARGS = dict(sep=",", id_col=0, rule_col=2, aliases_cols=[1], header=0)

BACKENDS = ["cobra", "reframed"]
FORMS = ["unprefixed", "prefixed"]
SBML_PREFIX = {"reaction": "R_", "metabolite": "M_", "gene": "G_"}
OBJ_TOL = 1e-6
FLUX_TOL = 1e-4
MIN_GROWTH = 0.01  # the default threshold of Simulator.essential_genes() / essential_reactions()

BIOMASS = "Biomass_Ecoli_core"
PRODUCT = "EX_succ_e"
UPTAKE = "EX_glc__D_e"
# FRD7 and SUCDi form a loop that only the +-1000 bounds limit; REFRAMED reads those bounds as +-inf
LOOP = {"FRD7", "SUCDi"}
FIXED_AEROBIC = ["EX_glc__D_e", "EX_o2_e", "EX_co2_e", "PGI", "CS", "ATPS4r", BIOMASS]
FIXED_ANAEROBIC = ["EX_ac_e", "EX_etoh_e", "EX_for_e", BIOMASS]
# docs/germ.md:861-862, e_coli_core with e_coli_core_trn.csv
RFBA_DOCS_WT = 0.8739215069684303
SRFBA_DOCS_WT = 0.8739215066796725


# ----------------------------------------------------------------------------------------------------------------------
# Models, ids and oracles
# ----------------------------------------------------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _cobra_core():
    from cobra.io import read_sbml_model

    return read_sbml_model(EC_CORE)


@lru_cache(maxsize=None)
def _reframed_core():
    from reframed.io.sbml import load_cbmodel

    return load_cbmodel(EC_CORE)


def load_model(backend):
    return _cobra_core().copy() if backend == "cobra" else deepcopy(_reframed_core())


def load_sim(backend, **kwargs):
    return get_simulator(load_model(backend), **kwargs)


def oracle_model(anaerobic=False):
    """A COBRApy model solved by cobra itself, not through MEWpy."""
    model = _cobra_core().copy()
    if anaerobic:
        model.reactions.EX_o2_e.lower_bound = 0
    return model


def native(backend, base, kind="reaction"):
    return SBML_PREFIX[kind] + base if backend == "reframed" else base


def given(base, form, kind="reaction"):
    """An id as the user writes it."""
    return SBML_PREFIX[kind] + base if form == "prefixed" else base


def canonical(identifier):
    return identifier[2:] if identifier[:2] in ("R_", "M_", "G_") else identifier


def finite(value, bound=1000.0):
    """Maps REFRAMED's +-inf back to the SBML's +-1000."""
    return max(-bound, min(bound, value))


def no_growth(value):
    return value is None or math.isnan(value) or value < MIN_GROWTH


@lru_cache(maxsize=None)
def oracle_fva(obj_frac, anaerobic=False):
    from cobra.flux_analysis import flux_variability_analysis

    return flux_variability_analysis(oracle_model(anaerobic), fraction_of_optimum=obj_frac, processes=1)


@lru_cache(maxsize=None)
def oracle_pfba(anaerobic=False):
    """cobra's pFBA at the fraction of the optimum that every MEWpy backend keeps (ModelConstants.PFBA_OBJ_FRAC)."""
    from cobra.flux_analysis import pfba

    from mewpy.util.constants import ModelConstants

    return pfba(oracle_model(anaerobic), fraction_of_optimum=ModelConstants.PFBA_OBJ_FRAC)


@lru_cache(maxsize=None)
def oracle_deletions(kind):
    """{gene or reaction: growth, or None if infeasible} from cobra's single deletions."""
    from cobra.flux_analysis import single_gene_deletion, single_reaction_deletion

    run = single_gene_deletion if kind == "gene" else single_reaction_deletion
    table = run(oracle_model(), processes=1)
    rows = zip(table["ids"], table["growth"], table["status"])
    return {next(iter(ids)): (growth if status == "optimal" else None) for ids, growth, status in rows}


@lru_cache(maxsize=None)
def oracle_lmoma(knockout):
    from cobra.flux_analysis import moma

    model = oracle_model()
    model.reactions.get_by_id(knockout).knock_out()
    return moma(model, solution=oracle_pfba(), linear=True).objective_value


def wt_reference(form, anaerobic=False):
    """COBRApy's wild-type pFBA fluxes, keyed by ids in the given form."""
    return {given(r_id, form): value for r_id, value in oracle_pfba(anaerobic).fluxes.items()}


def assert_fixed_fluxes(fluxes, backend, fva, reactions):
    for r_id in reactions:
        low, high = fva.loc[r_id, "minimum"], fva.loc[r_id, "maximum"]
        assert high - low < 1e-6, f"{r_id} is not fixed at the optimum"
        assert fluxes[native(backend, r_id)] == pytest.approx(low, abs=FLUX_TOL), r_id


def assert_growth(value, expected, label):
    if expected is None:
        assert value is None or math.isnan(value), label
    else:
        assert value == pytest.approx(expected, abs=OBJ_TOL), label


@pytest.fixture(params=BACKENDS)
def backend(request):
    return request.param


@pytest.fixture(params=FORMS)
def form(request):
    return request.param


def test_no_canonical_id_looks_prefixed():
    """canonical() strips R_/M_/G_; it is only safe if no COBRApy id starts with one."""
    model = oracle_model()
    for container in (model.reactions, model.metabolites, model.genes):
        assert not [x.id for x in container if canonical(x.id) != x.id]


# ----------------------------------------------------------------------------------------------------------------------
# 1. Simulation
# ----------------------------------------------------------------------------------------------------------------------
def test_fba_objective_and_fixed_fluxes_match_cobra(backend, form):
    result = load_sim(backend).simulate(objective={given(BIOMASS, form): 1})
    assert result.status == SStatus.OPTIMAL
    assert result.objective_value == pytest.approx(oracle_model().slim_optimize(), abs=OBJ_TOL)
    assert_fixed_fluxes(result.fluxes, backend, oracle_fva(1.0), FIXED_AEROBIC)


def _pfba(backend, form):
    return load_sim(backend).simulate(method="pFBA", objective={given(BIOMASS, form): 1})


@pytest.mark.parametrize("backend", ["cobra", "reframed"])
def test_pfba_total_flux_matches_cobra(backend, form):
    result, expected = _pfba(backend, form), oracle_pfba()
    assert result.objective_value == pytest.approx(expected.objective_value, abs=OBJ_TOL)
    assert sum(abs(v) for v in result.fluxes.values()) == pytest.approx(expected.objective_value, abs=FLUX_TOL)
    assert result.fluxes[native(backend, BIOMASS)] == pytest.approx(expected.fluxes[BIOMASS], abs=OBJ_TOL)


def test_pfba_within_the_obj_frac_gap(backend, form):
    """pFBA growth is within [0.999, 1] x the oracle's pFBA growth, and total flux within 0.2%."""
    result, expected = _pfba(backend, form), oracle_pfba()
    growth, optimum = result.fluxes[native(backend, BIOMASS)], expected.fluxes[BIOMASS]
    assert 0.999 * optimum - OBJ_TOL <= growth <= optimum + OBJ_TOL
    assert result.objective_value == pytest.approx(expected.objective_value, rel=2e-3)


@pytest.mark.parametrize("knockout", ["PGI", "PFK"])
def test_lmoma_knockout_matches_cobra_moma(backend, form, knockout):
    sim = load_sim(backend)
    result = sim.simulate(method="lMOMA", constraints={given(knockout, form): 0}, reference=wt_reference(form))
    assert result.status == SStatus.OPTIMAL
    assert result.objective_value == pytest.approx(oracle_lmoma(knockout), abs=OBJ_TOL)


@pytest.mark.parametrize(
    "backend",
    ["cobra", "reframed"],
)
def test_lmoma_default_reference_matches_cobra_moma(backend):
    result = load_sim(backend).simulate(method="lMOMA", constraints={native(backend, "PGI"): 0})
    assert result.objective_value == pytest.approx(oracle_lmoma("PGI"), abs=OBJ_TOL)


def test_reaction_knockouts_match_cobra(backend, form):
    sim = load_sim(backend)
    expected = oracle_deletions("reaction")
    assert len(expected) == len(sim.reactions)
    for r_id, growth in expected.items():
        assert_growth(sim.simulate(constraints={given(r_id, form): 0}, slim=True), growth, r_id)


def test_gene_knockouts_match_cobra(backend, form):
    """Gene knockout via the Simulator API: close what evaluate_gprs() no longer returns."""
    sim = load_sim(backend)
    expected = oracle_deletions("gene")
    genes = sorted(expected)
    assert len(genes) == len(sim.genes)
    for gene in genes:
        active = [given(g, form, "gene") for g in genes if g != gene]
        inactive = set(sim.reactions) - set(sim.evaluate_gprs(active))
        assert_growth(
            sim.simulate(constraints={r_id: 0 for r_id in inactive}, slim=True),
            expected[gene],
            gene,
        )


def test_anaerobic_envcond_matches_cobra(backend, form):
    sim = load_sim(backend, envcond={given("EX_o2_e", form): (0, 1000)})
    assert [canonical(r_id) for r_id in sim.environmental_conditions] == ["EX_o2_e"]
    result = sim.simulate()
    assert result.objective_value == pytest.approx(oracle_model(anaerobic=True).slim_optimize(), abs=OBJ_TOL)
    assert_fixed_fluxes(result.fluxes, backend, oracle_fva(1.0, anaerobic=True), FIXED_ANAEROBIC)


# maintenance far above what 10 mmol/gDW/h of glucose can supply
INFEASIBLE = ("ATPM", (500, 1000))


def test_infeasible_case_is_infeasible_for_cobra():
    model = oracle_model()
    model.reactions.get_by_id(INFEASIBLE[0]).bounds = INFEASIBLE[1]
    assert model.optimize().status == "infeasible"


@pytest.mark.parametrize("method", ["FBA", "pFBA", "lMOMA"])
def test_infeasible_status_and_objective(backend, form, method):
    sim = load_sim(backend)
    constraints = {given(INFEASIBLE[0], form): INFEASIBLE[1]}
    result = sim.simulate(method=method, constraints=constraints, reference=wt_reference(form))
    assert result.status == SStatus.INFEASIBLE
    assert result.objective_value is None


def _infeasible_fba(backend, slim):
    return load_sim(backend).simulate(constraints={native(backend, INFEASIBLE[0]): INFEASIBLE[1]}, slim=slim)


def test_infeasible_slim_value_is_the_same_on_both_backends():
    assert repr(_infeasible_fba("cobra", True)) == repr(_infeasible_fba("reframed", True))


def test_infeasible_fluxes_are_the_same_on_both_backends():
    assert (_infeasible_fba("cobra", False).fluxes is None) == (_infeasible_fba("reframed", False).fluxes is None)


# ----------------------------------------------------------------------------------------------------------------------
# 2. FVA
# ----------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("obj_frac", [0.9, 1.0])
def test_fva_matches_cobra(backend, form, obj_frac):
    expected = oracle_fva(obj_frac)
    reactions = [r_id for r_id in expected.index if r_id not in LOOP]
    ranges = load_sim(backend).FVA(reactions=[given(r_id, form) for r_id in reactions], obj_frac=obj_frac)
    assert sorted(canonical(r_id) for r_id in ranges) == sorted(reactions)
    for r_id, (low, high) in ranges.items():
        base = canonical(r_id)
        bounds = [expected.loc[base, "minimum"], expected.loc[base, "maximum"]]
        assert [finite(low), finite(high)] == pytest.approx(bounds, abs=FLUX_TOL), base


def test_fva_dataframe_has_the_same_layout_on_both_backends():
    frames = {b: load_sim(b).FVA(reactions=[native(b, "PGI")], format="df") for b in BACKENDS}
    assert list(frames["cobra"].columns) == list(frames["reframed"].columns)
    assert [canonical(str(i)) for i in frames["cobra"].index] == [canonical(str(i)) for i in frames["reframed"].index]


# ----------------------------------------------------------------------------------------------------------------------
# 3. Essentiality
# ----------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize("kind", ["gene", "reaction"])
def test_essential_sets_match_cobra(backend, kind):
    sim = load_sim(backend)
    found = sim.essential_genes() if kind == "gene" else sim.essential_reactions()
    expected = {i for i, growth in oracle_deletions(kind).items() if no_growth(growth)}
    assert expected
    assert {canonical(i) for i in found} == expected


# ----------------------------------------------------------------------------------------------------------------------
# 4. Model queries
# ----------------------------------------------------------------------------------------------------------------------
def test_get_reaction_matches_cobra(backend, form):
    sim = load_sim(backend)
    for rxn in oracle_model().reactions:
        r_id = given(rxn.id, form)
        info = sim.get_reaction(r_id)
        assert canonical(info["id"]) == rxn.id
        assert info["name"] == rxn.name
        assert (finite(info["lb"]), finite(info["ub"])) == rxn.bounds, rxn.id
        assert tuple(finite(b) for b in sim.get_reaction_bounds(r_id)) == rxn.bounds, rxn.id
        assert {canonical(m): c for m, c in info["stoichiometry"].items()} == {
            m.id: c for m, c in rxn.metabolites.items()
        }


def test_get_gene_matches_cobra(backend, form):
    sim = load_sim(backend)
    for gene in oracle_model().genes:
        info = sim.get_gene(given(gene.id, form, "gene"))
        assert canonical(info["id"]) == gene.id
        assert info["name"] == gene.name
        assert {canonical(r_id) for r_id in info["reactions"]} == {r.id for r in gene.reactions}


def _gpr_genes(gpr):
    tokens = re.findall(r"\w+", str(gpr) if gpr else "")
    return {canonical(t) for t in tokens if t.lower() not in ("and", "or", "not")}


def test_gpr_matches_cobra(backend, form):
    sim = load_sim(backend)
    model = oracle_model()
    for rxn in model.reactions:
        assert _gpr_genes(sim.get_gpr(given(rxn.id, form))) == {g.id for g in rxn.genes}, rxn.id
    genes = sorted(g.id for g in model.genes)
    for step in (2, 3, 5):  # knock out every step-th gene
        knocked = set(genes[::step])
        active = [given(g, form, "gene") for g in genes if g not in knocked]
        expected = {r.id for r in model.reactions if not r.genes or r.gpr.eval(knockouts=knocked)}
        assert {canonical(r_id) for r_id in sim.evaluate_gprs(active)} == expected, step


def test_set_bounds_matches_cobra(backend, form):
    sim = load_sim(backend)
    sim.set_reaction_bounds(given(UPTAKE, form), -5, 1000)
    assert sim.get_reaction_bounds(given(UPTAKE, form)) == (-5, 1000)
    model = oracle_model()
    model.reactions.get_by_id(UPTAKE).lower_bound = -5
    assert sim.simulate().objective_value == pytest.approx(model.slim_optimize(), abs=OBJ_TOL)


@pytest.mark.parametrize("maximize", [True, False])
@pytest.mark.parametrize("backend", ["cobra", "reframed"])
def test_set_objective_matches_cobra(backend, form, maximize):
    sim = load_sim(backend)
    sim.set_objective(given("EX_o2_e", form))
    assert {canonical(r_id): v for r_id, v in sim.objective.items() if v} == {"EX_o2_e": 1.0}
    model = oracle_model()
    model.objective = "EX_o2_e"
    model.objective_direction = "max" if maximize else "min"
    assert sim.simulate(maximize=maximize).objective_value == pytest.approx(model.slim_optimize(), abs=OBJ_TOL)


def test_objective_setter_matches_cobra(backend, form):
    sim = load_sim(backend)
    sim.objective = given("EX_o2_e", form)
    assert {canonical(r_id): v for r_id, v in sim.objective.items() if v} == {"EX_o2_e": 1.0}
    model = oracle_model()
    model.objective = "EX_o2_e"
    assert sim.simulate().objective_value == pytest.approx(model.slim_optimize(), abs=OBJ_TOL)


def test_metabolite_producers_and_consumers_match_cobra(backend, form):
    sim = load_sim(backend)
    for met in oracle_model().metabolites:
        m_id = given(met.id, form, "metabolite")
        coefficients = {r.id: r.metabolites[met] for r in met.reactions}
        assert {canonical(r) for r in sim.get_metabolite_producers(m_id)} == {
            r for r, c in coefficients.items() if c > 0
        }
        assert {canonical(r) for r in sim.get_metabolite_consumers(m_id)} == {
            r for r, c in coefficients.items() if c < 0
        }


# ----------------------------------------------------------------------------------------------------------------------
# 5. Omics: one deterministic profile; the COBRApy run with unprefixed gene ids is the reference
# ----------------------------------------------------------------------------------------------------------------------
# GIMME's default growth_frac; iMAT is run with growth_frac=0.5
GROWTH_FLOOR = {"GIMME": 0.9, "iMAT": 0.5}


def small_expression(form):
    """Every second gene of e_coli_core (sorted), at levels 1..23."""
    import numpy as np

    from mewpy.omics import ExpressionSet

    genes = sorted(g.id for g in _cobra_core().genes)[::2]
    levels = np.array([[float((7 * i) % 23 + 1)] for i in range(len(genes))])
    return ExpressionSet([given(g, form, "gene") for g in genes], ["c1"], levels)


@lru_cache(maxsize=None)
def _omics(backend, form, method):
    from mewpy.omics import GIMME, Preprocessing, eFlux, iMAT

    sim, expr = load_sim(backend), small_expression(form)
    n_expressed = len(Preprocessing(sim, expr).reactions_expression(0))
    if method == "eFlux":
        result = eFlux(sim, expr)
    elif method == "GIMME":
        result = GIMME(sim, expr)
    else:
        result = iMAT(sim, expr, growth_frac=GROWTH_FLOOR["iMAT"])
    return (
        n_expressed,
        result.status,
        result.objective_value,
        result.fluxes[native(backend, BIOMASS)],
    )


@pytest.mark.parametrize("method", ["eFlux", "GIMME", "iMAT"])
def test_omics_objective_and_growth(backend, form, method):
    n_ref, status_ref, objective_ref, _ = _omics("cobra", "unprefixed", method)
    n_expressed, status, objective, growth = _omics(backend, form, method)
    assert n_ref > 0
    assert n_expressed == n_ref
    assert status == status_ref == SStatus.OPTIMAL
    assert objective == pytest.approx(objective_ref, abs=OBJ_TOL)
    if method == "eFlux":
        # eFlux maximises growth under the expression-scaled bounds
        assert growth == pytest.approx(objective, abs=OBJ_TOL)
    else:
        # GIMME and iMAT growth is not unique: only the floor they impose is checked
        assert growth >= GROWTH_FLOOR[method] * oracle_model().slim_optimize() - OBJ_TOL


# ----------------------------------------------------------------------------------------------------------------------
# 6. Problems: anaerobic succinate, evaluated with FBA; expected fitness from COBRApy
# ----------------------------------------------------------------------------------------------------------------------
WYIELD_ALPHA, WYIELD_OBJ_FRAC = 0.3, 0.99  # WYIELD defaults
RKO_TARGETS = ["ACALD", "ACKr", "ALCD2x", "LDH_D", "PFL", "PTAr"]
RKO_CANDIDATE = ["ALCD2x", "LDH_D"]
GKO_TARGETS = ["b0351", "b0902", "b1241", "b1380", "b2133"]
GKO_CANDIDATE = ["b0351", "b1241", "b1380", "b2133"]  # closes ACALD and LDH_D
ROU_TARGETS = ["ACKr", "ALCD2x", "LDH_D", "PTAr"]
ROU_CANDIDATE = {"ALCD2x": 0, "LDH_D": 0, "ACKr": 0.5}
EA_TARGETS = ["ACALD", "ACKr", "ALCD2x", "FUM", "LDH_D", "MDH", "ME1", "ME2", "PFL", "PPC"]
EA_TARGETS += ["PPCK", "PTAr"]


def evaluation_functions(form):
    from mewpy.optimization.evaluation import BPCY, WYIELD, TargetFlux

    biomass, product = given(BIOMASS, form), given(PRODUCT, form)
    return [
        BPCY(biomass, product, uptake=given(UPTAKE, form), method="FBA"),
        WYIELD(biomass, product),
        TargetFlux(product, biomass=biomass, method="FBA"),
    ]


def oracle_fitness(model):
    """[BPCY, WYIELD, TargetFlux] for a COBRApy model that carries the candidate."""
    from cobra.flux_analysis import flux_variability_analysis

    growth = model.slim_optimize()
    optimum = flux_variability_analysis(model, [PRODUCT, UPTAKE], fraction_of_optimum=1.0, processes=1)
    for r_id in (PRODUCT, UPTAKE):
        assert optimum.loc[r_id, "maximum"] - optimum.loc[r_id, "minimum"] < 1e-6, f"{r_id} is not fixed"
    product, uptake = optimum.loc[PRODUCT, "minimum"], optimum.loc[UPTAKE, "minimum"]
    near = flux_variability_analysis(model, [PRODUCT], fraction_of_optimum=WYIELD_OBJ_FRAC, processes=1)
    wyield = WYIELD_ALPHA * near.loc[PRODUCT, "maximum"] + (1 - WYIELD_ALPHA) * near.loc[PRODUCT, "minimum"]
    return [growth * product / abs(uptake), wyield, product]


@pytest.mark.parametrize("kind", ["RKO", "GKO"])
def test_ko_problem_fitness_matches_cobra(backend, form, kind):
    from mewpy.problems import GKOProblem, RKOProblem

    model = oracle_model(anaerobic=True)
    if kind == "RKO":
        cls, id_kind, targets, candidate = RKOProblem, "reaction", RKO_TARGETS, RKO_CANDIDATE
        for r_id in candidate:
            model.reactions.get_by_id(r_id).knock_out()
    else:
        cls, id_kind, targets, candidate = GKOProblem, "gene", GKO_TARGETS, GKO_CANDIDATE
        for g_id in candidate:
            model.genes.get_by_id(g_id).knock_out()
    problem = cls(
        load_model(backend),
        evaluation_functions(form),
        target=[given(t, form, id_kind) for t in targets],
        envcond={given("EX_o2_e", form): (0, 1000)},
    )
    problem.pre_process()
    fitness = problem.evaluate_solution(problem.encode([given(c, form, id_kind) for c in candidate]))
    expected = oracle_fitness(model)
    assert expected[0] > 0
    assert fitness == pytest.approx(expected, abs=OBJ_TOL)


def _rou(backend, form, reference):
    from mewpy.problems import ROUProblem

    kwargs = {} if reference is None else {"reference": reference}
    problem = ROUProblem(
        load_model(backend),
        evaluation_functions(form),
        target=[given(t, form) for t in ROU_TARGETS],
        envcond={given("EX_o2_e", form): (0, 1000)},
        **kwargs,
    )
    problem.pre_process()
    solution = problem.encode({given(r_id, form): level for r_id, level in ROU_CANDIDATE.items()})
    constraints = problem.solution_to_constraints(problem.decode(solution))
    return constraints, problem.evaluate_solution(solution)


@lru_cache(maxsize=None)
def rou_oracle():
    """The ROU constraints derived from COBRApy's anaerobic pFBA, and their COBRApy fitness."""
    acetate = oracle_pfba(anaerobic=True).fluxes["ACKr"]
    assert acetate < 0  # level 0.5 of a negative reference flux gives (0.5 * flux, 0)
    constraints = {"ALCD2x": (0, 0), "LDH_D": (0, 0), "ACKr": (0.5 * acetate, 0)}
    model = oracle_model(anaerobic=True)
    for r_id, bounds in constraints.items():
        model.reactions.get_by_id(r_id).bounds = bounds
    return constraints, oracle_fitness(model)


def test_rou_problem_fitness_matches_cobra(backend, form):
    expected_constraints, expected_fitness = rou_oracle()
    constraints, fitness = _rou(backend, form, wt_reference(form, anaerobic=True))
    assert sorted(canonical(r_id) for r_id in constraints) == sorted(expected_constraints)
    for r_id, bounds in constraints.items():
        assert list(bounds) == pytest.approx(list(expected_constraints[canonical(r_id)]), abs=FLUX_TOL), r_id
    assert expected_fitness[0] > 0
    assert fitness == pytest.approx(expected_fitness, abs=OBJ_TOL)


@pytest.mark.parametrize(
    "backend",
    ["cobra", "reframed"],
)
def test_rou_default_reference_matches_cobra(backend):
    _, fitness = _rou(backend, "unprefixed", None)
    assert fitness == pytest.approx(rou_oracle()[1], abs=OBJ_TOL)


@pytest.mark.parametrize("engine", ["inspyred", "jmetal"])
def test_seeded_ea_gives_the_same_best_fitness(engine):
    from mewpy.optimization import EA, get_default_engine, set_default_engine
    from mewpy.optimization.evaluation import WYIELD
    from mewpy.problems import RKOProblem

    previous = get_default_engine()
    set_default_engine(engine)
    best = {}
    try:
        for backend in BACKENDS:
            problem = RKOProblem(
                load_model(backend),
                [WYIELD(BIOMASS, PRODUCT)],
                target=EA_TARGETS,
                envcond={"EX_o2_e": (0, 1000)},
                candidate_max_size=3,
            )
            ea = EA(problem, max_generations=2, mp=False, population_size=6, seed=3, algorithm="GA")
            best[backend] = max(solution.fitness[0] for solution in ea.run(simplify=False))
    finally:
        set_default_engine(previous)
    assert best["cobra"] > 0
    assert best["reframed"] == pytest.approx(best["cobra"], abs=OBJ_TOL)


# ----------------------------------------------------------------------------------------------------------------------
# 7. GERM: RegulatoryExtension.from_sbml(flavor='cobra') vs flavor='reframed'
# ----------------------------------------------------------------------------------------------------------------------
@pytest.fixture(params=BACKENDS)
def flavor(request):
    return request.param


def regulatory_model(flavor, trn=EC_CORE_TRN):
    from mewpy.germ.models import RegulatoryExtension

    return RegulatoryExtension.from_sbml(EC_CORE, str(trn), regulatory_format="csv", flavor=flavor, **TRN_KWARGS)


def test_germ_wild_type_matches_docs(flavor):
    from mewpy.germ.analysis import RFBA, SRFBA

    model = regulatory_model(flavor)
    assert RFBA(model).build().optimize().objective_value == pytest.approx(RFBA_DOCS_WT, abs=OBJ_TOL)
    assert SRFBA(model).build().optimize().objective_value == pytest.approx(SRFBA_DOCS_WT, abs=OBJ_TOL)


@pytest.mark.parametrize("method", ["rfba", "srfba"])
def test_germ_regulator_deletions_match_across_flavors(method):
    from mewpy.germ.analysis import isingle_regulator_deletion

    tables = {}
    for flavor in BACKENDS:
        table = isingle_regulator_deletion(regulatory_model(flavor), method=method)
        tables[flavor] = {canonical(regulator): growth or 0.0 for regulator, growth in table["growth"].items()}
    assert sorted(tables["reframed"]) == sorted(tables["cobra"])
    for regulator, growth in tables["cobra"].items():
        assert tables["reframed"][regulator] == pytest.approx(growth, abs=OBJ_TOL), regulator
    wild_type = RFBA_DOCS_WT if method == "rfba" else SRFBA_DOCS_WT
    # most knockouts change nothing; some lower growth
    assert statistics.median(tables["cobra"].values()) == pytest.approx(wild_type, abs=OBJ_TOL)
    assert min(tables["cobra"].values()) < wild_type - 1e-3


def test_germ_prom_and_coregflux_match_across_flavors():
    from mewpy.germ.analysis import PROM, CoRegFlux

    results = {}
    for flavor in BACKENDS:
        model = regulatory_model(flavor)
        probabilities = {}
        for _, interaction in model.yield_interactions():
            for regulator in interaction.yield_regulators():
                probabilities[(interaction.target.id, regulator.id)] = 0.0
        prom = PROM(model).build()
        prom = prom.optimize(initial_state=probabilities, regulators=["b4401", "b3261"])
        coregflux = CoRegFlux(model).build()
        coregflux = coregflux.optimize(initial_state={g_id: 10.0 for g_id in model.genes})
        results[flavor] = [prom.solutions[f"ko_{r}"].objective_value for r in ("b4401", "b3261")]
        results[flavor].append(coregflux.objective_value)
    assert max(results["cobra"]) < oracle_model().slim_optimize() - 1e-3  # all below the FBA WT
    assert results["reframed"] == pytest.approx(results["cobra"], abs=OBJ_TOL)


def test_germ_fva_matches_cobra(flavor):
    from mewpy.germ.analysis import fva

    expected = oracle_fva(0.9)
    reactions = [r_id for r_id in expected.index if r_id not in LOOP]
    result = fva(
        regulatory_model(flavor),
        fraction=0.9,
        reactions=[native(flavor, r_id) for r_id in reactions],
    )
    assert sorted(canonical(r_id) for r_id in result.index) == sorted(reactions)
    for r_id in result.index:
        found = [finite(result.loc[r_id, "minimum"]), finite(result.loc[r_id, "maximum"])]
        base = canonical(r_id)
        assert found == pytest.approx([expected.loc[base, "minimum"], expected.loc[base, "maximum"]], abs=FLUX_TOL)


@pytest.fixture(scope="module")
def conflicting_trn(tmp_path_factory):
    """The TRN with two essential genes, b1779 (GAPD) and b2926 (PGK), repressed by Fis (b3261)."""
    with open(EC_CORE_TRN, newline="") as handle:
        rows = list(csv.reader(handle))
    path = tmp_path_factory.mktemp("trn").joinpath("trn_conflicts.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows([row[0], row[1], "NOT b3261"] if row[0] in ("b1779", "b2926") else row for row in rows)
    return path


def test_germ_find_conflicts_match_across_flavors(flavor, conflicting_trn):
    from mewpy.germ.analysis import find_conflicts

    genes, reactions = find_conflicts(regulatory_model(flavor))
    assert genes.empty and reactions.empty
    essential = {g for g, growth in oracle_deletions("gene").items() if no_growth(growth)}
    assert {"b1779", "b2926"} <= essential
    genes, reactions = find_conflicts(regulatory_model(flavor, conflicting_trn))
    assert {canonical(g_id) for g_id in genes.index} == {"b1779", "b2926"}
    assert {canonical(r_id) for r_id in reactions.index} == {"GAPD", "PGK"}


def _germ_case(case, model, form):
    from mewpy.germ.analysis import fva, slim_rfba, slim_srfba

    if case == "gene_knockout":
        return slim_srfba(model, initial_state={given("b1779", form, "gene"): 0.0})
    if case == "constraint":
        return slim_srfba(model, constraints={given("EX_o2_e", form): (0, 1000)})
    if case == "regulator_knockout":  # 'o2_e' is the TRN's id for the oxygen regulator
        return slim_rfba(model, initial_state={given("o2_e", form, "metabolite"): 0.0})
    table = fva(model, fraction=0.9, reactions=[given("PGI", form)])
    return [table.iloc[0]["minimum"], table.iloc[0]["maximum"]]


@pytest.mark.parametrize("case", ["gene_knockout", "constraint", "regulator_knockout", "fva"])
@pytest.mark.parametrize(
    "flavor, form",
    [
        ("cobra", "prefixed"),
        ("reframed", "unprefixed"),
        ("reframed", "prefixed"),
    ],
)
def test_germ_user_ids_in_either_form(case, flavor, form):
    """Expected: COBRApy's own result for the gene knockout and FVA; otherwise the COBRApy flavor
    with its native (unprefixed) ids."""
    if case == "gene_knockout":  # b1779 is GAPD's only gene
        expected = oracle_deletions("gene")["b1779"]
    elif case == "fva":
        expected = [oracle_fva(0.9).loc["PGI", "minimum"], oracle_fva(0.9).loc["PGI", "maximum"]]
    else:
        expected = _germ_case(case, regulatory_model("cobra"), "unprefixed")
    assert _germ_case(case, regulatory_model(flavor), form) == pytest.approx(expected, abs=FLUX_TOL)


# ----------------------------------------------------------------------------------------------------------------------
# 8. Community: two e_coli_core copies; the cobra+cobra community is the reference
# ----------------------------------------------------------------------------------------------------------------------
PAIRS = {
    "cobra+cobra": ("cobra", "cobra"),
    "reframed+reframed": ("reframed", "reframed"),
    "cobra+reframed": BACKENDS,
}


def community(pair, knockout=None, flavor="reframed"):
    """Members m1 and m2; ``knockout`` closes a reaction of m1."""
    from mewpy.com import CommunityModel

    models = []
    for org_id, backend in zip(("m1", "m2"), PAIRS[pair]):
        model = load_model(backend)
        if knockout and org_id == "m1":
            get_simulator(model).set_reaction_bounds(native(backend, knockout), 0, 0)
        model.id = org_id
        models.append(model)
    return CommunityModel(models, verbose=False, flavor=flavor, merge_biomasses=False)


def glucose_medium(com):
    """e_coli_core's medium (glucose at 10) on the community exchanges; the others only secrete."""
    from mewpy.simulation import Environment

    constraints = {r_id: (0, 1000) for r_id in com.merged_model.get_exchange_reactions()}
    for r_id, bounds in Environment.from_model(oracle_model()).simplify().items():
        constraints[native(com.flavor, r_id)] = bounds
    return constraints


@lru_cache(maxsize=None)
def _steadycom(pair, flavor="reframed"):
    """m1 cannot take up glucose (GLCpts closed); SteadyComVA(obj_frac=1) pins it to within 3e-5."""
    from mewpy.com.steadycom import SteadyCom

    com = community(pair, knockout="GLCpts", flavor=flavor)
    solution = SteadyCom(com, constraints=glucose_medium(com))
    return solution.growth, dict(solution.abundance)


@pytest.mark.parametrize("flavor", ["reframed", "cobra"])
@pytest.mark.parametrize("pair", list(PAIRS))
def test_steadycom_growth_and_abundance_match_across_backends(pair, flavor):
    expected_growth, expected_abundance = _steadycom("cobra+cobra")
    growth, abundance = _steadycom(pair, flavor)
    assert growth > 0.1
    assert growth == pytest.approx(expected_growth, abs=OBJ_TOL)
    assert sum(abundance.values()) == pytest.approx(1.0, abs=OBJ_TOL)
    assert abundance == pytest.approx(expected_abundance, abs=FLUX_TOL)


@lru_cache(maxsize=None)
def _steadycomva(pair):
    from mewpy.com.steadycom import SteadyComVA

    com = community(pair)
    return {org_id: list(bounds) for org_id, bounds in SteadyComVA(com, 0.9, glucose_medium(com)).items()}


@pytest.mark.parametrize("pair", ["reframed+reframed", "cobra+reframed"])
def test_steadycomva_ranges_match_across_backends(pair):
    expected = _steadycomva("cobra+cobra")
    found = _steadycomva(pair)
    assert sorted(found) == sorted(expected)
    for org_id, bounds in expected.items():
        assert found[org_id] == pytest.approx(bounds, abs=FLUX_TOL), org_id


@pytest.mark.parametrize("pair", ["reframed+reframed", "cobra+reframed"])
def test_steadycom_without_constraints_growth_matches_across_backends(pair):
    from mewpy.com.steadycom import SteadyCom

    expected = SteadyCom(community("cobra+cobra")).growth
    assert SteadyCom(community(pair)).growth == pytest.approx(expected, abs=OBJ_TOL)


# ----------------------------------------------------------------------------------------------------------------------
# 9. Environment
# ----------------------------------------------------------------------------------------------------------------------
MEDIUM = ["glc__D", "o2", "nh4", "pi", "h2o", "h", "co2"]


@pytest.mark.parametrize("anaerobic", [False, True])
@pytest.mark.parametrize("target", ["simulator", "model", "constraints"])
def test_environment_from_compounds_matches_cobra_medium(backend, target, anaerobic):
    """from_compounds(model=...); apply() to a simulator, a raw model, or with inplace=False."""
    from mewpy.simulation import Environment

    compounds = [c for c in MEDIUM if not (anaerobic and c == "o2")]
    model = load_model(backend)
    if target == "model":
        env = Environment.from_compounds(compounds, model=model)
        env.apply(model)
        sim = get_simulator(model)
    else:
        sim = get_simulator(model)
        env = Environment.from_compounds(compounds, model=sim)
        if target == "simulator":
            env.apply(sim)
        else:
            sim = load_sim(backend, envcond=env.apply(sim, inplace=False))
    assert [canonical(r_id) for r_id in env] == [f"EX_{c}_e" for c in compounds]
    oracle = oracle_model()
    oracle.medium = {f"EX_{c}_e": 10.0 for c in compounds}
    assert sim.simulate().objective_value == pytest.approx(oracle.slim_optimize(), abs=OBJ_TOL)
