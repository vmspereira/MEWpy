"""
Regression tests for the GERM analysis/algebra review findings (F1-F12 and two confirmed extras).

Each test pins a value, not just "it returns something": FVA ranges against COBRApy, essential-gene
knockouts under SRFBA giving zero growth, PROM/CoRegFlux knockouts lowering growth, and no mutation
of caller-owned dictionaries.
"""

import csv
from pathlib import Path

import pytest

DATA = Path(__file__).parent.joinpath("data")
EC_CORE = DATA.joinpath("e_coli_core.xml")
EC_CORE_TRN = DATA.joinpath("e_coli_core_trn.csv")
TRN_KWARGS = dict(sep=",", id_col=0, rule_col=2, aliases_cols=[1], header=0)
FVA_REACTIONS = ["PGI", "EX_glc__D_e", "ACKr", "Biomass_Ecoli_core"]
TOL = 1e-4


# ----------------------------------------------------------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def cobra_model():
    cobra = pytest.importorskip("cobra")
    return cobra.io.read_sbml_model(str(EC_CORE))


@pytest.fixture(scope="module")
def wt_growth(cobra_model):
    return cobra_model.slim_optimize()


@pytest.fixture(scope="module")
def legacy_metabolic():
    from mewpy.io import Engines, Reader, read_model

    model = read_model(Reader(Engines.MetabolicSBML, EC_CORE))
    model.objective = {"Biomass_Ecoli_core": 1}
    return model


@pytest.fixture(scope="module")
def legacy_integrated():
    from mewpy.io import Engines, Reader, read_model

    met = Reader(Engines.MetabolicSBML, EC_CORE)
    reg = Reader(Engines.BooleanRegulatoryCSV, EC_CORE_TRN, **TRN_KWARGS)
    model = read_model(reg, met)
    model.objective = {"Biomass_Ecoli_core": 1}
    return model


def _regulatory_extension(trn_path):
    from mewpy.germ.models import RegulatoryExtension

    return RegulatoryExtension.from_sbml(str(EC_CORE), str(trn_path), regulatory_format="csv", **TRN_KWARGS)


@pytest.fixture(scope="module")
def re_model():
    return _regulatory_extension(EC_CORE_TRN)


@pytest.fixture(scope="module")
def re_reduced_trn(tmp_path_factory):
    """The e_coli_core TRN without its constitutive ('TRUE') rows: unregulated genes are absent from the TRN."""
    rows = list(csv.reader(open(EC_CORE_TRN)))
    path = tmp_path_factory.mktemp("trn").joinpath("trn_reduced.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(rows[0])
        writer.writerows(row for row in rows[1:] if row[2].strip().upper() != "TRUE")
    return _regulatory_extension(path)


def _all_zero_probabilities(model):
    probabilities = {}
    for item in model.yield_interactions():
        interaction = item[1] if isinstance(item, tuple) else item
        for regulator in interaction.yield_regulators():
            probabilities[(interaction.target.id, regulator.id)] = 0.0
    return probabilities


def _gene_ids(model):
    return list(model.genes.keys()) if isinstance(model.genes, dict) else list(model.genes)


# ----------------------------------------------------------------------------------------------------------------------
# F1 - SRFBA MILP rows must reach the solver
# ----------------------------------------------------------------------------------------------------------------------
def test_f1_srfba_build_flushes_pending_rows(re_model):
    from mewpy.germ.analysis import SRFBA

    srfba = SRFBA(re_model).build()
    assert not getattr(srfba.solver, "_cached_vars", [])
    assert not getattr(srfba.solver, "_cached_constrs", [])


def test_f1_srfba_essential_gene_ko_is_lethal_regulatory_extension(re_model, wt_growth):
    from mewpy.germ.analysis import SRFBA

    wt = SRFBA(re_model).build().optimize(to_solver=True)
    assert wt.fobj == pytest.approx(wt_growth, abs=TOL)
    # G_b1779 is the only gene of GAPD (essential in e_coli_core)
    ko = SRFBA(re_model).build().optimize(initial_state={"G_b1779": 0.0}, to_solver=True)
    assert ko.fobj == pytest.approx(0.0, abs=TOL)


def test_f1_srfba_essential_gene_ko_is_lethal_legacy(legacy_integrated, wt_growth):
    from mewpy.germ.analysis import SRFBA

    ko = SRFBA(legacy_integrated).build().optimize(initial_state={"b1779": 0.0}, to_solver=True)
    assert ko.fobj == pytest.approx(0.0, abs=TOL)


# ----------------------------------------------------------------------------------------------------------------------
# F2 - per-reaction objectives must be honoured (FVA, iFVA, slim_*)
# ----------------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def cobra_fva(cobra_model):
    from cobra.flux_analysis import flux_variability_analysis

    return flux_variability_analysis(cobra_model, reaction_list=FVA_REACTIONS, fraction_of_optimum=0.9)


def _assert_fva_equal(result, expected):
    for rxn in FVA_REACTIONS:
        assert result.loc[rxn, "minimum"] == pytest.approx(expected.loc[rxn, "minimum"], abs=TOL), rxn
        assert result.loc[rxn, "maximum"] == pytest.approx(expected.loc[rxn, "maximum"], abs=TOL), rxn


def test_f2_fva_matches_cobra(legacy_metabolic, cobra_fva):
    from mewpy.germ.analysis import fva

    _assert_fva_equal(fva(legacy_metabolic, fraction=0.9, reactions=FVA_REACTIONS), cobra_fva)


@pytest.mark.parametrize("method", ["rfba", "srfba"])
def test_f2_ifva_matches_cobra_on_metabolic_model(legacy_metabolic, cobra_fva, method):
    from mewpy.germ.analysis import ifva

    _assert_fva_equal(ifva(legacy_metabolic, fraction=0.9, reactions=FVA_REACTIONS, method=method), cobra_fva)


@pytest.mark.parametrize("slim", ["slim_rfba", "slim_srfba"])
def test_f2_slim_methods_honour_objective(legacy_metabolic, cobra_model, slim):
    import mewpy.germ.analysis as analysis

    with cobra_model as tmp:
        tmp.objective = "PGI"
        expected = tmp.slim_optimize()
    assert getattr(analysis, slim)(legacy_metabolic, objective="PGI") == pytest.approx(expected, abs=TOL)


# ----------------------------------------------------------------------------------------------------------------------
# F3 - PROM knockouts must change reaction bounds
# ----------------------------------------------------------------------------------------------------------------------
def test_f3_prom_ko_lowers_growth_legacy(legacy_integrated, wt_growth):
    from mewpy.germ.analysis import PROM

    probabilities = _all_zero_probabilities(legacy_integrated)
    solution = PROM(legacy_integrated).build().optimize(initial_state=probabilities, regulators=["b4401"])
    growth = next(iter(solution.solutions.values())).objective_value
    assert growth < wt_growth - 1e-3


def test_f3_prom_ko_lowers_growth_regulatory_extension(re_model, wt_growth):
    """Needs TRN/model gene-id reconciliation in RegulatoryExtension (b4401 targets are 'b...' vs 'G_b...')."""
    from mewpy.germ.analysis import PROM

    probabilities = _all_zero_probabilities(re_model)
    solution = PROM(re_model).build().optimize(initial_state=probabilities, regulators=["b4401"])
    growth = next(iter(solution.solutions.values())).objective_value
    assert growth < wt_growth - 1e-3


# ----------------------------------------------------------------------------------------------------------------------
# F4 - CoRegFlux must apply the gene state (legacy models: reaction keys were objects)
# ----------------------------------------------------------------------------------------------------------------------
def test_f4_coregflux_legacy_applies_gene_state(legacy_integrated, re_model, wt_growth):
    from mewpy.germ.analysis import CoRegFlux

    legacy = (
        CoRegFlux(legacy_integrated).build().optimize(initial_state={g: 10.0 for g in _gene_ids(legacy_integrated)})
    )
    re = CoRegFlux(re_model).build().optimize(initial_state={g: 10.0 for g in _gene_ids(re_model)})
    assert legacy.objective_value < wt_growth - 1e-3
    assert legacy.objective_value == pytest.approx(re.objective_value, abs=TOL)


# ----------------------------------------------------------------------------------------------------------------------
# F5 - continuous GPR evaluation must find genes in nested GPRs
# ----------------------------------------------------------------------------------------------------------------------
def test_f5_continuous_gpr_handles_nested_gprs(re_model):
    from mewpy.germ.analysis.analysis_utils import continuous_gpr

    with_gpr = [r for r in re_model.reactions if not re_model.get_parsed_gpr(r).is_none]
    states = continuous_gpr(re_model, {g: 0.5 for g in _gene_ids(re_model)})
    assert set(states) == set(with_gpr)
    # ATPS4r = (8 genes AND) OR (9 genes AND)
    assert states["R_ATPS4r"] == pytest.approx(0.5)


# ----------------------------------------------------------------------------------------------------------------------
# F6 - RFBA: genes absent from the TRN are not regulated, so they stay active
# ----------------------------------------------------------------------------------------------------------------------
def test_f6_rfba_genes_absent_from_trn_stay_active(re_reduced_trn):
    from mewpy.germ.analysis import RFBA

    rfba = RFBA(re_reduced_trn).build()
    state = rfba.decode_metabolic_state(rfba.initial_state())
    all_active = {g: 1.0 for g in _gene_ids(re_reduced_trn)}
    all_active.update(state)
    assert rfba.decode_constraints(state) == rfba.decode_constraints(all_active)


# ----------------------------------------------------------------------------------------------------------------------
# F7 - slim_coregflux / slim_prom
# ----------------------------------------------------------------------------------------------------------------------
def test_f7_slim_coregflux_returns_objective(re_model):
    from mewpy.germ.analysis import CoRegFlux, slim_coregflux

    state = {g: 10.0 for g in _gene_ids(re_model)}
    expected = CoRegFlux(re_model).build().optimize(initial_state=state).objective_value
    assert expected > 0
    assert slim_coregflux(re_model, initial_state=state) == pytest.approx(expected, abs=TOL)


def test_f7_slim_coregflux_accepts_objective(re_model):
    from mewpy.germ.analysis import slim_coregflux

    state = {g: 10.0 for g in _gene_ids(re_model)}
    assert slim_coregflux(re_model, initial_state=state, objective="R_PGI") > 0


def test_f7_slim_prom_accepts_objective(legacy_integrated):
    from mewpy.germ.analysis import slim_prom

    probabilities = _all_zero_probabilities(legacy_integrated)
    default = slim_prom(legacy_integrated, initial_state=probabilities, regulator="b4401")
    explicit = slim_prom(
        legacy_integrated, initial_state=probabilities, regulator="b4401", objective="Biomass_Ecoli_core"
    )
    assert explicit == pytest.approx(default, abs=TOL)


# ----------------------------------------------------------------------------------------------------------------------
# F8 - CoRegFlux must not report an infeasible step as optimal
# ----------------------------------------------------------------------------------------------------------------------
def test_f8_coregflux_reports_infeasible_status(re_model):
    from mewpy.germ.analysis import CoRegFlux
    from mewpy.solvers.solution import Status

    solution = CoRegFlux(re_model).build().optimize(initial_state={g: 0.5 for g in _gene_ids(re_model)})
    assert solution.status != Status.OPTIMAL


# ----------------------------------------------------------------------------------------------------------------------
# F9 - PROM max rates: FVA at 99% of wild-type growth, without leaking into caller constraints
# ----------------------------------------------------------------------------------------------------------------------
def test_f9_prom_max_rates_use_99pct_growth(re_model, cobra_model):
    from cobra.flux_analysis import flux_variability_analysis

    from mewpy.germ.analysis import PROM

    expected = flux_variability_analysis(cobra_model, reaction_list=["PDH"], fraction_of_optimum=0.99)
    prom = PROM(re_model).build()
    prom.solver.solve(linear=prom._linear_objective, minimize=False)
    rates = prom._max_rates(solver_kwargs={"get_values": True})
    assert rates["R_PDH"] == pytest.approx(expected.loc["PDH", "maximum"], abs=1e-3)


def test_f9_prom_does_not_mutate_caller_constraints(re_model):
    from mewpy.germ.analysis import PROM

    user = {"R_EX_glc__D_e": (-10.0, 1000.0)}
    PROM(re_model).build().optimize(
        initial_state=_all_zero_probabilities(re_model), regulators=["b4401"], solver_kwargs={"constraints": user}
    )
    assert user == {"R_EX_glc__D_e": (-10.0, 1000.0)}


# ----------------------------------------------------------------------------------------------------------------------
# F10 - caller-owned constraint dicts must not be mutated
# ----------------------------------------------------------------------------------------------------------------------
def test_f10_fva_does_not_mutate_constraints(legacy_metabolic):
    from mewpy.germ.analysis import fva

    user = {"EX_glc__D_e": (-10.0, 1000.0)}
    fva(legacy_metabolic, reactions=["PGI"], constraints=user)
    assert user == {"EX_glc__D_e": (-10.0, 1000.0)}


def test_f10_ifva_does_not_mutate_constraints(legacy_metabolic):
    from mewpy.germ.analysis import ifva

    user = {"EX_glc__D_e": (-10.0, 1000.0)}
    ifva(legacy_metabolic, reactions=["PGI"], constraints=user)
    assert user == {"EX_glc__D_e": (-10.0, 1000.0)}


def test_f10_isingle_gene_deletion_does_not_accumulate_knockouts(re_model, wt_growth):
    from mewpy.germ.analysis import isingle_gene_deletion

    user = {"R_EX_glc__D_e": (-10.0, 1000.0)}
    result = isingle_gene_deletion(re_model, genes=["G_b1779", "G_b0008"], constraints=user)
    assert user == {"R_EX_glc__D_e": (-10.0, 1000.0)}
    assert result.loc["G_b1779", "growth"] == pytest.approx(0.0, abs=TOL)
    # G_b0008 (TALA, has an isozyme) is not essential; it must not inherit the previous GAPD knockout
    assert result.loc["G_b0008", "growth"] == pytest.approx(wt_growth, abs=TOL)


# ----------------------------------------------------------------------------------------------------------------------
# F11 - crashes on realistic inputs
# ----------------------------------------------------------------------------------------------------------------------
def test_f11_soft_plus_does_not_overflow():
    from mewpy.germ.analysis.analysis_utils import lb_soft_plus, ub_soft_plus

    assert ub_soft_plus(0, 800.0) == pytest.approx(800.0)
    assert lb_soft_plus(0, 800.0) == pytest.approx(-800.0)
    assert ub_soft_plus(0, 1.0) == pytest.approx(1.3132616875182228)


def test_f11_coregflux_scale_true(re_model):
    from mewpy.germ.analysis import CoRegFlux

    solution = CoRegFlux(re_model).build().optimize(initial_state={g: 10.0 for g in _gene_ids(re_model)}, scale=True)
    assert solution.objective_value is not None


def test_f11_coregflux_without_initial_state(re_model, wt_growth):
    from mewpy.germ.analysis import CoRegFlux

    assert CoRegFlux(re_model).build().optimize().objective_value == pytest.approx(wt_growth, abs=TOL)


def test_f11_fva_default_objective_regulatory_extension(re_model, cobra_fva):
    from mewpy.germ.analysis import fva

    result = fva(re_model, fraction=0.9, reactions=["R_PGI"])
    assert result.loc["R_PGI", "minimum"] == pytest.approx(cobra_fva.loc["PGI", "minimum"], abs=TOL)
    assert result.loc["R_PGI", "maximum"] == pytest.approx(cobra_fva.loc["PGI", "maximum"], abs=TOL)


def test_f11_pfba_regulatory_extension(re_model, wt_growth):
    from mewpy.germ.analysis.pfba import _pFBA

    from mewpy.util.constants import ModelConstants

    solution = _pFBA(re_model).optimize()
    fraction = getattr(ModelConstants, "PFBA_OBJ_FRAC", 0.999)
    assert solution.values["R_Biomass_Ecoli_core"] == pytest.approx(fraction * wt_growth, abs=1e-6)


def test_f11_single_gene_deletion_regulatory_extension(re_model, wt_growth):
    from mewpy.germ.analysis import single_gene_deletion

    result = single_gene_deletion(re_model, genes=["G_b1779", "G_b0008"])
    assert result.loc["G_b1779", "growth"] == pytest.approx(0.0, abs=TOL)
    assert result.loc["G_b0008", "growth"] == pytest.approx(wt_growth, abs=TOL)


def test_f11_find_conflicts_regulatory_extension(re_model):
    from mewpy.germ.analysis import find_conflicts

    repressed_genes, repressed_reactions = find_conflicts(re_model)
    assert repressed_genes is not None and repressed_reactions is not None


# ----------------------------------------------------------------------------------------------------------------------
# Extras confirmed during review
# ----------------------------------------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "expression, values, expected",
    [
        ("a and b", {"a": 2, "b": 1}, 1),
        ("a and b", {"a": 0.5, "b": 1}, 1),
        ("a and b", {"a": 2, "b": 0}, 0),
        ("a or b", {"a": 0.5, "b": 0}, 1),
        ("a or b", {"a": 0, "b": 0}, 0),
        ("not a", {"a": 2}, 0),
        ("not a", {"a": 0.5}, 0),
        ("not a", {"a": 0}, 1),
        ("not a", {"a": 1}, 0),
    ],
)
def test_boolean_operators_on_non_binary_values(expression, values, expected):
    from mewpy.germ.algebra import parse_expression

    assert parse_expression(expression).evaluate(values=values) == expected


@pytest.mark.parametrize("minimize, expected_attr", [(True, "REACTION_LOWER_BOUND"), (False, "REACTION_UPPER_BOUND")])
def test_unbounded_objective_sign_follows_direction(minimize, expected_attr):
    from mewpy.germ.analysis.analysis_utils import run_method_and_decode
    from mewpy.solvers.solution import Solution, Status
    from mewpy.util.constants import ModelConstants

    class _Unbounded:
        def optimize(self, **kwargs):
            return Solution(status=Status.UNBOUNDED)

    value, _ = run_method_and_decode(_Unbounded(), objective={"r": 1}, minimize=minimize)
    assert value == getattr(ModelConstants, expected_attr)


# ----------------------------------------------------------------------------------------------------------------------
# Found while fixing: caller constraints overridden, RFBA knockouts ignored, PROM without regulators
# ----------------------------------------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def low_glucose_growth(cobra_model):
    with cobra_model as tmp:
        tmp.reactions.get_by_id("EX_glc__D_e").lower_bound = -5.0
        return tmp.slim_optimize()


def test_prom_keeps_caller_constraints_in_ko_solves(re_model, low_glucose_growth):
    from mewpy.germ.analysis import PROM

    # probability 1 everywhere: the knockout changes nothing, so the KO growth is the constrained wild type
    probabilities = {key: 1.0 for key in _all_zero_probabilities(re_model)}
    solution = (
        PROM(re_model)
        .build()
        .optimize(
            initial_state=probabilities,
            regulators=["b4401"],
            solver_kwargs={"constraints": {"R_EX_glc__D_e": (-5.0, 1000.0)}},
        )
    )
    growth = next(iter(solution.solutions.values())).objective_value
    assert growth == pytest.approx(low_glucose_growth, abs=TOL)


def test_coregflux_keeps_caller_constraints(re_model, low_glucose_growth):
    from mewpy.germ.analysis import CoRegFlux

    solution = (
        CoRegFlux(re_model)
        .build()
        .optimize(initial_state={}, solver_kwargs={"constraints": {"R_EX_glc__D_e": (-5.0, 1000.0)}})
    )
    assert solution.objective_value == pytest.approx(low_glucose_growth, abs=TOL)


def test_rfba_gene_knockout_in_initial_state_is_applied(legacy_metabolic):
    from mewpy.germ.analysis import isingle_gene_deletion

    result = isingle_gene_deletion(legacy_metabolic, genes=["b1779"], method="rfba")
    assert result.loc["b1779", "growth"] == pytest.approx(0.0, abs=TOL)


def test_prom_without_regulators_knocks_out_every_regulator(re_model):
    from mewpy.germ.analysis import PROM

    solution = PROM(re_model).build().optimize(initial_state=_all_zero_probabilities(re_model))
    assert len(solution.solutions) == len(re_model.regulators)


def test_build_metabolites_legacy_model(legacy_metabolic):
    """Legacy stoichiometries are keyed by Metabolite objects; concentrations were silently dropped."""
    from mewpy.germ.analysis.analysis_utils import build_metabolites

    metabolites = build_metabolites(legacy_metabolic, {"glc__D_e": 10.0})
    assert metabolites["glc__D_e"].exchange == "EX_glc__D_e"


def test_pfba_infeasible_returns_status(legacy_metabolic):
    """An infeasible pFBA returns a non-optimal Solution, like the other methods, instead of raising."""
    from mewpy.germ.analysis.pfba import _pFBA
    from mewpy.solvers.solution import Status

    # biomass forced above its maximum (0.874)
    solution = _pFBA(legacy_metabolic).optimize(solver_kwargs={"constraints": {"Biomass_Ecoli_core": (1.0, 1000.0)}})
    assert solution.status in (Status.INFEASIBLE, Status.INF_OR_UNB)
    assert solution.fobj is None


@pytest.mark.parametrize("minimize", [False, True])
def test_pfba_honours_requested_objective(legacy_metabolic, cobra_model, minimize):
    """The first (FBA) step of pFBA must optimize the requested objective and direction, not model.objective."""
    from mewpy.germ.analysis.pfba import _pFBA

    with cobra_model as tmp:
        tmp.objective = "PGI"
        tmp.objective_direction = "min" if minimize else "max"
        expected = tmp.slim_optimize()

    from mewpy.util.constants import ModelConstants

    solution = _pFBA(legacy_metabolic).optimize(solver_kwargs={"linear": {"PGI": 1.0}, "minimize": minimize})
    # pFBA keeps the objective within the shared fraction (PFBA_OBJ_FRAC) of its optimum
    margin = (1 - getattr(ModelConstants, "PFBA_OBJ_FRAC", 0.999)) * abs(expected)
    assert solution.values["PGI"] == pytest.approx(expected + margin if minimize else expected - margin, abs=TOL)


# ----------------------------------------------------------------------------------------------------------------------
# attach=True: attached analyses are notified of model changes and rebuild on the next optimize
# ----------------------------------------------------------------------------------------------------------------------
def _fresh_legacy_integrated():
    from mewpy.io import Engines, Reader, read_model

    model = read_model(
        Reader(Engines.BooleanRegulatoryCSV, EC_CORE_TRN, **TRN_KWARGS), Reader(Engines.MetabolicSBML, EC_CORE)
    )
    model.objective = {"Biomass_Ecoli_core": 1}
    return model


def _optimize(analysis):
    from mewpy.germ.analysis import PROM

    if isinstance(analysis, PROM):
        # probability 1: the regulator knockout changes nothing, only the model change matters
        probabilities = {key: 1.0 for key in _all_zero_probabilities(analysis.model)}
        solution = analysis.optimize(initial_state=probabilities, regulators=["b4401"])
        return next(iter(solution.solutions.values())).objective_value or 0.0
    return analysis.optimize().objective_value or 0.0


@pytest.mark.parametrize("method", ["_FBA", "_pFBA", "RFBA", "SRFBA", "PROM", "CoRegFlux"])
def test_attached_analysis_rebuilds_after_model_change(method, wt_growth):
    import mewpy.germ.analysis as analysis_module
    from mewpy.germ.analysis.fba import _FBA
    from mewpy.germ.analysis.pfba import _pFBA

    cls = {"_FBA": _FBA, "_pFBA": _pFBA}.get(method) or getattr(analysis_module, method)
    model = _fresh_legacy_integrated() if method in ("PROM", "SRFBA", "CoRegFlux") else _fresh_legacy_metabolic()
    attached = cls(model, attach=True).build()
    assert attached in model.simulators
    wild_type = _optimize(attached)

    with model:
        model.reactions["GAPD"].ko()
        knockout = _optimize(attached)
        assert knockout == pytest.approx(_optimize(cls(model).build()), abs=TOL)
        assert knockout < wild_type - 1e-3

    assert _optimize(attached) == pytest.approx(wild_type, abs=TOL)


def _fresh_legacy_metabolic():
    from mewpy.io import Engines, Reader, read_model

    model = read_model(Reader(Engines.MetabolicSBML, EC_CORE))
    model.objective = {"Biomass_Ecoli_core": 1}
    return model


def test_attached_srfba_sees_regulator_knockout_regulatory_extension():
    from mewpy.germ.analysis import SRFBA

    model = _regulatory_extension(EC_CORE_TRN)
    attached = SRFBA(model, attach=True).build()
    wild_type = attached.optimize().objective_value
    # b3261 (Fis) deletion lowers SRFBA growth to ~0.21 (see the regulator deletion analysis)
    model.get_regulator("b3261").ko()
    fresh = SRFBA(model).build().optimize().objective_value or 0.0
    assert fresh < wild_type - 1e-3
    assert (attached.optimize().objective_value or 0.0) == pytest.approx(fresh, abs=TOL)


# ----------------------------------------------------------------------------------------------------------------------
# RFBA wild type (docs/germ.md setup): regulators default to the medium, TF targets are resolved before the genes
# ----------------------------------------------------------------------------------------------------------------------
EXAMPLES = Path(__file__).parent.parent.joinpath("examples", "models", "germ")
# documented RFBA wild-type growth for e_coli_core + e_coli_core_trn (docs/germ.md)
RFBA_DOCS_WT = 0.8739215069684303


def test_rfba_wild_type_matches_documented_growth_legacy():
    from mewpy.germ.analysis import RFBA
    from mewpy.io import Engines, Reader, read_model

    gem = Reader(Engines.MetabolicSBML, str(EXAMPLES.joinpath("e_coli_core.xml")))
    trn = Reader(Engines.BooleanRegulatoryCSV, str(EXAMPLES.joinpath("e_coli_core_trn.csv")), **TRN_KWARGS)
    model = read_model(gem, trn)
    rfba = RFBA(model).build()
    state = rfba.initial_state()
    # aerobic glucose medium: O2 and glucose present, acetate and glutamate absent
    assert state["o2_e"] > 0 and state["glc__D_e"] > 0
    assert state["ac_e"] == 0 and state["glu__L_e"] == 0
    assert rfba.optimize().objective_value == pytest.approx(RFBA_DOCS_WT, abs=TOL)


def test_rfba_wild_type_matches_documented_growth_regulatory_extension():
    from mewpy.germ.analysis import RFBA
    from mewpy.germ.models import RegulatoryExtension

    model = RegulatoryExtension.from_sbml(
        str(EXAMPLES.joinpath("e_coli_core.xml")),
        str(EXAMPLES.joinpath("e_coli_core_trn.csv")),
        regulatory_format="csv",
        **TRN_KWARGS,
    )
    assert RFBA(model).build().optimize().objective_value == pytest.approx(RFBA_DOCS_WT, abs=TOL)


@pytest.mark.parametrize("rule", ["b2415", "b2415 OR surplusFDP"])
def test_srfba_self_regulated_target_is_unconstrained(tmp_path, wt_growth, rule):
    """
    A rule 'x = x' (autoregulation) is a tautology; it must not force x to 0 (pre-fix: SRFBA infeasible).
    'x = x OR y' links x to a separate OR variable, so the target and its rule result never share a key.
    """
    from mewpy.germ.analysis import SRFBA

    rows = list(csv.reader(open(EC_CORE_TRN)))
    path = tmp_path.joinpath("trn_self_loop.csv")
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(rows[0])
        # b2415 (ptsH, glucose PTS) regulates itself
        writer.writerows([row[0], row[1], rule] if row[0] == "b2415" else row for row in rows[1:])
    model = _regulatory_extension(path)
    assert SRFBA(model).build().optimize().objective_value == pytest.approx(wt_growth, abs=TOL)


@pytest.mark.parametrize("method", ["_FBA", "_pFBA", "RFBA", "SRFBA", "CoRegFlux"])
def test_legacy_gene_knockout_closes_its_reactions(method):
    """gene.ko() sets the coefficients to 0; the analysis must close the reactions whose GPR is then off."""
    import mewpy.germ.analysis as analysis_module
    from mewpy.germ.analysis.fba import _FBA
    from mewpy.germ.analysis.pfba import _pFBA

    cls = {"_FBA": _FBA, "_pFBA": _pFBA}.get(method) or getattr(analysis_module, method)
    model = _fresh_legacy_metabolic()
    with model:
        # b1779 is the only gene of GAPD (essential)
        model.genes["b1779"].ko()
        assert model.reactions["GAPD"].bounds != (0, 0)
        solution = cls(model).build().optimize()
        assert (solution.values or {}).get("Biomass_Ecoli_core", 0.0) == pytest.approx(0.0, abs=TOL)


# ======================================================================================================================
# Round 3: SRFBA ko() semantics, metabolite regulators, tolerances, solver argument, ordering, user ids, pFBA fraction
# ======================================================================================================================
def _fba_growth(model):
    from mewpy.germ.analysis.fba import _FBA

    return _FBA(model).build().optimize().objective_value or 0.0


def _srfba(model, **kwargs):
    from mewpy.germ.analysis import SRFBA

    return SRFBA(model).build().optimize(to_solver=True, **kwargs)


@pytest.mark.parametrize("variable", ["b0721", "b1779", "b0113", "b3261"])
def test_srfba_ko_overrides_the_regulatory_rule_like_initial_state(variable):
    """USER DECISION: ko() on a TRN target (or regulator) overrides its rule exactly like initial_state={x: 0}."""
    expected = _srfba(_fresh_legacy_integrated(), initial_state={variable: 0.0})
    model = _fresh_legacy_integrated()
    with model:
        model.get(variable).ko()
        knockout = _srfba(model)
        assert knockout.status == expected.status
        assert (knockout.fobj or 0.0) == pytest.approx(expected.fobj or 0.0, abs=TOL)
        if variable == "b0721":  # sdhC, a target gene: same as FBA with the gene knocked out (0.8143)
            assert knockout.fobj == pytest.approx(_fba_growth(model), abs=TOL)


def test_srfba_metabolite_regulators_follow_their_exchange_flux():
    """USER DECISION: a metabolite regulator is 1 iff its exchange carries uptake beyond the tolerance."""
    model = _fresh_legacy_integrated()
    wild_type = _srfba(model)
    assert wild_type.values["EX_glc__D_e"] < -1 and wild_type.values["glc__D_e"] == pytest.approx(1.0)
    assert wild_type.values["EX_o2_e"] < -1 and wild_type.values["o2_e"] == pytest.approx(1.0)
    # acetate is not in the medium: no uptake, so the regulator is off
    assert wild_type.values["EX_ac_e"] >= -1e-6 and wild_type.values["ac_e"] == pytest.approx(0.0)

    anaerobic = _srfba(model, initial_state={"o2_e": 0.0})
    assert anaerobic.values["EX_o2_e"] == pytest.approx(0.0, abs=1e-5)


def test_srfba_indicator_tolerance_is_above_the_solver_feasibility_tolerance(re_model):
    from mewpy.germ.analysis import SRFBA

    srfba = SRFBA(re_model).build()
    # rows of the 'x > 0' indicators: a*(c + tol - ub) + x <= c + tol, with c = 0 for every '> 0' in the TRN
    rows = {name: data for name, data in srfba.solver._constr_data.items() if name.startswith("rel_")}
    upper = [rhs for name, (_, sense, rhs) in rows.items() if name.endswith("_1") and ">0)" in name and sense == "<"]
    assert upper and min(upper) >= 1e-6


def test_solver_argument_is_used(legacy_metabolic):
    from mewpy.germ.analysis.fba import _FBA
    from mewpy.solvers.optlang_solver import OptLangSolver

    fba = _FBA(legacy_metabolic, solver="optlang").build()
    assert isinstance(fba.solver, OptLangSolver)
    assert fba.optimize().objective_value == pytest.approx(0.8739215, abs=TOL)


def test_find_conflicts_order_is_deterministic():
    from mewpy.germ.analysis import find_conflicts
    from mewpy.io import Engines, Reader, read_model

    model = read_model(
        Reader(Engines.BooleanRegulatoryCSV, DATA.joinpath("SampleRegNet.csv"), sep=",", id_col=0, rule_col=1),
        Reader(Engines.MetabolicSBML, DATA.joinpath("SampleNet.xml")),
    )
    model.objective = {"r11": 1}
    model.get("pH").coefficients = (0, 14)
    genes, reactions = find_conflicts(model, initial_state={"g21": 0, "g23": 0})
    assert list(genes.index) == sorted(genes.index) and len(genes) > 1
    assert list(reactions.index) == sorted(reactions.index) and len(reactions) > 1
    assert list(genes.columns) == ["interaction"] + sorted(c for c in genes.columns if c != "interaction")
    assert list(reactions.columns) == ["gpr"] + sorted(c for c in reactions.columns if c != "gpr")


@pytest.fixture(scope="module")
def re_cobra():
    from mewpy.germ.models import RegulatoryExtension

    return RegulatoryExtension.from_sbml(
        str(EC_CORE), str(EC_CORE_TRN), regulatory_format="csv", flavor="cobra", **TRN_KWARGS
    )


@pytest.mark.parametrize("flavor, form", [("reframed", "unprefixed"), ("cobra", "prefixed")])
def test_user_ids_in_the_other_form_are_resolved(re_model, re_cobra, cobra_fva, flavor, form):
    from mewpy.germ.analysis import fva, slim_rfba, slim_srfba

    model = re_model if flavor == "reframed" else re_cobra

    def given(base, prefix):
        return prefix + base if form == "prefixed" else base

    assert slim_srfba(model, initial_state={given("b1779", "G_"): 0.0}) == pytest.approx(0.0, abs=TOL)
    table = fva(model, fraction=0.9, reactions=[given("PGI", "R_")])
    assert table.iloc[0]["minimum"] == pytest.approx(cobra_fva.loc["PGI", "minimum"], abs=TOL)
    assert table.iloc[0]["maximum"] == pytest.approx(cobra_fva.loc["PGI", "maximum"], abs=TOL)
    native = slim_srfba(model, constraints={("R_EX_o2_e" if flavor == "reframed" else "EX_o2_e"): (0, 1000)})
    assert slim_srfba(model, constraints={given("EX_o2_e", "R_"): (0, 1000)}) == pytest.approx(native, abs=TOL)
    assert native < 0.5  # anaerobic
    # 'o2_e' is the TRN's id of the oxygen regulator
    native_o2 = "M_o2_e" if flavor == "reframed" else "o2_e"
    expected = slim_rfba(model, initial_state={native_o2: 0.0})
    assert slim_rfba(model, initial_state={given("o2_e", "M_"): 0.0}) == pytest.approx(expected, abs=TOL)


def test_pfba_keeps_the_shared_fraction_of_the_optimum(legacy_metabolic, wt_growth):
    from mewpy.germ.analysis.pfba import _pFBA
    from mewpy.util.constants import ModelConstants

    fraction = getattr(ModelConstants, "PFBA_OBJ_FRAC", 0.999)
    solution = _pFBA(legacy_metabolic).optimize()
    assert solution.values["Biomass_Ecoli_core"] == pytest.approx(fraction * wt_growth, abs=1e-6)
    exact = _pFBA(legacy_metabolic).optimize(fraction=1.0)
    assert exact.values["Biomass_Ecoli_core"] == pytest.approx(wt_growth, abs=1e-6)


def test_prom_max_rates_only_for_regulated_reactions(re_model):
    """slim_prom/PROM with one regulator: FVA only over the reactions of that regulator's targets."""
    from mewpy.germ.analysis import PROM

    prom = PROM(re_model).build()
    calls = []
    original = prom._max_rates

    def spy(solver_kwargs, *args, **kwargs):
        rates = original(solver_kwargs, *args, **kwargs)
        calls.append(len(rates))
        return rates

    prom._max_rates = spy
    solution = prom.optimize(initial_state=_all_zero_probabilities(re_model), regulators=["b4401"])
    assert calls and calls[0] < len(re_model.reactions) / 2
    assert next(iter(solution.solutions.values())).objective_value < 0.874 - 1e-3
