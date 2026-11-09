"""Regression tests for omics, community, util and visualization defects (2026-11-09 review).

Each test reproduces one finding; the finding id is in the test name or docstring.
"""

import math
import subprocess
import sys
import warnings
from copy import deepcopy
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from mewpy.simulation import SStatus, get_simulator

EC_CORE_MODEL = str(Path(__file__).parent.joinpath("data", "e_coli_core.xml.gz"))


@lru_cache(maxsize=None)
def _cobra_core():
    from cobra.io import read_sbml_model

    return read_sbml_model(EC_CORE_MODEL)


@lru_cache(maxsize=None)
def _reframed_core():
    from reframed import load_cbmodel

    return load_cbmodel(EC_CORE_MODEL)


def cobra_model():
    return _cobra_core().copy()


def reframed_model():
    return deepcopy(_reframed_core())


def random_expression(sim, seed=0, scale=1.0):
    from mewpy.omics import ExpressionSet

    genes = list(sim.genes)
    rng = np.random.default_rng(seed)
    return ExpressionSet(genes, ["c1"], rng.random((len(genes), 1)) * scale)


def two_member_community(loader=reframed_model, knockouts=(), **kwargs):
    """Two copies of E. coli core with ATPM off; ``knockouts`` are (organism, reaction) pairs."""
    from mewpy.com import CommunityModel

    models = []
    for org_id in ("m1", "m2"):
        model = loader()
        sim = get_simulator(model)
        prefix = sim._r_prefix
        sim.set_reaction_bounds(prefix + "ATPM", 0, 0)
        for org, r_id in knockouts:
            if org == org_id:
                sim.set_reaction_bounds(prefix + r_id, 0, 0)
        model.id = org_id
        models.append(model)
    return CommunityModel(models, verbose=False, **kwargs)


# ---------------------------------------------------------------------------------------------------------------------
# F1: iMAT big-M constraints
# ---------------------------------------------------------------------------------------------------------------------
def test_f1_imat_full_expression_is_feasible():
    from mewpy.omics import iMAT

    sim = get_simulator(cobra_model())
    res = iMAT(sim, random_expression(sim))
    assert res.status == SStatus.OPTIMAL


def test_f1_imat_reversible_high_reaction_is_active():
    from mewpy.omics import iMAT

    sim = get_simulator(cobra_model())
    res = iMAT(sim, ({}, {"PGI": 1.0}))
    assert res.status == SStatus.OPTIMAL
    assert abs(res.fluxes["PGI"]) >= 1 - 1e-6


def test_f1_imat_low_reaction_is_not_forced_inactive():
    from mewpy.omics import iMAT

    sim = get_simulator(cobra_model())
    res = iMAT(sim, ({"ATPS4r": 1.0}, {}), constraints={"ATPS4r": (40, 1000)})
    assert res.status == SStatus.OPTIMAL
    assert res.fluxes["ATPS4r"] >= 40 - 1e-6


def test_f1_imat_handles_infinite_bounds():
    from mewpy.omics import iMAT

    sim = get_simulator(reframed_model())
    assert any(math.isinf(sim.get_reaction_bounds(r)[1]) for r in sim.reactions)
    res = iMAT(sim, random_expression(sim))
    assert res.status == SStatus.OPTIMAL


# ---------------------------------------------------------------------------------------------------------------------
# F2: SteadyCom must disable the merged biomass equation
# ---------------------------------------------------------------------------------------------------------------------
def test_f2_steadycom_default_community_abundances_sum_to_one():
    from mewpy.com.steadycom import SteadyCom

    sol = SteadyCom(two_member_community())
    ref = SteadyCom(two_member_community(merge_biomasses=False))
    assert sum(sol.abundance.values()) == pytest.approx(1.0, abs=1e-6)
    assert sol.growth == pytest.approx(ref.growth, rel=1e-3)


# ---------------------------------------------------------------------------------------------------------------------
# F3: GIMME build_model=True
# ---------------------------------------------------------------------------------------------------------------------
def test_f3_gimme_build_model_returns_context_specific_model():
    from mewpy.omics import GIMME

    model = cobra_model()
    sim = get_simulator(model)
    n_before = len(sim.reactions)
    wt = sim.simulate().objective_value
    biomass = sim.biomass_reaction
    expr = random_expression(sim, scale=100)

    res, new_sim = GIMME(model, expr, build_model=True)

    assert res.status == SStatus.OPTIMAL
    assert len(get_simulator(model).reactions) == n_before
    assert res.fluxes[biomass] >= 0.9 * wt - 1e-6
    assert biomass in new_sim.reactions
    assert new_sim.simulate().objective_value >= 0.9 * wt - 1e-6
    n_irreversible = n_before + sum(1 for r in sim.reactions if sim.get_reaction_bounds(r)[0] < 0)
    assert len(new_sim.reactions) < n_irreversible


# ---------------------------------------------------------------------------------------------------------------------
# F4: NaN and negative expression values
# ---------------------------------------------------------------------------------------------------------------------
def _nan_expression(sim):
    from mewpy.omics import ExpressionSet

    genes = list(sim.genes)
    single = [sim.get_gpr(r).strip("() ") for r in sim.reactions if sim.get_gpr(r)]
    single = [g for g in single if " " not in g]
    rng = np.random.default_rng(0)
    values = rng.random((len(genes), 1))
    values[genes.index(single[0]), 0] = np.nan
    return ExpressionSet(genes, ["c1"], values)


def test_f4_reaction_expression_ignores_nan_genes():
    from mewpy.omics import Preprocessing

    sim = get_simulator(cobra_model())
    rxn_exp = Preprocessing(sim, _nan_expression(sim)).reactions_expression(0)
    assert rxn_exp
    assert not any(math.isnan(v) for v in rxn_exp.values())


def test_f4_nan_and_operand_order_does_not_matter():
    from mewpy.omics import gene_to_reaction_expression

    sim = get_simulator(cobra_model())
    gpr = sim.get_gpr("SUCDi")
    genes = sorted(set(gpr.replace("(", " ").replace(")", " ").split()) - {"and", "or"})
    for nan_gene in genes:
        values = {g: 5.0 for g in genes}
        values[nan_gene] = float("nan")
        assert gene_to_reaction_expression(sim, values)["SUCDi"] == 5.0


def test_f4_percentile_threshold_is_finite_with_nan():
    from mewpy.omics import Preprocessing

    sim = get_simulator(cobra_model())
    coeffs, threshold = Preprocessing(sim, _nan_expression(sim)).percentile(0, cutoff=25)
    assert math.isfinite(threshold)
    assert len(coeffs) > 0


def test_f4_eflux_with_nan_expression_is_feasible():
    from mewpy.omics import eFlux

    sim = get_simulator(cobra_model())
    assert eFlux(sim, _nan_expression(sim)).status == SStatus.OPTIMAL


def test_f4_eflux_rejects_negative_expression():
    from mewpy.omics import eFlux

    sim = get_simulator(cobra_model())
    expr = random_expression(sim, seed=1, scale=4)
    expr.apply()  # log2, so values below 1 become negative
    assert expr.get_condition(0).min() < 0
    with pytest.raises(ValueError, match="non-negative"):
        eFlux(sim, expr)


def test_f4_quantile_pipeline_imputes_nan_by_default():
    pytest.importorskip("sklearn")
    from mewpy.omics import ExpressionSet

    rows = [[1.0, 2.0, 3.0], [np.nan, 2.5, 3.1], [0.5, 0.7, 0.9]]
    rows += [[2.0, 3.0, 4.0], [1.1, 1.2, 1.3], [5.0, 6.0, 7.0]]
    values = np.array(rows)
    expr = ExpressionSet([f"g{i}" for i in range(6)], ["a", "b", "c"], values)
    quantile, binary = expr.quantile_pipeline()
    assert not quantile.isna().any().any()
    assert not binary.isna().any().any()


# ---------------------------------------------------------------------------------------------------------------------
# F6: SMETANA MIP score
# ---------------------------------------------------------------------------------------------------------------------
def _core_medium():
    from mewpy.simulation import Environment

    return Environment.from_model(reframed_model()).simplify()


def test_f6_mip_noninteracting_community_cannot_cross_feed():
    """m1 cannot take up glucose: it grows only on m2 secretions, which non-interaction forbids."""
    from mewpy.com.analysis import mip_score

    community = two_member_community(knockouts=[("m1", "GLCpts")])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        score, extras = mip_score(community, environment=_core_medium())
    assert score is None and extras is None


def test_f6_mip_counts_cross_feeding():
    """m1 lacks glucose uptake, m2 fructose uptake: only non-interacting organisms need fructose."""
    from mewpy.com.analysis import mip_score

    env = _core_medium()
    env["R_EX_fru_e"] = (-10, math.inf)
    knockouts = [("m1", "GLCpts"), ("m2", "FRUpts2")]
    score, extras = mip_score(two_member_community(knockouts=knockouts), environment=env)
    assert score == 1
    assert "fru" in extras["noninteracting_medium"]
    assert "fru" not in extras["interacting_medium"]
    score, _ = mip_score(two_member_community(), environment=env)
    assert score == 0


def test_f6_mip_exclude_removes_compounds():
    from mewpy.com.analysis import mip_score

    community = two_member_community()
    score, extras = mip_score(community)
    assert "pi" in extras["interacting_medium"]
    score, extras = mip_score(community, exclude={"pi"})
    assert "pi" not in extras["interacting_medium"]
    assert "pi" not in extras["noninteracting_medium"]


# ---------------------------------------------------------------------------------------------------------------------
# F7: SMETANA MP score medium exclusion and metabolite trimming
# ---------------------------------------------------------------------------------------------------------------------
def test_f7_mp_score_excludes_environment_compounds():
    from mewpy.com.analysis import mp_score
    from mewpy.simulation import Environment

    community = two_member_community()
    env = Environment.from_model(reframed_model())
    assert env["R_EX_glc__D_e"][0] < 0
    scores = mp_score(community, environment=env)
    assert "M_glc__D_e" not in scores["m1"]
    assert "M_ac_e" in scores["m1"]


def test_f7_mro_names_compounds_by_shared_metabolite():
    """Guards the F7 trim change: organism-specific ids must not leak into MRO media."""
    from mewpy.com.analysis import mro_score

    score, extras = mro_score(two_member_community(), environment=_core_medium())
    media = extras["individual_media"]
    assert media["m1"] and media["m2"]
    assert not any(c.endswith(("_m1", "_m2")) for medium in media.values() for c in medium)
    assert "pi" in media["m1"] & media["m2"]
    assert score > 0


def test_f7_trim_keeps_stereoisomers_apart():
    from mewpy.com.analysis import _trim_metabolite_prefix

    sim = two_member_community().merged_model
    assert _trim_metabolite_prefix(sim, "M_lac__D_e") == "lac__D"
    assert _trim_metabolite_prefix(sim, "M_lac__L_e") == "lac__L"


# ---------------------------------------------------------------------------------------------------------------------
# F8: CommunitySolution exchanges and mass flow
# ---------------------------------------------------------------------------------------------------------------------
def _assert_exchanges_balance(sol):
    sim = sol.community.merged_model
    exchanges = sol.compute_exchanges()
    assert exchanges
    for m_id in sol.community.ext_mets:
        total = sum(rate for (_, met), rate in exchanges.items() if met == m_id)
        for r_id in set(sim.get_metabolite_reactions(m_id)) & set(sim.get_exchange_reactions()):
            total += sim.get_reaction(r_id).stoichiometry[m_id] * sol.values[r_id]
        assert total == pytest.approx(0, abs=1e-6), m_id


def test_f8_compute_exchanges_for_cobra_organisms():
    from mewpy.com.steadycom import SteadyCom

    sol = SteadyCom(two_member_community(loader=cobra_model, merge_biomasses=False))
    _assert_exchanges_balance(sol)
    assert len(sol.cross_feeding()) > 0


def test_f8_compute_exchanges_for_reframed_organisms():
    from mewpy.com.steadycom import SteadyCom

    sol = SteadyCom(two_member_community(merge_biomasses=False))
    _assert_exchanges_balance(sol)


def test_f8_mass_flow_runs():
    from mewpy.com.steadycom import SteadyCom

    sol = SteadyCom(two_member_community(merge_biomasses=False))
    flow = sol.mass_flow(element="C")
    assert isinstance(flow, dict)
    df = sol.mass_flow(as_df=True)
    assert list(df.columns) == ["donor", "receiver", "flow"]


# ---------------------------------------------------------------------------------------------------------------------
# F9: importing mewpy must not reconfigure the host application's logging
# ---------------------------------------------------------------------------------------------------------------------
def _run_python(code):
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)


def test_f9_import_keeps_host_logging():
    code = (
        "import io, logging\n"
        "buf = io.StringIO()\n"
        "logging.basicConfig(stream=buf, level=logging.WARNING)\n"
        "host = logging.StreamHandler(buf)\n"
        "logging.getLogger('mewpy').addHandler(host)\n"
        "import mewpy\n"
        "lg = logging.getLogger('mewpy')\n"
        "assert host in lg.handlers, 'host handler removed'\n"
        "assert lg.propagate, 'propagation disabled'\n"
        "logging.getLogger('mewpy.simulation').warning('library-warning')\n"
        "assert 'library-warning' in buf.getvalue(), 'host did not receive the record'\n"
    )
    proc = _run_python(code)
    assert proc.returncode == 0, proc.stderr[-2000:]


def test_f9_console_logging_is_opt_in():
    code = (
        "import logging, mewpy\n"
        "logging.getLogger('mewpy.x').info('quiet-by-default')\n"
        "mewpy.configure_logging('INFO')\n"
        "logging.getLogger('mewpy.x').info('now-visible')\n"
    )
    proc = _run_python(code)
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "quiet-by-default" not in proc.stderr
    assert "now-visible" in proc.stderr


# ---------------------------------------------------------------------------------------------------------------------
# F12 and F-VIZ-2: visualization
# ---------------------------------------------------------------------------------------------------------------------
def test_f12_flux_envelope_flip_x_keeps_pairs():
    import matplotlib

    matplotlib.use("Agg", force=True)
    from mewpy.visualization.envelope import flux_envelope, plot_flux_envelope

    model = reframed_model()
    r_x, r_y = "R_EX_glc__D_e", "R_BIOMASS_Ecoli_core_w_GAM"
    x, _, ymax = flux_envelope(model, r_x, r_y, steps=3)
    ax = plot_flux_envelope(model, r_x, r_y, steps=3, flip_x=True)
    xs, ys = ax.lines[1].get_data()
    expected = np.array(sorted(zip((-x).tolist(), ymax.tolist())))
    plotted = np.array(sorted(zip(list(xs), list(ys))))
    np.testing.assert_allclose(plotted, expected, atol=1e-9)


def test_fviz2_pause_works_on_current_matplotlib():
    import matplotlib

    matplotlib.use("Agg", force=True)
    from mewpy.visualization.plot import pause

    pause(0.001)


# ---------------------------------------------------------------------------------------------------------------------
# Lower-severity items
# ---------------------------------------------------------------------------------------------------------------------
def test_differences_uses_consecutive_pair_p_values():
    from mewpy.omics import ExpressionSet

    values = np.array([[1.0, 2.0, 1.0]])
    # p-value columns are the pairs (t0,t1), (t0,t2), (t1,t2)
    expr = ExpressionSet(["g"], ["t0", "t1", "t2"], values, p_values=np.array([[0.9, 0.001, 0.9]]))
    assert expr.differences() == {"g": [0, 0]}
    expr = ExpressionSet(["g"], ["t0", "t1", "t2"], values, p_values=np.array([[0.9, 0.9, 0.001]]))
    assert expr.differences() == {"g": [0, -1]}


def test_from_dataframe_does_not_mutate_input():
    from mewpy.omics import ExpressionSet

    df = pd.DataFrame({0: [1.0, 2.0, 3.0], 1: [4.0, 5.0, 6.0]}, index=["g1", "g1", "g2"])
    original = df.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ExpressionSet.from_dataframe(df)
    pd.testing.assert_frame_equal(df, original)


def test_from_dataframe_suffix_avoids_existing_ids():
    from mewpy.omics import ExpressionSet

    df = pd.DataFrame({"a": [1.0, 2.0, 3.0]}, index=["A", "A", "A_2"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expr = ExpressionSet.from_dataframe(df)
    ids = expr.dataframe.index.tolist()
    assert len(ids) == 3 and len(set(ids)) == 3
    assert ids[0] == "A" and ids[2] == "A_2"


def test_elements_accepts_fractional_counts():
    from mewpy.util.utilities import elements, molecular_weight

    assert elements("C1.5H3O2") == {"C": 1.5, "H": 3, "O": 2}
    assert elements("C6H12O6") == {"C": 6, "H": 12, "O": 6}
    assert molecular_weight("C1.5H3O2") > 0


def test_create_metabolic_graph_undirected():
    import networkx as nx

    from mewpy.util.graph import create_metabolic_graph

    model = reframed_model()
    directed = create_metabolic_graph(model)
    undirected = create_metabolic_graph(model, directed=False)
    assert isinstance(undirected, nx.Graph) and not undirected.is_directed()
    assert set(undirected.nodes) == set(directed.nodes)
    assert undirected.number_of_edges() > 0


def test_balance_exchanges_false_restores_stoichiometry():
    community = two_member_community(abundances=[0.3, 0.7], balance_exchange=True)
    sim = community.merged_model
    trp = community.reaction_map[("m1", "R_EX_glc__D_e")]
    assert sim.get_reaction(trp).stoichiometry["M_glc__D_e"] == pytest.approx(0.3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        community.balance_exchanges = False
    sim = community.merged_model
    trp = community.reaction_map[("m1", "R_EX_glc__D_e")]
    assert sim.get_reaction(trp).stoichiometry["M_glc__D_e"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------------------------------------------------
# Built context-specific models (upstream issues #45 and #48), on iMM904 with real expression data
# ---------------------------------------------------------------------------------------------------------------------
IMM904_DIR = Path(__file__).parent.parent.joinpath("examples", "models", "germ")


@lru_cache(maxsize=None)
def _imm904_data():
    from cobra.io import read_sbml_model

    from mewpy.omics import ExpressionSet

    df = pd.read_csv(IMM904_DIR.joinpath("iMM904_gene_expression.csv"), sep=";", index_col=0)
    condition = df.columns[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expr = ExpressionSet.from_dataframe(df[[condition]])
    return read_sbml_model(str(IMM904_DIR.joinpath("iMM904.xml"))), expr, condition


def _imm904():
    model, expr, condition = _imm904_data()
    return get_simulator(model.copy()), expr, condition


def _assert_built_model_grows_and_saves(built, tmp_path):
    from cobra.io import write_sbml_model

    sol = built.simulate()
    assert sol.status == SStatus.OPTIMAL
    assert sol.objective_value > 1e-6
    out = tmp_path / "built.xml"
    write_sbml_model(built.model, str(out))
    assert out.stat().st_size > 0
    return sol


requires_imm904 = pytest.mark.skipif(not IMM904_DIR.joinpath("iMM904.xml").exists(), reason="iMM904 not available")


@requires_imm904
def test_imat_built_model_grows_and_saves(tmp_path):
    from mewpy.omics import iMAT

    sim, expr, condition = _imm904()
    n_before = len(sim.reactions)
    res, built = iMAT(sim, expr, condition=condition, build_model=True)
    assert res.status == SStatus.OPTIMAL
    assert len(sim.reactions) == n_before
    assert sim.biomass_reaction in built.reactions
    assert len(built.reactions) < n_before
    _assert_built_model_grows_and_saves(built, tmp_path)


@requires_imm904
def test_imat_reports_milp_objective():
    from mewpy.omics import iMAT

    sim, expr, condition = _imm904()
    res = iMAT(sim, expr, condition=condition)
    assert res.status == SStatus.OPTIMAL
    assert res.objective_value is not None and res.objective_value > 0


@requires_imm904
def test_eflux_built_model_reproduces_solution(tmp_path):
    from mewpy.omics import eFlux

    sim, expr, condition = _imm904()
    fba = sim.simulate().objective_value
    res, built = eFlux(sim, expr, condition=condition, build_model=True)
    assert res.status == SStatus.OPTIMAL
    assert res.objective_value < 0.5 * fba
    sol = _assert_built_model_grows_and_saves(built, tmp_path)
    assert sol.objective_value == pytest.approx(res.objective_value, rel=1e-4)
    assert sim.simulate().objective_value == pytest.approx(fba, rel=1e-6)


@requires_imm904
def test_gimme_built_model_still_grows_and_saves(tmp_path):
    from mewpy.omics import GIMME

    sim, expr, condition = _imm904()
    res, built = GIMME(sim, expr, condition=condition, build_model=True)
    assert res.status == SStatus.OPTIMAL
    _assert_built_model_grows_and_saves(built, tmp_path)


# ---------------------------------------------------------------------------------------------------------------------
# Backend parity: expression ids without the SBML 'G_' prefix must match REFRAMED genes too
# ---------------------------------------------------------------------------------------------------------------------
def _core_expression_unprefixed():
    from mewpy.omics import ExpressionSet

    genes = sorted(get_simulator(cobra_model()).genes)  # 'b1241', as in typical expression files
    assert not any(g.startswith("G_") for g in genes)
    rng = np.random.default_rng(0)
    return ExpressionSet(genes, ["c1"], rng.random((len(genes), 1)) * 100)


_BACKENDS = {"cobra": cobra_model, "reframed": reframed_model}


@lru_cache(maxsize=None)
def _parity_run(backend, method):
    from mewpy.omics import GIMME, Preprocessing, eFlux, iMAT

    sim = get_simulator(_BACKENDS[backend]())
    expr = _core_expression_unprefixed()
    n_rxn_exp = len(Preprocessing(sim, expr).reactions_expression(0))
    run = {"eFlux": eFlux, "GIMME": GIMME, "iMAT": iMAT}[method]
    res = run(sim, expr)
    return n_rxn_exp, res.status, res.objective_value


@pytest.mark.parametrize("method", ["eFlux", "GIMME", "iMAT"])
@pytest.mark.parametrize("backend", ["cobra", "reframed"])
def test_backend_parity_omics(backend, method):
    n_cobra, status_cobra, obj_cobra = _parity_run("cobra", method)
    n, status, obj = _parity_run(backend, method)
    assert n_cobra > 0
    assert n == n_cobra
    assert status == status_cobra == SStatus.OPTIMAL
    assert obj == pytest.approx(obj_cobra, rel=1e-4, abs=1e-6)


@requires_imm904
def test_backend_parity_imm904_reaction_expression():
    from reframed import load_cbmodel

    from mewpy.omics import Preprocessing

    sim, expr, condition = _imm904()
    reframed_sim = get_simulator(load_cbmodel(str(IMM904_DIR.joinpath("iMM904.xml"))))
    n_cobra = len(Preprocessing(sim, expr).reactions_expression(condition))
    n_reframed = len(Preprocessing(reframed_sim, expr).reactions_expression(condition))
    assert n_cobra > 0
    assert n_reframed == n_cobra


@lru_cache(maxsize=None)
def _imm904_reframed_model():
    from reframed import load_cbmodel

    return load_cbmodel(str(IMM904_DIR.joinpath("iMM904.xml")))


@requires_imm904
@pytest.mark.parametrize("method", ["eFlux", "GIMME", "iMAT"])
def test_backend_parity_imm904_methods(method):
    """COBRApy reads iMM904's 'unbounded' reactions as +-999999 and REFRAMED as +-inf."""
    from mewpy.omics import GIMME, eFlux, iMAT

    run = {"eFlux": eFlux, "GIMME": GIMME, "iMAT": iMAT}[method]
    sim, expr, condition = _imm904()
    reframed_sim = get_simulator(deepcopy(_imm904_reframed_model()))
    res_cobra = run(sim, expr, condition=condition)
    res_reframed = run(reframed_sim, expr, condition=condition)
    assert res_cobra.status == res_reframed.status == SStatus.OPTIMAL
    assert res_reframed.objective_value == pytest.approx(res_cobra.objective_value, rel=1e-4, abs=1e-6)


# ---------------------------------------------------------------------------------------------------------------------
# Prefix handling through the Simulator API: mixed COBRApy/REFRAMED communities and model comparisons
# ---------------------------------------------------------------------------------------------------------------------
def _named(loader, org_id):
    model = loader()
    model.id = org_id
    return model


def _community_of(loaders, **kwargs):
    from mewpy.com import CommunityModel

    return CommunityModel([_named(loader, f"m{i + 1}") for i, loader in enumerate(loaders)], verbose=False, **kwargs)


def _community_ids(community):
    sim = community.merged_model
    return {kind: set(getattr(sim, kind)) for kind in ("reactions", "metabolites", "genes")}


@pytest.mark.parametrize("flavor", ["reframed", "cobra"])
@pytest.mark.parametrize("merge_biomasses", [True, False])
def test_mixed_backend_community_ids_match_same_backend(flavor, merge_biomasses):
    kwargs = {"flavor": flavor, "merge_biomasses": merge_biomasses}
    mixed = _community_ids(_community_of((cobra_model, reframed_model), **kwargs))
    all_cobra = _community_ids(_community_of((cobra_model, cobra_model), **kwargs))
    all_reframed = _community_ids(_community_of((reframed_model, reframed_model), **kwargs))
    assert mixed == all_cobra == all_reframed


@pytest.mark.parametrize("flavor", ["reframed", "cobra"])
def test_community_ids_unchanged(flavor):
    """Guard: ids built from organism ids keep today's form."""
    ids = _community_ids(_community_of((reframed_model, reframed_model), flavor=flavor))
    if flavor == "reframed":
        expected = {
            "reactions": {"R_PGI_m1", "R_EX_glc__D_e_m2", "R_EX_glc__D_e", "community_growth"},
            "metabolites": {"M_glc__D_e", "M_glc__D_e_m1", "M_g6p_c_m2"},
            "genes": {"G_b1241_m1"},
        }
    else:
        expected = {
            "reactions": {"PGI_m1", "EX_glc__D_e_m2", "EX_glc__D_e", "community_growth"},
            "metabolites": {"glc__D_e", "glc__D_e_m1", "g6p_c_m2"},
            "genes": {"b1241_m1"},
        }
    for kind, sample in expected.items():
        assert sample <= ids[kind], (kind, sample - ids[kind])


def _glucose_constraints(community):
    from mewpy.simulation import Environment

    sim = community.merged_model
    medium = Environment.from_model(reframed_model()).simplify()
    constraints = {r_id: (0, 1000) for r_id in sim.get_exchange_reactions()}
    constraints.update({sim.resolve(r_id, "reaction"): bounds for r_id, bounds in medium.items()})
    return constraints


def test_mixed_backend_steadycom_matches_same_backend():
    from mewpy.com.steadycom import SteadyCom

    solutions = {}
    for name, loaders in (("mixed", (cobra_model, reframed_model)), ("reframed", (reframed_model, reframed_model))):
        community = _community_of(loaders, merge_biomasses=False)
        solutions[name] = SteadyCom(community, constraints=_glucose_constraints(community))
    mixed, same = solutions["mixed"], solutions["reframed"]
    assert mixed.growth == pytest.approx(same.growth, rel=1e-4)
    assert sum(mixed.abundance.values()) == pytest.approx(1.0, abs=1e-6)
    _assert_exchanges_balance(mixed)
    assert {m for _, m in mixed.compute_exchanges()} <= set(same.community.ext_mets)
    assert len(mixed.cross_feeding()) > 0


def test_similarity_across_backends(tmp_path):
    from mewpy.com.similarity import (
        jaccard_similarity,
        resource_overlap,
        write_out_common_metabolites,
        write_out_common_reactions,
    )

    n_reactions = len(get_simulator(cobra_model()).reactions)
    n_metabolites = len(get_simulator(cobra_model()).metabolites)
    assert jaccard_similarity(cobra_model(), reframed_model()) == (1.0, 1.0)
    assert resource_overlap(cobra_model(), reframed_model()) == 1.0
    common_mets = write_out_common_metabolites([cobra_model(), reframed_model()], prefix=str(tmp_path / "m.csv"))
    assert len(common_mets) == n_metabolites
    common_rxns = write_out_common_reactions([cobra_model(), reframed_model()], prefix=str(tmp_path / "r.csv"))
    assert len(common_rxns) == n_reactions


def test_escher_reaction_data_uses_canonical_ids(monkeypatch, tmp_path):
    import types

    fake = types.ModuleType("escher")
    fake.Builder = lambda **kwargs: kwargs
    fake.list_available_maps = lambda: [{"map_name": "e_coli_core.Core metabolism"}]
    monkeypatch.setitem(sys.modules, "escher", fake)
    monkeypatch.chdir(tmp_path)
    from mewpy.visualization.escher import build_escher

    assert build_escher(model=reframed_model(), fluxes={"R_PGI": 1.0})["reaction_data"] == {"PGI": 1.0}
    assert build_escher(model=cobra_model(), fluxes={"PGI": 1.0})["reaction_data"] == {"PGI": 1.0}
    assert build_escher(fluxes={"R_PGI": 1.0})["reaction_data"] == {"PGI": 1.0}


def test_reaction_expression_exact_gene_id_wins():
    from mewpy.omics import gene_to_reaction_expression

    sim = get_simulator(reframed_model())
    gpr = sim.get_gpr("R_PGI")
    assert gpr.strip("() ") == "G_b4025"
    assert gene_to_reaction_expression(sim, {"b4025": 1.0, "G_b4025": 7.0})["R_PGI"] == 7.0
    assert gene_to_reaction_expression(sim, {"G_b4025": 7.0, "b4025": 1.0})["R_PGI"] == 7.0
    assert gene_to_reaction_expression(sim, {"b4025": 3.0})["R_PGI"] == 3.0
