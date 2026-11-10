"""The Simulator owns SBML id prefixes: ids are accepted with or without R_/M_/G_ on every backend.

REFRAMED keeps the SBML prefixes (R_PGI, M_glc__D_e, G_b4025); COBRApy and GERM drop them.
Return values stay backend-native.
"""

from pathlib import Path

import pytest

EC_CORE_XML = str(Path(__file__).parent.joinpath("data", "e_coli_core.xml"))
WT_GROWTH = 0.8739
PGI_KO_GROWTH = 0.8632
BACKENDS = ["cobra", "reframed"]


def load(backend, **kwargs):
    from mewpy.simulation import get_simulator

    if backend == "cobra":
        from cobra.io import read_sbml_model

        return get_simulator(read_sbml_model(EC_CORE_XML), **kwargs)
    from reframed.io.sbml import load_cbmodel

    return get_simulator(load_cbmodel(EC_CORE_XML), **kwargs)


def native(backend, base, kind="reaction"):
    prefixes = {"reaction": "R_", "metabolite": "M_", "gene": "G_"}
    return prefixes[kind] + base if backend == "reframed" else base


@pytest.fixture(params=BACKENDS)
def backend(request):
    return request.param


@pytest.fixture
def sim(backend):
    return load(backend)


def both_forms(base, kind="reaction"):
    prefixes = {"reaction": "R_", "metabolite": "M_", "gene": "G_"}
    return [base, prefixes[kind] + base]


def test_api(sim, backend):
    expected_prefix = {"reaction": "R_", "metabolite": "M_", "gene": "G_"} if backend == "reframed" else None
    for kind, base in (("reaction", "PGI"), ("metabolite", "glc__D_e"), ("gene", "b4025")):
        nid = native(backend, base, kind)
        assert sim.prefix(kind) == (expected_prefix[kind] if expected_prefix else "")
        assert [sim.resolve(i, kind) for i in both_forms(base, kind)] == [nid, nid]
        assert sim.resolve_many(both_forms(base, kind), kind) == [nid, nid]
        assert sim.to_native(base, kind) == nid
        assert sim.to_canonical(nid, kind) == base
    expected = {native(backend, "PGI"): 0, native(backend, "ACKr"): (0, 1)}
    assert sim.resolve_keys({"PGI": 0, "R_ACKr": (0, 1)}) == expected
    assert sim._r_prefix == sim.prefix("reaction")
    with pytest.raises(ValueError):
        sim.prefix("compartment")


def test_fba_objective_both_forms(sim, backend):
    for objective in both_forms("Biomass_Ecoli_core"):
        assert sim.simulate(objective={objective: 1}).objective_value == pytest.approx(WT_GROWTH, abs=1e-4)
        assert sim.simulate(objective=objective).objective_value == pytest.approx(WT_GROWTH, abs=1e-4)
    for objective in both_forms("EX_ac_e"):
        sim.set_objective(objective)
        assert sim.objective[native(backend, "EX_ac_e")] == 1.0
        sim.objective = {objective: 1}
        assert {k: v for k, v in sim.objective.items() if v} == {native(backend, "EX_ac_e"): 1.0}


def test_knockout_both_forms(sim, backend):
    # covers the unprefixed id on REFRAMED and the prefixed id on COBRApy
    for ko in both_forms("PGI"):
        res = sim.simulate(constraints={ko: 0})
        assert res.objective_value == pytest.approx(PGI_KO_GROWTH, abs=1e-4)
        assert res.fluxes[native(backend, "PGI")] == pytest.approx(0.0, abs=1e-9)
        assert set(res.simulation_constraints) == {native(backend, "PGI")}
        assert sim.is_essential_reaction(ko) is False
        assert sim.simulate(constraints={ko: 0}, slim=True) == pytest.approx(PGI_KO_GROWTH, abs=1e-4)
    for method in ("pFBA", "lMOMA"):
        values = [sim.simulate(method=method, constraints={ko: 0}).objective_value for ko in both_forms("PGI")]
        assert values[0] == pytest.approx(values[1], abs=1e-6)


def test_fva_both_forms(sim, backend):
    expected = {"PGI": (-14.2990, 9.8388), "ACKr": (-3.8136, 0.0)}
    for names in (["PGI", "ACKr"], ["R_PGI", "R_ACKr"]):
        ranges = sim.FVA(reactions=names, obj_frac=0.9)
        assert set(ranges) == {native(backend, "PGI"), native(backend, "ACKr")}
        for base, (lo, hi) in expected.items():
            assert ranges[native(backend, base)] == pytest.approx([lo, hi], abs=1e-3)
    for ko in both_forms("PGI"):
        ranges = sim.FVA(reactions="ACKr", constraints={ko: 0}, obj_frac=0.9)
        assert list(ranges) == [native(backend, "ACKr")]
    blocked = [sim.blocked_reactions(reactions=[r]) for r in both_forms("PGI")]
    assert blocked == [[], []]


def test_genes_both_forms(sim, backend):
    essential = sim.essential_genes()
    assert len(essential) == 7
    for gene in both_forms("b0720", "gene"):
        assert sim.is_essential_gene(gene) is True
        assert sim.get_gene(gene)["id"] == native(backend, "b0720", "gene")
    active = [g for g in sim.genes if not g.endswith("b4025")]
    for genes in ([sim.to_canonical(g, "gene") for g in active], ["G_" + sim.to_canonical(g, "gene") for g in active]):
        assert native(backend, "PGI") not in sim.evaluate_gprs(genes)


def test_bounds_round_trip(sim, backend):
    for set_id, get_id in (("PGI", "R_PGI"), ("R_PGI", "PGI")):
        sim.set_reaction_bounds(set_id, -5, 5)
        assert sim.get_reaction_bounds(get_id) == (-5, 5)
        assert sim.get_reaction_bounds(native(backend, "PGI")) == (-5, 5)
        sim.set_reaction_bounds(set_id, -1000, 1000)


def test_getters_both_forms(sim, backend):
    pgi, glc = native(backend, "PGI"), native(backend, "glc__D_e", "metabolite")
    for r in both_forms("PGI"):
        assert sim.get_reaction(r)["id"] == pgi
        assert sim.get_gpr(r) == native(backend, "b4025", "gene")
        assert set(sim.get_reaction_metabolites(r)) == set(sim.get_reaction_metabolites(pgi))
        assert sim.get_substrates(r) == sim.get_substrates(pgi)
        assert sim.get_products(r) == sim.get_products(pgi)
        assert sim.reverse_reaction(r) == pgi
    for m in both_forms("glc__D_e", "metabolite"):
        assert sim.get_metabolite(m)["id"] == glc
        assert sorted(sim.get_metabolite_reactions(m)) == sorted(sim.get_metabolite_reactions(glc))
        assert sim.get_metabolite_producers(m) == sim.get_metabolite_producers(glc)
        assert sim.get_metabolite_consumers(m) == sim.get_metabolite_consumers(glc)
        assert sim.metabolite_elements(m) == sim.metabolite_elements(glc)
    result = sim.simulate()
    for m in both_forms("glc__D_e", "metabolite"):
        assert result.get_metabolite(m, format="dict") == result.get_metabolite(glc, format="dict")


def test_constructor_inputs_both_forms(backend):
    for ex in both_forms("EX_glc__D_e"):
        sim = load(backend, envcond={ex: (-1, 1000)}, constraints={both_forms("PGI")[1]: (-1000, 1000)})
        assert list(sim.environmental_conditions) == [native(backend, "EX_glc__D_e")]
        assert sim.get_reaction_bounds("EX_glc__D_e") == (-1, 1000)
        assert sim.simulate().objective_value == pytest.approx(0.0484, abs=1e-3)
    reference = {sim.to_canonical(k, "reaction"): v for k, v in sim.simulate(method="pFBA").fluxes.items()}
    sim2 = load(backend, reference=reference)
    assert set(sim2.reference) == set(sim2.reactions)


def test_environment_apply_other_form(backend):
    # the non-native form: unprefixed ids on REFRAMED, prefixed ids on COBRApy
    from mewpy.simulation import Environment

    sim = load(backend)
    exchanges = ["EX_glc__D_e", "EX_o2_e", "EX_nh4_e", "EX_pi_e", "EX_h2o_e", "EX_h_e"]
    env = Environment.from_reactions(exchanges if backend == "reframed" else ["R_" + r for r in exchanges])
    env.apply(sim, exclusive=False)
    assert sim.get_reaction_bounds(native(backend, "EX_glc__D_e"))[0] == -10.0


def test_unknown_id_unchanged_and_old_error(sim, backend):
    assert sim.resolve("NOPE", "reaction") == "NOPE"
    assert sim.resolve("R_NOPE", "reaction") == "R_NOPE"
    assert sim.resolve_keys({"NOPE": 0}) == {"NOPE": 0}
    expected_error = KeyError if backend == "cobra" else ValueError
    with pytest.raises(expected_error):
        sim.get_reaction("NOPE")
    with pytest.raises(KeyError):
        sim.simulate(constraints={"NOPE": 0})


def test_ambiguous_id_raises():
    from reframed.core.cbmodel import CBModel, CBReaction
    from reframed.core.model import Compartment, Metabolite

    from mewpy.simulation import get_simulator

    model = CBModel("tiny")
    model.add_compartment(Compartment("c"))
    model.add_metabolite(Metabolite("M_a_c", compartment="c"))
    model.add_reaction(CBReaction("X", stoichiometry={"M_a_c": 1}, lb=0, ub=1))
    model.add_reaction(CBReaction("R_R_X", stoichiometry={"M_a_c": -1}, lb=0, ub=1))
    sim = get_simulator(model)
    assert sim.resolve("X") == "X"
    assert sim.resolve("R_R_X") == "R_R_X"
    with pytest.raises(ValueError, match="Ambiguous"):
        sim.resolve("R_X")


def test_cache_refreshes_after_add_and_remove(sim, backend):
    new = sim.to_native("NEWRX", "reaction")
    assert sim.resolve("NEWRX") == "NEWRX"  # builds the cache before the change
    sim.add_reaction(new, stoichiometry={"glc__D_e": -1, "M_g6p_c": 1}, lb=0, ub=1)
    assert sim.resolve("NEWRX") == new
    assert sim.resolve("R_NEWRX") == new
    assert set(sim.get_reaction_metabolites("NEWRX")) == {
        native(backend, "glc__D_e", "metabolite"),
        native(backend, "g6p_c", "metabolite"),
    }
    sim.remove_reactions(["NEWRX"])
    assert new not in sim.reactions
    assert sim.resolve("NEWRX") == "NEWRX"


def test_cache_detects_direct_model_change(sim, backend):
    new = sim.to_native("DIRECT", "reaction")
    assert sim.resolve("DIRECT") == "DIRECT"
    if backend == "cobra":
        from cobra import Reaction

        sim.model.add_reactions([Reaction(new)])
    else:
        from reframed.core.cbmodel import CBReaction

        sim.model.add_reaction(CBReaction(new, stoichiometry={}))
    assert sim.resolve("R_DIRECT") == new


def test_germ_accepts_prefixed_ids():
    from mewpy.io import read_sbml
    from mewpy.simulation import get_simulator

    sim = get_simulator(read_sbml(EC_CORE_XML, regulatory=False, warnings=False))
    assert sim.resolve("R_PGI") == "PGI"
    res = sim.simulate(constraints={"R_PGI": 0})
    assert res.objective_value == pytest.approx(PGI_KO_GROWTH, abs=1e-4)
    assert sim.get_reaction_bounds("R_PGI") == sim.get_reaction_bounds("PGI")
    assert set(sim.FVA(reactions=["R_PGI"], obj_frac=0.9)) == {"PGI"}
    assert sim.essential_genes() and all(not g.startswith("G_") for g in sim.essential_genes())
    assert sim.get_gpr("R_PGI") == sim.get_gpr("PGI")


MEDIUM = ["glc__D", "o2", "nh4", "pi", "h2o", "h", "co2"]


def test_environment_compounds_use_simulator_prefix():
    from mewpy.simulation import Environment

    growth, compounds = {}, {}
    for backend in BACKENDS:
        sim = load(backend)
        env = Environment.from_compounds(MEDIUM, model=sim)
        assert list(env) == [native(backend, f"EX_{c}_e") for c in MEDIUM]
        compounds[backend] = env.get_compounds(model=sim)
        assert Environment.from_compounds(MEDIUM, model=sim.model).keys() == env.keys()  # a model works too
        env.apply(sim)
        growth[backend] = sim.simulate().objective_value
    assert compounds["cobra"] == compounds["reframed"] == MEDIUM
    assert growth["cobra"] == pytest.approx(growth["reframed"], abs=1e-6)
    assert growth["cobra"] > 0.1  # the medium (uptake capped at 10 for each compound) supports growth


def test_environment_compounds_backward_compatible():
    # passes before the change too: it guards the defaults kept for callers that pass no model
    from mewpy.simulation import Environment

    assert list(Environment.from_compounds(["glc__D"])) == ["EX_glc__D_e"]
    assert list(Environment.from_compounds(["glc__D"], prefix="R_")) == ["R_EX_glc__D_e"]
    assert list(Environment.from_compounds(["glc__D"], fmt_func="EX_{}_b")) == ["EX_glc__D_b"]
    from_model = {backend: Environment.from_model(load(backend)).get_compounds() for backend in BACKENDS}
    assert from_model["cobra"] == from_model["reframed"]


# ----------------------------------------------------------------------------------------------------------------------
# Round 3: one behaviour on every backend (COBRApy, REFRAMED and, where it applies, GERM)
# ----------------------------------------------------------------------------------------------------------------------
def load_germ():
    from mewpy.io import read_sbml
    from mewpy.simulation import get_simulator

    return get_simulator(read_sbml(EC_CORE_XML, regulatory=False, warnings=False))


def all_backends():
    return {"cobra": load("cobra"), "reframed": load("reframed"), "germ": load_germ()}


def test_pfba_keeps_the_shared_fraction_on_every_backend():
    from mewpy.util.constants import ModelConstants

    assert ModelConstants.PFBA_OBJ_FRAC == 0.999
    results = {name: sim.simulate(method="pFBA") for name, sim in all_backends().items()}
    growth = {
        name: res.fluxes[native("reframed" if name == "reframed" else "cobra", "Biomass_Ecoli_core")]
        for name, res in results.items()
    }
    for name in results:
        assert growth[name] == pytest.approx(ModelConstants.PFBA_OBJ_FRAC * 0.8739215, abs=1e-6), name
        assert results[name].objective_value == pytest.approx(results["cobra"].objective_value, abs=1e-5), name


INFEASIBLE = {"ATPM": (500, 1000)}


def test_infeasible_results_follow_one_convention():
    for name, sim in all_backends().items():
        constraints = {native(name if name == "reframed" else "cobra", k): v for k, v in INFEASIBLE.items()}
        assert sim.simulate(constraints=constraints, slim=True) is None, name
        res = sim.simulate(constraints=constraints)
        assert (res.objective_value, res.fluxes) == (None, None), name
        assert isinstance(sim.simulate(slim=True), float), name
        if name != "germ":
            assert isinstance(sim.simulate(method="pFBA", slim=True), float), name
            assert sim.simulate(method="pFBA", constraints=constraints, slim=True) is None, name


def test_fva_layout_is_the_same_on_every_backend():
    frames, dicts = {}, {}
    for name, sim in all_backends().items():
        ids = [native(name if name == "reframed" else "cobra", r) for r in ("PGI", "ACKr")]
        dicts[name] = {sim.to_canonical(k, "reaction"): v for k, v in sim.FVA(reactions=ids).items()}
        frames[name] = sim.FVA(reactions=ids, format="df")
    for name, frame in frames.items():
        assert list(frame.columns) == ["Reaction ID", "Minimum", "Maximum"], name
        assert list(frame.index) == [0, 1], name
        assert not frame.isna().any().any(), name
        assert sorted(dicts[name]) == ["ACKr", "PGI"], name
        for r_id, (low, high) in dicts[name].items():
            assert [low, high] == pytest.approx(dicts["cobra"][r_id], abs=1e-4), (name, r_id)


def test_set_objective_replaces_the_objective(backend):
    sim = load(backend)
    sim.set_objective("EX_o2_e")
    assert {k: v for k, v in sim.objective.items() if v} == {native(backend, "EX_o2_e"): 1.0}
    assert sim.simulate(maximize=False).objective_value == pytest.approx(-60.0, abs=1e-6)


def test_result_reports_the_sense_that_was_optimised(backend):
    sim = load(backend)
    assert sim.simulate().maximize is True
    for method in ("pFBA", "lMOMA"):
        res = sim.simulate(method=method, constraints={"PGI": 0})
        assert res.maximize is False, method
        assert "Objective (minimize)" in repr(res), method
    if backend == "reframed":  # ROOM on COBRApy/GLPK takes minutes
        wt = sim.simulate(method="pFBA").fluxes
        reference = {r: wt[r] for r in ("R_PGI", "R_PFK", "R_CS")}
        assert sim.simulate(method="ROOM", constraints={"PGI": 0}, reference=reference).maximize is False


def test_transport_is_defined_once():
    counts, genes = {}, {}
    for name, sim in all_backends().items():
        transport = sim.get_transport_reactions()
        assert not set(transport) & set(sim.get_exchange_reactions()), name
        counts[name] = sorted(sim.to_canonical(r, "reaction") for r in transport)
        genes[name] = sorted(sim.to_canonical(g, "gene") for g in sim.get_transport_genes())
    assert counts["cobra"] == counts["reframed"] == counts["germ"]
    assert genes["cobra"] == genes["reframed"] == genes["germ"]
    assert "GLCpts" in counts["cobra"] and "NADTRHD" not in counts["cobra"]
    assert "ACt2r" in counts["cobra"]  # acetate-proton symport moves acetate
    assert not {"ATPS4r", "NADH16", "CYTBD", "THD2"} & set(counts["cobra"])  # proton pumps are enzymes
    assert not {"b1602", "b1603"} & set(genes["cobra"])  # pntAB also catalyse the internal NADTRHD
