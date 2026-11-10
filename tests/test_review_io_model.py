"""Regression tests for io/, model/, cobra/ and util/parsing defects found in the 2026-11-09 review.

Each test reproduces one finding (F1, F3-F12, P-a..P-c, U-*); see the docstrings.
"""

import subprocess
import sys
import warnings
from pathlib import Path

import pytest

MODELS_PATH = Path(__file__).parent.joinpath("data")
EC_CORE = str(MODELS_PATH.joinpath("e_coli_core.xml"))
EC_CORE_TRN = str(MODELS_PATH.joinpath("e_coli_core_trn.csv"))
CHASSAGNOLE = str(MODELS_PATH.joinpath("chassagnole2002.xml"))
PARSING_PY = Path(__file__).parents[1].joinpath("src", "mewpy", "util", "parsing.py")


def _read_ecoli():
    from mewpy.io import read_sbml

    return read_sbml(EC_CORE, regulatory=False, warnings=False)


def _read_ecoli_trn():
    from mewpy.io import Engines, Reader, read_model

    return read_model(
        Reader(Engines.MetabolicSBML, EC_CORE),
        Reader(Engines.BooleanRegulatoryCSV, EC_CORE_TRN, sep=",", id_col=0, rule_col=2, aliases_cols=[1], header=0),
        warnings=False,
    )


def _kinetic_sbml(path, laws, amount=False, function=None, rules=None, extra_params=None, no_law=False):
    """Writes a two-species kinetic SBML (A -> B per law) and returns its path.

    A=2, B=1 (concentrations, or amounts in a compartment of size 2 if amount=True), k=0.5.
    """
    import libsbml

    doc = libsbml.SBMLDocument(3, 1)
    m = doc.createModel()
    m.setId("t")
    c = m.createCompartment()
    c.setId("c")
    c.setSize(2.0 if amount else 1.0)
    c.setConstant(True)
    for sid, val in (("A", 2.0), ("B", 1.0)):
        s = m.createSpecies()
        s.setId(sid)
        s.setCompartment("c")
        if amount:
            s.setInitialAmount(val)
        else:
            s.setInitialConcentration(val)
        s.setConstant(False)
        s.setBoundaryCondition(False)
        s.setHasOnlySubstanceUnits(False)
    params = {"k": (0.5, True)}
    params.update(extra_params or {})
    for pid, (val, const) in params.items():
        p = m.createParameter()
        p.setId(pid)
        p.setValue(val)
        p.setConstant(const)
    if function:
        fid, fname, lam = function
        fd = m.createFunctionDefinition()
        fd.setId(fid)
        fd.setName(fname)
        fd.setMath(libsbml.parseL3Formula(lam))
    for var, formula in (rules or {}).items():
        ar = m.createAssignmentRule()
        ar.setVariable(var)
        ar.setMath(libsbml.parseL3Formula(formula))
    for i, law in enumerate(laws):
        r = m.createReaction()
        r.setId(f"R{i}")
        r.setReversible(False)
        sr = r.createReactant()
        sr.setSpecies("A")
        sr.setStoichiometry(1)
        sr.setConstant(True)
        sp = r.createProduct()
        sp.setSpecies("B")
        sp.setStoichiometry(1)
        sp.setConstant(True)
        kl = r.createKineticLaw()
        kl.setMath(libsbml.parseL3Formula(law))
    if no_law:
        r = m.createReaction()
        r.setId("NOLAW")
        r.setReversible(False)
        sr = r.createReactant()
        sr.setSpecies("A")
        sr.setStoichiometry(1)
        sr.setConstant(True)
    libsbml.writeSBMLToFile(doc, str(path))
    return str(path)


def _ode_and_deriv(model, y):
    import numpy as np

    return np.array(model.get_ode()(0, y), dtype=float), np.array(model.deriv(0, y), dtype=float)


# -----------------------------------------------------------------------------
# F1: SBML writer gave the upper-bound parameter the lower-bound id
# -----------------------------------------------------------------------------
def test_f1_sbml_roundtrip_preserves_bounds(tmp_path):
    from mewpy.io import Engines, Writer, read_sbml, write_model

    model = _read_ecoli()
    model.reactions["ATPM"].bounds = (-5.0, 7.0)
    out = str(tmp_path / "rt.xml")
    write_model(Writer(Engines.MetabolicSBML, io=out, model=model), warnings=False)
    model2 = read_sbml(out, regulatory=False, warnings=False)
    changed = {
        k: (model.reactions[k].bounds, model2.reactions[k].bounds)
        for k in model.reactions
        if model.reactions[k].bounds != model2.reactions[k].bounds
    }
    assert changed == {}


# -----------------------------------------------------------------------------
# F3: pFBA crashed on non-enzyme models and ignored a user objective in the pre-solve
# -----------------------------------------------------------------------------
def test_f3_pfba_plain_model_default_reactions():
    from reframed import load_cbmodel

    from mewpy.cobra import pFBA
    from mewpy.simulation import get_simulator
    from mewpy.solvers.solution import Status

    sim = get_simulator(load_cbmodel(EC_CORE))
    solution = pFBA(sim)
    assert solution.status == Status.OPTIMAL
    assert not any(k.endswith("_p") or k.endswith("_n") for k in solution.values)


def test_f3_pfba_custom_objective():
    from reframed import load_cbmodel

    from mewpy.cobra import pFBA
    from mewpy.simulation import get_simulator
    from mewpy.solvers.solution import Status

    sim = get_simulator(load_cbmodel(EC_CORE))
    solution = pFBA(sim, objective={"R_ATPM": 1}, reactions=list(sim.reactions))
    assert solution.status == Status.OPTIMAL
    assert solution.values["R_ATPM"] > 100


# -----------------------------------------------------------------------------
# F4: SBML-qual reader turned constant (TRUE) function terms into empty expressions
# -----------------------------------------------------------------------------
def test_f4_regulatory_sbml_roundtrip_keeps_constant_rules(tmp_path):
    from mewpy.io import Engines, Reader, Writer, read_model, write_model

    model = _read_ecoli_trn()
    out = str(tmp_path / "reg.xml")
    write_model(Writer(Engines.RegulatorySBML, io=out, model=model), warnings=False)
    model2 = read_model(Reader(Engines.RegulatorySBML, out), warnings=False)

    lost = [
        k
        for k, i in model.interactions.items()
        if any(not e.is_none for e in i.regulatory_events.values())
        and all(e.is_none for e in model2.interactions[k].regulatory_events.values())
    ]
    assert lost == []
    table = model2.interactions["b0008_interaction"].regulatory_truth_table
    assert sorted(table["result"].tolist()) == [0, 1]


def test_f4_qual_species_level_range_check():
    from mewpy.io.dto import DataTransferObject, VariableRecord
    from mewpy.io.engines.regulatory_sbml import RegulatorySBML

    engine = RegulatorySBML(io="unused.xml", config={})
    engine._dto = DataTransferObject()
    engine.dto.variables["g"] = VariableRecord(id="g", coefficients={0.0, 1.0})
    engine._update_qual_species_coefficients("g", 0.5)
    assert engine.warnings == []
    engine._update_qual_species_coefficients("g", 5)
    assert len(engine.warnings) == 1
    assert engine.dto.variables["g"].coefficients == {0.0, 0.5, 1.0, 5.0}


# -----------------------------------------------------------------------------
# F5: cobra/reframed engines stored objective coefficients and charges as strings
# -----------------------------------------------------------------------------
def test_f5_cobra_engine_objective_is_numeric_and_signed(tmp_path):
    import cobra

    from mewpy.io import Engines, Writer, read_cbmodel, write_model

    cobra_model = cobra.io.read_sbml_model(EC_CORE)
    germ = read_cbmodel(cobra_model, warnings=False)
    assert {r.id: c for r, c in germ.objective.items()} == {"Biomass_Ecoli_core": 1.0}
    write_model(Writer(Engines.MetabolicSBML, io=str(tmp_path / "c.xml"), model=germ), warnings=False)

    cobra_model.objective = {cobra_model.reactions.ATPM: 2.5}
    cobra_model.objective_direction = "min"
    germ2 = read_cbmodel(cobra_model, warnings=False)
    assert {r.id: c for r, c in germ2.objective.items()} == {"ATPM": -2.5}


def test_f5_cobra_engine_model_simulates():
    import cobra

    from mewpy.io import read_cbmodel
    from mewpy.simulation import get_simulator

    germ = read_cbmodel(cobra.io.read_sbml_model(EC_CORE), warnings=False)
    assert get_simulator(germ).simulate().objective_value == pytest.approx(0.8739, abs=1e-3)


def test_f5_reframed_engine_charge_and_compartments(tmp_path):
    from reframed import load_cbmodel

    from mewpy.io import Engines, Writer, read_cbmodel, write_model

    germ = read_cbmodel(load_cbmodel(EC_CORE), cobrapy=False, reframed=True, warnings=False)
    assert germ.compartments == {"c": "cytosol", "e": "extracellular"}
    assert all(m.charge is None or isinstance(m.charge, int) for m in germ.metabolites.values())
    write_model(Writer(Engines.MetabolicSBML, io=str(tmp_path / "r.xml"), model=germ), warnings=False)


# -----------------------------------------------------------------------------
# F6 / P-b: kinetic laws with e-notation, unary minus and function definitions
# -----------------------------------------------------------------------------
def test_f6_parsing_scientific_notation_and_unary_minus():
    import numexpr as ne

    from mewpy.util.parsing import Arithmetic, build_tree

    cases = {
        "2.5e-05 * A": 2.5e-05 * 2.0,
        "A * 1e-3 + 1E+2": 2.0e-3 + 100.0,
        "exp(-A/k)": 0.0183156,
        "2 * -A": -4.0,
        "-A^2": -4.0,
        "(-A) * k": -1.0,
        "k1e - 3": -2.0,
    }
    for expr, expected in cases.items():
        tree = build_tree(expr, Arithmetic)
        infix = tree.replace({"A": "2.0", "k": "0.5", "k1e": "1.0"}).to_infix()
        assert float(ne.evaluate(infix, local_dict={})) == pytest.approx(expected, rel=1e-5), (expr, infix)
    assert build_tree("2.5e-05 * A", Arithmetic).get_parameters() == {"A"}


@pytest.mark.parametrize(
    "law,expected_rate",
    [("2.5e-05 * A", 5e-05), ("k * exp(-A/k)", 0.5 * 0.01831563888873418)],
)
def test_f6_kinetic_rate_laws(tmp_path, law, expected_rate):
    import numpy as np

    from mewpy.io.sbml import load_ODEModel

    model = load_ODEModel(_kinetic_sbml(tmp_path / "k.xml", [law]))
    ode, der = _ode_and_deriv(model, [2.0, 1.0])
    np.testing.assert_allclose(ode, [-expected_rate, expected_rate], rtol=1e-6)
    np.testing.assert_allclose(der, ode, rtol=1e-6)


def test_f6_kinetic_function_definition_by_id_with_arguments(tmp_path):
    import numpy as np

    from mewpy.io.sbml import load_ODEModel

    function = ("MM", "Michaelis Menten", "lambda(S, Km, 2*S/(Km+S))")
    path = _kinetic_sbml(tmp_path / "f.xml", ["MM(A, k)"], function=function)
    model = load_ODEModel(path)
    assert list(model.function_definition) == ["MM"]
    rate = 2 * 2.0 / (0.5 + 2.0)
    ode, der = _ode_and_deriv(model, [2.0, 1.0])
    np.testing.assert_allclose(ode, [-rate, rate], rtol=1e-6)
    np.testing.assert_allclose(der, ode, rtol=1e-6)


# -----------------------------------------------------------------------------
# F8: assignment rules ignored, initialAmount read as NaN, reaction without kinetic law crashed
# -----------------------------------------------------------------------------
def test_f8_assignment_rules_are_applied(tmp_path):
    import numpy as np

    from mewpy.io.sbml import load_ODEModel

    path = _kinetic_sbml(
        tmp_path / "ar.xml",
        ["kr * A"],
        rules={"kr": "2 * kq", "kq": "A"},
        extra_params={"kr": (0.0, False), "kq": (0.0, False)},
    )
    model = load_ODEModel(path)
    ode, der = _ode_and_deriv(model, [3.0, 1.0])
    np.testing.assert_allclose(ode, [-18.0, 18.0], rtol=1e-6)
    np.testing.assert_allclose(der, ode, rtol=1e-6)


def test_f8_initial_amount_converted_to_concentration(tmp_path):
    from mewpy.io.sbml import load_ODEModel

    model = load_ODEModel(_kinetic_sbml(tmp_path / "am.xml", ["k * A"], amount=True))
    assert dict(model.concentrations) == {"A": 1.0, "B": 0.5}


def test_f8_reaction_without_kinetic_law_is_skipped(tmp_path):
    from mewpy.io.sbml import load_ODEModel

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = load_ODEModel(_kinetic_sbml(tmp_path / "nl.xml", ["k * A"], no_law=True))
    assert list(model.ratelaws) == ["R0"]
    assert any("NOLAW" in str(w.message) for w in caught)


# -----------------------------------------------------------------------------
# F9: GeckoModel crashes
# -----------------------------------------------------------------------------
def test_f9_gecko_accepts_protein_properties_dataframe():
    from reframed import load_cbmodel

    from mewpy.model.gecko import GeckoModel, ModelList

    properties = ModelList().protein_properties()
    model = GeckoModel(load_cbmodel(EC_CORE), protein_properties=properties)
    assert model.protein_properties is properties


@pytest.fixture(scope="module")
def multi_pool_limited():
    import pandas as pd

    from mewpy.model.gecko import GeckoModel

    in_model = {"P00549": 0.1, "P31373": 0.1, "P31382": 0.1, "P39708": 0.1, "P39714": 0.1, "P39726": 0.1, "Q01574": 0.1}
    not_in_model = {"P10591": 0.1, "P31383": 0.1, "P32471": 0.1}
    model = GeckoModel("multi-pool")
    model.limit_proteins(fractions=pd.concat([pd.Series(in_model), pd.Series(not_in_model)]))
    return model, in_model


def test_f9_gecko_multi_pool_limit_proteins(multi_pool_limited):
    from mewpy.simulation.reframed import GeckoSimulation

    model, in_model = multi_pool_limited
    assert len(model.proteins) - len(model.pool_proteins) - len(in_model) == 0
    assert all(model.reactions[rxn].ub > 0 for rxn in model.individual_protein_exchanges)
    assert GeckoSimulation(model).simulate().objective_value > 0.05


def test_f9_gecko_adjust_pool_bounds():
    """A protein pool that limits growth by a margin far above the solver tolerances.

    In the bundled multi-pool model the measured pools change growth by ~1e-3 and the bound
    corrections are ~1e-7, below SCIP's feasibility tolerance (1e-6): whether they are reported
    depends on the column order (PYTHONHASHSEED). Here the pool supplies 1 unit of protein per
    unit of biomass, so reaching a growth of 0.5 needs a correction of 0.4.
    """
    from reframed import load_cbmodel
    from reframed.core.cbmodel import CBReaction
    from reframed.core.model import Metabolite

    from mewpy.model.gecko import GeckoModel
    from mewpy.simulation.reframed import GeckoSimulation

    cbmodel = load_cbmodel(EC_CORE)
    cbmodel.add_metabolite(Metabolite("prot_P1_c", "prot_P1", "c"))
    cbmodel.add_reaction(
        CBReaction("prot_P1_exchange", reversible=False, stoichiometry={"prot_P1_c": 1.0}, lb=0, ub=0.1)
    )
    cbmodel.reactions["R_Biomass_Ecoli_core"].stoichiometry["prot_P1_c"] = -1.0
    model = GeckoModel(cbmodel, biomass_reaction_id="R_Biomass_Ecoli_core")
    assert model.individual_protein_exchanges == {"prot_P1_exchange"}
    assert GeckoSimulation(model).simulate().objective_value == pytest.approx(0.1, rel=1e-6)

    result = model.adjust_pool_bounds(min_objective=0.5, inplace=True)
    assert list(result.columns) == ["reaction", "addition", "original", "adjusted"]
    assert result["reaction"].tolist() == ["prot_P1_exchange"]
    assert result["addition"].iloc[0] == pytest.approx(0.4, rel=1e-6)
    assert model.reactions["prot_P1_exchange"].ub == pytest.approx(0.5, rel=1e-6)
    assert GeckoSimulation(model).simulate().objective_value == pytest.approx(0.5, rel=1e-6)


# -----------------------------------------------------------------------------
# F10: cobra/util transformations
# -----------------------------------------------------------------------------
def test_f10_convert_gpr_to_dnf_updates_model():
    from reframed import load_cbmodel
    from reframed.io.sbml import parse_gpr_rule

    from mewpy.cobra.util import convert_gpr_to_dnf
    from mewpy.simulation import get_simulator

    model = load_cbmodel(EC_CORE)
    model.reactions["R_PFK"].gpr = parse_gpr_rule("(G_b3916 or G_b1723) and G_b0001")
    sim = get_simulator(model)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        convert_gpr_to_dnf(sim)
    gpr = sim.get_reaction("R_PFK").gpr
    terms = sorted(gpr.replace("(", "").replace(")", "").split(" or "))
    assert terms == ["G_b0001 and G_b1723", "G_b0001 and G_b3916"] or terms == [
        "G_b1723 and G_b0001",
        "G_b3916 and G_b0001",
    ]


def test_f10_convert_to_irreversible_forced_reverse_flux():
    from reframed import load_cbmodel

    from mewpy.cobra.util import convert_to_irreversible
    from mewpy.simulation import SStatus, get_simulator

    model = load_cbmodel(EC_CORE)
    model.reactions["R_EX_glc__D_e"].lb = -10
    model.reactions["R_EX_glc__D_e"].ub = -10
    irr, mapping = convert_to_irreversible(model)
    assert isinstance(mapping, dict) and mapping["R_EX_glc__D_e"] == "R_EX_glc__D_e_REV"
    assert irr.get_reaction_bounds("R_EX_glc__D_e") == (0, 0)
    assert irr.get_reaction_bounds("R_EX_glc__D_e_REV") == (10, 10)
    original = get_simulator(model).simulate()
    result = irr.simulate()
    assert result.status == SStatus.OPTIMAL
    assert result.objective_value == pytest.approx(original.objective_value, rel=1e-6)


def test_f10_enzyme_constraints_use_original_reaction_kcats():
    from reframed import load_cbmodel

    from mewpy.cobra.util import add_enzyme_constraints
    from mewpy.simulation import get_simulator

    model = load_cbmodel(EC_CORE)
    sim = get_simulator(model)
    prot_mw = {g: {"protein": g[2:], "mw": 1} for g in sim.genes}
    enz_kcats = {"G_b1723": {"R_PFK": {"protein": "b1723", "kcat": 100.0}}}
    ec = add_enzyme_constraints(model, prot_mw=prot_mw, enz_kcats=enz_kcats)
    coeffs = {r: ec.get_reaction(r).stoichiometry.get("prot_b1723_c") for r in ec.reactions if r.startswith("R_PFK")}
    assert pytest.approx(-0.01) in [c for c in coeffs.values() if c is not None], coeffs


def _ecoli_cobra_and_reframed():
    """e_coli_core loaded by COBRApy and by REFRAMED. REFRAMED reads the +/-1000 default bounds as +/-inf:
    they are set back to +/-1000 so that both backends start from the same model."""
    import cobra
    from reframed import load_cbmodel

    cb_model = load_cbmodel(EC_CORE)
    for rxn in cb_model.reactions.values():
        rxn.lb = max(rxn.lb, -1000.0)
        rxn.ub = min(rxn.ub, 1000.0)
    return cobra.io.read_sbml_model(EC_CORE), cb_model


def _canonical_model(sim):
    """Reactions (bounds, stoichiometry, GPR) and metabolites of a simulator, backend prefixes stripped."""
    import re

    from mewpy.cobra.util import gpr_to_dnf

    def gpr(rule):
        # COBRApy gives '' and REFRAMED None for no GPR, and the backends order the genes differently:
        # compare the GPRs as sets of isozymes (sets of genes)
        if not rule:
            return frozenset()
        rule = re.sub(r"[\w.]+", lambda m: sim.to_canonical(m.group(0), "gene"), rule)
        terms = gpr_to_dnf(rule).split(" or ")
        return frozenset(frozenset(term.strip("()").split(" and ")) for term in terms)

    reactions = {}
    for r_id in sim.reactions:
        rxn = sim.get_reaction(r_id)
        stoichiometry = {sim.to_canonical(m, "metabolite"): round(v, 9) for m, v in rxn.stoichiometry.items()}
        # split_isozymes numbers the isozymes (_No1, _No2, ...) in the GPR gene order of each backend:
        # isozyme reactions are matched by their base id and GPR
        base_id = re.sub(r"_No\d+$", "", sim.to_canonical(r_id, "reaction"))
        reactions[(base_id, gpr(rxn.gpr))] = (rxn.lb, rxn.ub, stoichiometry)
    return reactions, sorted(sim.to_canonical(m, "metabolite") for m in sim.metabolites)


def test_prefix_parity_convert_to_irreversible():
    from mewpy.cobra.util import convert_to_irreversible

    results = []
    for model in _ecoli_cobra_and_reframed():
        sim, mapping = convert_to_irreversible(model)
        canonical_map = {sim.to_canonical(k, "reaction"): sim.to_canonical(v, "reaction") for k, v in mapping.items()}
        results.append((_canonical_model(sim), canonical_map, sim.simulate().objective_value))
    (cobra_model, cobra_map, cobra_growth), (cb_model, cb_map, cb_growth) = results
    assert cobra_model == cb_model
    assert cobra_map == cb_map
    assert cobra_growth == pytest.approx(cb_growth, rel=1e-9)


def test_prefix_parity_add_enzyme_constraints():
    from mewpy.cobra.util import add_enzyme_constraints

    results = []
    for model in _ecoli_cobra_and_reframed():
        sim = add_enzyme_constraints(model)
        assert sim.resolve("prot_pool_exchange", "reaction") in sim.reactions
        assert sim.resolve("prot_pool_c", "metabolite") in sim.metabolites
        results.append((_canonical_model(sim), sim.simulate().objective_value))
    (cobra_model, cobra_growth), (cb_model, cb_growth) = results
    assert cobra_model == cb_model
    assert cobra_growth is not None and cobra_growth > 0
    assert cobra_growth == pytest.approx(cb_growth, rel=1e-9)


def test_enzyme_constraints_model_is_feasible():
    """The draw_prot reactions (whose GPR is their gene) were given the enzyme usage of a catalysed reaction,
    which replaced their protein production by a consumption: the model could not grow."""
    from reframed import load_cbmodel

    from mewpy.cobra.util import add_enzyme_constraints
    from mewpy.simulation import SStatus, get_simulator

    model = load_cbmodel(EC_CORE)
    growth = get_simulator(model).simulate().objective_value
    sim = add_enzyme_constraints(model)
    assert sim.get_reaction("draw_prot_b0351").stoichiometry["prot_b0351_c"] == 1
    result = sim.simulate()
    assert result.status == SStatus.OPTIMAL
    # MW = kcat = 1 and an unbounded protein pool: the enzyme constraints do not limit growth
    assert result.objective_value == pytest.approx(growth, rel=1e-6)


def test_enzyme_constraints_unprefixed_gene_in_prefixed_backend():
    """A REFRAMED gene id without the G_ prefix kept its first two characters chopped off (b9999 -> 9999)."""
    from reframed import load_cbmodel
    from reframed.core.cbmodel import Gene
    from reframed.io.sbml import parse_gpr_rule

    from mewpy.cobra.util import add_enzyme_constraints

    model = load_cbmodel(EC_CORE)
    model.add_gene(Gene("b9999"))
    model.reactions["R_PFK"].gpr = parse_gpr_rule("G_b3916 or b9999")
    sim = add_enzyme_constraints(model)
    assert "prot_b9999_c" in sim.metabolites
    assert "draw_prot_b9999" in sim.reactions


def test_smoment_enzyme_prefix_with_or_without_sbml_prefix():
    from reframed import load_cbmodel
    from reframed.core.cbmodel import CBReaction

    from mewpy.model.smoment import SMomentModel

    model = load_cbmodel(EC_CORE)
    for protein in ("b0001", "b0002"):
        model.add_reaction(CBReaction(f"R_ENZYME_DELIVERY_{protein}", stoichiometry={}, lb=0, ub=1))
    assert sorted(SMomentModel(model).proteins) == ["b0001", "b0002"]
    assert sorted(SMomentModel(model, enzyme_reaction_prefix="ENZYME_DELIVERY_").proteins) == ["b0001", "b0002"]
    assert SMomentModel(model, enzyme_reaction_prefix="R__TG_ER_").proteins == []


# -----------------------------------------------------------------------------
# F11: repeated species references overwrote each other's stoichiometry
# -----------------------------------------------------------------------------
def test_f11_sbml_repeated_species_stoichiometry(tmp_path):
    import libsbml

    from mewpy.io import read_sbml

    doc = libsbml.SBMLDocument(libsbml.SBMLNamespaces(3, 1, "fbc", 2))
    doc.setPackageRequired("fbc", False)
    m = doc.createModel()
    m.setId("dup")
    m.getPlugin("fbc").setStrict(True)
    c = m.createCompartment()
    c.setId("c")
    c.setConstant(True)
    for s in ("A", "B", "C"):
        sp = m.createSpecies()
        sp.setId("M_" + s)
        sp.setCompartment("c")
        sp.setConstant(False)
        sp.setBoundaryCondition(False)
        sp.setHasOnlySubstanceUnits(False)
    for pid, val in (("lb", -1000.0), ("ub", 1000.0)):
        p = m.createParameter()
        p.setId(pid)
        p.setValue(val)
        p.setConstant(True)
    r = m.createReaction()
    r.setId("R_1")
    r.setReversible(True)
    r.setFast(False)
    for side, species in ((r.createReactant, ("A", "A", "B")), (r.createProduct, ("B", "C"))):
        for s in species:
            ref = side()
            ref.setSpecies("M_" + s)
            ref.setStoichiometry(1)
            ref.setConstant(True)
    fbc = r.getPlugin("fbc")
    fbc.setLowerFluxBound("lb")
    fbc.setUpperFluxBound("ub")
    out = str(tmp_path / "dup.xml")
    libsbml.writeSBMLToFile(doc, out)
    model = read_sbml(out, regulatory=False, warnings=False)
    assert {k.id: v for k, v in model.reactions["1"].stoichiometry.items()} == {"A": -2.0, "C": 1.0}


# -----------------------------------------------------------------------------
# F12: Reader dropped regulator_col; read_csv/read_cbmodel default flags forced a second reader
# -----------------------------------------------------------------------------
def test_f12_reader_regulator_col_and_flags(tmp_path):
    from reframed import load_cbmodel

    from mewpy.io import Engines, Reader, read_cbmodel, read_csv, read_model

    csv = tmp_path / "tr.csv"
    csv.write_text("R1,g1\nR2,g1\nR1,g2\n")
    model = read_model(
        Reader(Engines.TargetRegulatorRegulatoryCSV, str(csv), target_col=1, regulator_col=0), warnings=False
    )
    assert sorted(model.targets) == ["g1", "g2"]
    assert sorted(model.regulators) == ["R1", "R2"]

    csv2 = tmp_path / "tr2.csv"
    csv2.write_text("g1,R1\ng1,R2\ng2,R1\n")
    assert sorted(read_csv(str(csv2), target_regulator=True, warnings=False).targets) == ["g1", "g2"]
    assert len(read_cbmodel(load_cbmodel(EC_CORE), reframed=True, warnings=False).reactions) == 95


# -----------------------------------------------------------------------------
# P-a: Boolean AND/OR precedence; P-c: recursion limit set on import
# -----------------------------------------------------------------------------
def test_pa_boolean_and_binds_tighter_than_or():
    from mewpy.util.parsing import evaluate_expression

    assert evaluate_expression("a or b and c", ["a"]) is True
    assert evaluate_expression("b and c or a", ["a"]) is True
    assert evaluate_expression("not a and b", ["b"]) is True
    assert evaluate_expression("a and not b or c", ["a", "b"]) is False


def test_pc_recursion_limit_not_raised_on_import():
    code = (
        "import importlib.util, sys\n"
        f"spec = importlib.util.spec_from_file_location('parsing', {str(PARSING_PY)!r})\n"
        "mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)\n"
        "assert sys.getrecursionlimit() < 100000, sys.getrecursionlimit()\n"
        "gpr = ' or '.join('g%d' % i for i in range(3000))\n"
        "assert mod.evaluate_expression(gpr, ['g2999']) is True\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-2000:]


# -----------------------------------------------------------------------------
# U-*: confirmed uncertain items
# -----------------------------------------------------------------------------
def test_u_sbml_write_over_existing_non_fbc_file(tmp_path):
    import libsbml

    from mewpy.io import Engines, Writer, read_sbml, write_model

    path = str(tmp_path / "l2.xml")
    doc = libsbml.SBMLDocument(2, 4)
    doc.createModel().setId("x")
    libsbml.writeSBMLToFile(doc, path)
    write_model(Writer(Engines.MetabolicSBML, io=path, model=_read_ecoli()), warnings=False)
    assert len(read_sbml(path, regulatory=False, warnings=False).reactions) == 95


def test_u_minimal_medium_validation_infeasible_warns():
    from reframed import load_cbmodel

    from mewpy.cobra.medium import validate_solution
    from mewpy.simulation import get_simulator

    model = load_cbmodel(EC_CORE)
    exchanges = get_simulator(model).get_exchange_reactions()
    with pytest.warns(UserWarning, match="invalid"):
        validate_solution(model, set(), exchanges, -1, 1, 10)


# -----------------------------------------------------------------------------
# B-CSV: a self-regulated target was built as two objects (target and regulator) with the same id
# -----------------------------------------------------------------------------
def test_boolean_csv_self_regulated_target_is_one_object(tmp_path):
    from mewpy.io import read_csv

    csv = tmp_path / "self.csv"
    csv.write_text("b0008,b0008\nb0114,NOT b0008\n")
    model = read_csv(str(csv), warnings=False)

    target = model.targets["b0008"]
    assert model.regulators["b0008"] is target
    assert model.interactions["b0008_interaction"].target is target
    assert model.interactions["b0008_interaction"].regulators["b0008"] is target
    assert model.interactions["b0114_interaction"].regulators["b0008"] is target
    assert target.is_target() and target.is_regulator()


def test_target_regulator_csv_self_regulated_target_is_one_object(tmp_path):
    from mewpy.io import read_csv

    csv = tmp_path / "tr_self.csv"
    csv.write_text("g1,g1\ng2,g1\n")
    model = read_csv(str(csv), target_regulator=True, warnings=False)

    target = model.targets["g1"]
    assert model.regulators["g1"] is target
    assert model.interactions["g1_interaction"].regulators["g1"] is target
    assert model.interactions["g2_interaction"].regulators["g1"] is target


def test_co_expression_csv_self_regulated_target_is_one_object(tmp_path):
    from mewpy.io import read_csv

    csv = tmp_path / "co_self.csv"
    csv.write_text("g1,g1,\ng2,g1,\n")
    model = read_csv(str(csv), co_expression=True, warnings=False)

    target = model.targets["g1"]
    assert model.regulators["g1"] is target
    assert model.interactions["g1_interaction"].regulators["g1"] is target
    assert model.interactions["g2_interaction"].regulators["g1"] is target


# -----------------------------------------------------------------------------
# N-1: convert_to_irreversible printed a progress bar that could not be turned off
# -----------------------------------------------------------------------------
def test_convert_to_irreversible_is_quiet_by_default(capsys):
    from reframed import load_cbmodel

    from mewpy.cobra.util import convert_to_irreversible

    sim, mapping = convert_to_irreversible(load_cbmodel(EC_CORE))
    assert "Converting to irreversible" not in capsys.readouterr().err
    assert mapping and all(rev == f"{r_id}_REV" for r_id, rev in mapping.items())

    convert_to_irreversible(load_cbmodel(EC_CORE), verbose=True)
    assert "Converting to irreversible" in capsys.readouterr().err


# -----------------------------------------------------------------------------
# N-2: ODEModel.get_ode changed numpy's global floating-point error handling
# -----------------------------------------------------------------------------
def test_get_ode_keeps_numpy_error_state(tmp_path):
    import numpy as np

    from mewpy.io.sbml import load_ODEModel

    model = load_ODEModel(_kinetic_sbml(tmp_path / "div.xml", ["k / A"]))
    previous = np.seterr(divide="warn", invalid="warn")
    try:
        ode = model.get_ode()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            dxdt = ode(0, np.array([0.0, 1.0]))
        assert np.isinf(dxdt[0])
        # divisions by zero in the rate laws are silenced while the ODE is evaluated...
        assert not [w for w in caught if issubclass(w.category, RuntimeWarning)]
        # ...but numpy's global error handling is left as it was
        assert np.geterr()["divide"] == "warn" and np.geterr()["invalid"] == "warn"
    finally:
        np.seterr(**previous)
