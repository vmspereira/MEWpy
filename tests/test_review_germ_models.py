"""
Regression tests for the GERM models/variables/solution review findings (F1-F12 and minor items).
"""

import copy
import logging
import pickle
import warnings
from pathlib import Path

import pytest

MODELS = Path(__file__).parent.parent / "examples" / "models" / "germ"
ECOLI = str(MODELS / "e_coli_core.xml")
ECOLI_TRN = str(MODELS / "e_coli_core_trn.csv")
TRN_KW = dict(sep=",", id_col=0, rule_col=2, aliases_cols=[1], header=0)


def _legacy():
    from mewpy.io import Engines, Reader, read_model

    gem = Reader(Engines.MetabolicSBML, ECOLI)
    trn = Reader(Engines.BooleanRegulatoryCSV, ECOLI_TRN, **TRN_KW)
    return read_model(gem, trn, warnings=False)


def _extension():
    from mewpy.germ.models import RegulatoryExtension

    return RegulatoryExtension.from_sbml(ECOLI, ECOLI_TRN, **TRN_KW)


@pytest.fixture
def model():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _legacy()


@pytest.fixture
def ext():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _extension()


def _expr(rule, variables):
    from mewpy.germ.algebra import Expression, parse_expression

    return Expression(parse_expression(rule), variables)


# --------------------------------------------------------------------------------------------------
# F1 - bound setters
# --------------------------------------------------------------------------------------------------
def test_f1_bound_setters_accept_zero(model):
    rxn = model.reactions["EX_glc__D_e"]
    rxn.lower_bound = 0
    assert rxn.bounds == (0, 1000.0)
    rxn.upper_bound = 0
    assert rxn.bounds == (0, 0)


# --------------------------------------------------------------------------------------------------
# F2 - observer pattern
# --------------------------------------------------------------------------------------------------
def test_f2_attached_analysis_is_resynchronized_after_model_change(model):
    from mewpy.germ.analysis import SRFBA

    srfba = SRFBA(model).build()
    model.attach(srfba)
    assert srfba.optimize().objective_value > 0.5

    model.reactions["Biomass_Ecoli_core"].ko()
    assert srfba.optimize().objective_value == pytest.approx(0.0, abs=1e-6)


def test_f2_regulatory_extension_supports_attach(ext):
    class Observer:
        updates = 0

        def update(self):
            Observer.updates += 1

    observer = Observer()
    ext.attach(observer)
    assert observer in ext.simulators
    ext.notify()
    assert Observer.updates == 1
    ext.detach(observer)
    assert observer not in ext.simulators


def test_f2_attach_true_documented_workflow(model):
    from mewpy.germ.analysis import SRFBA

    srfba = SRFBA(model, attach=True).build()
    assert srfba in model.simulators
    with model:
        model.reactions["Biomass_Ecoli_core"].ko()
        assert srfba.optimize().objective_value == pytest.approx(0.0, abs=1e-6)
    assert srfba.optimize().objective_value > 0.5


# --------------------------------------------------------------------------------------------------
# F3 - TRN id reconciliation in RegulatoryExtension
# --------------------------------------------------------------------------------------------------
def test_f3_trn_ids_are_reconciled_with_simulator_ids(ext):
    genes = set(ext.genes)
    assert len(set(ext.targets) & genes) == 137
    assert ext.id_map["b0351"] == "G_b0351"
    assert ext.id_map["PFK"] == "R_PFK"
    assert "R_PFK" in ext.regulators and "PFK" not in ext.regulators
    # symbols inside regulatory expressions follow the renamed variables
    interaction = ext.get_interaction("surplusPYR_interaction")
    assert "R_PFK" in interaction.regulators
    for expression in interaction.yield_expressions():
        assert set(expression.symbols) == set(expression.variables)
    # TRN-only regulators keep their ids
    assert "surplusFDP" in ext.regulators


def test_f3_regulator_deletion_reaches_metabolism(ext):
    from mewpy.germ.analysis import isingle_regulator_deletion

    # the same regulator deletions reduce SRFBA growth as with the legacy model (pre-fix: none did)
    legacy = isingle_regulator_deletion(_legacy(), method="srfba")["growth"].astype(float)
    growth = isingle_regulator_deletion(ext, method="srfba")["growth"].astype(float)
    original_ids = {new: old for old, new in ext.id_map.items()}
    hits = sorted(original_ids.get(reg_id, reg_id) for reg_id in growth[growth < 0.87].index)
    assert "PFK" in hits
    assert hits == sorted(legacy[legacy < 0.87].index)
    assert growth["R_PFK"] == pytest.approx(legacy["PFK"], abs=1e-6)


def test_f3_variables_added_later_are_reconciled(ext):
    from mewpy.germ.models import RegulatoryExtension
    from mewpy.germ.variables import Interaction, Regulator, Target

    empty = RegulatoryExtension(ext.simulator)
    pfk, pdhr = Regulator("PFK"), Regulator("b0113")
    events = {1.0: _expr("PFK & ~b0113", {"PFK": pfk, "b0113": pdhr})}
    interaction = Interaction("b0351_interaction", target=Target("b0351"), regulatory_events=events)

    empty.add_interaction(interaction)
    assert interaction.target.id == "G_b0351"
    assert sorted(interaction.regulators) == ["R_PFK", "b0113"]
    expression = interaction.regulatory_events[1.0]
    assert set(expression.symbols) == set(expression.variables) == {"R_PFK", "b0113"}
    assert expression.evaluate(values={"R_PFK": 1, "b0113": 0}) == 1

    empty.add_target(Target("b0114"))
    empty.add_regulator(Regulator("glc__D_e"))
    assert "G_b0114" in empty.targets and "M_glc__D_e" in empty.regulators
    expected = {"b0351": "G_b0351", "PFK": "R_PFK", "b0114": "G_b0114", "glc__D_e": "M_glc__D_e"}
    assert empty.id_map == expected


def test_f3_self_regulated_target_is_renamed_once(ext, tmp_path):
    from mewpy.germ.models import RegulatoryExtension
    from mewpy.io import Engines, Reader, read_model

    path = tmp_path / "self_regulated_trn.csv"
    path.write_text("b0008,b0008\nb0114,NOT b0008\n")
    network = read_model(Reader(Engines.BooleanRegulatoryCSV, str(path), sep=","), warnings=False)
    small = RegulatoryExtension(ext.simulator, network)

    self_regulated = small.get_interaction("b0008_interaction")
    assert self_regulated.target.id == "G_b0008"
    assert set(self_regulated.regulators) == {"G_b0008"}
    assert set(small.get_interaction("b0114_interaction").regulators) == {"G_b0008"}
    assert "G_b0008" in small.targets and "G_b0008" in small.regulators
    assert "b0008" not in small.targets and "b0008" not in small.regulators
    for interaction in small.interactions.values():
        for expression in interaction.yield_expressions():
            assert set(expression.symbols) == set(expression.variables)
    assert small.id_map == {"b0008": "G_b0008", "b0114": "G_b0114"}


def test_f3_reconciliation_follows_the_simulator_prefixes():
    from mewpy.germ.models import RegulatoryExtension, RegulatoryModel
    from mewpy.germ.variables import Interaction, Regulator, Target
    from mewpy.simulation.simulation import Simulator

    class PrefixedSimulator(Simulator):
        """A backend whose native ids use non-SBML prefixes."""

        _r_prefix, _m_prefix, _g_prefix = "rxn_", "met_", "gene_"
        id = "prefixed"
        reactions = ["rxn_PFK"]
        metabolites = ["met_glc__D_e"]
        genes = ["gene_b0351", "gene_b0114"]

        def __init__(self):
            self.model = None

        def FVA(self, *args, **kwargs):
            raise NotImplementedError

        def get_reaction_bounds(self, *args, **kwargs):
            raise NotImplementedError

        def metabolite_reaction_lookup(self, *args, **kwargs):
            raise NotImplementedError

    pfk, glc = Regulator("PFK"), Regulator("glc__D_e")
    events = {1.0: _expr("PFK & glc__D_e", {"PFK": pfk, "glc__D_e": glc})}
    interaction = Interaction("b0351_interaction", target=Target("b0351"), regulatory_events=events)
    network = RegulatoryModel("net", interactions={"b0351_interaction": interaction})

    ext = RegulatoryExtension(PrefixedSimulator(), network)
    assert interaction.target.id == "gene_b0351"
    assert sorted(interaction.regulators) == ["met_glc__D_e", "rxn_PFK"]
    ext.add_target(Target("b0114"))
    expected = {
        "b0351": "gene_b0351",
        "PFK": "rxn_PFK",
        "glc__D_e": "met_glc__D_e",
        "b0114": "gene_b0114",
    }
    assert ext.id_map == expected


def test_f3_warns_when_trn_does_not_overlap_metabolism(caplog):
    from mewpy.germ.models import RegulatoryExtension, RegulatoryModel
    from mewpy.germ.variables import Interaction, Regulator, Target

    base = _extension()
    reg = Regulator("unknown_regulator")
    interaction = Interaction(
        "unknown_target_interaction",
        target=Target("unknown_target"),
        regulatory_events={1.0: _expr("unknown_regulator", {"unknown_regulator": reg})},
    )
    network = RegulatoryModel("net")
    network.add(interaction)
    with caplog.at_level(logging.WARNING, logger="mewpy.germ.models.regulatory_extension"):
        RegulatoryExtension(base.simulator, network)
    assert "no regulatory target matches a metabolic gene" in caplog.text


# --------------------------------------------------------------------------------------------------
# F4 - shallow copy must not share containers
# --------------------------------------------------------------------------------------------------
def test_f4_copy_does_not_share_containers(model):
    pfk = model.reactions["PFK"]
    shallow = model.copy()
    assert shallow.reactions["PFK"] is pfk
    shallow.remove(shallow.reactions["PFK"])
    assert "PFK" in model.reactions
    assert pfk.model is model


# --------------------------------------------------------------------------------------------------
# F5 - undo / redo
# --------------------------------------------------------------------------------------------------
def test_f5_context_add_of_existing_variable_keeps_it(model):
    pfk = model.reactions["PFK"]
    with model:
        model.add(pfk)
    assert "PFK" in model.reactions


def test_f5_context_remove_of_foreign_variable_does_not_add_it(model):
    from mewpy.germ.variables import Reaction

    foreign = Reaction("FOREIGN")
    with model:
        model.remove(foreign)
    assert "FOREIGN" not in model.reactions


def test_f5_add_regulatory_event_undo_restores_previous_event(model):
    interaction = model.interactions["b0080_interaction"]
    before = {k: v.to_string() for k, v in interaction.regulatory_events.items()}
    reg = model.regulators["b3261"]
    with model:
        interaction.add_regulatory_event(1, _expr("b3261", {"b3261": reg}))
        assert interaction.regulatory_events[1].to_string() == "b3261"
    assert {k: v.to_string() for k, v in interaction.regulatory_events.items()} == before


def test_f5_add_gpr_undo_restores_previous_gpr(model):
    rxn = model.reactions["PDH"]
    before = rxn.gene_protein_reaction_rule
    gene = model.genes["b0116"]
    with model:
        rxn.add_gpr(_expr("b0116", {"b0116": gene}))
        assert rxn.gene_protein_reaction_rule == "b0116"
    assert rxn.gene_protein_reaction_rule == before
    assert set(rxn.genes) == {"b0114", "b0115", "b0116"}


def test_f5_add_target_undo_restores_previous_target():
    from mewpy.germ.models import RegulatoryModel
    from mewpy.germ.variables import Interaction, Regulator, Target

    a, old_target, new_target = Regulator("A"), Target("T"), Target("T2")
    interaction = Interaction("i", target=old_target, regulatory_events={1.0: _expr("A", {"A": a})})
    network = RegulatoryModel("rm", interactions={"i": interaction})
    with network:
        interaction.add_target(new_target)
        assert old_target.interaction is None
    assert interaction.target is old_target
    assert old_target.interaction is interaction
    assert new_target.interaction is None
    assert sorted(network.targets) == ["T"]


def test_f5_undo_redo_undo(model):
    rxn = model.reactions["PFK"]
    old = rxn.bounds
    rxn.bounds = (0, 0)
    model.undo()
    assert rxn.bounds == old
    model.redo()
    assert rxn.bounds == (0, 0)
    model.undo()
    assert rxn.bounds == old


def test_f5_new_command_clears_redo_stack(model):
    rxn = model.reactions["PFK"]
    rxn.bounds = (0, 0)
    model.undo()
    assert len(model.history.redo_able_commands) == 1
    rxn.bounds = (1, 2)
    assert len(model.history.redo_able_commands) == 0


# --------------------------------------------------------------------------------------------------
# F6 - deepcopy must start with a clean history
# --------------------------------------------------------------------------------------------------
def test_f6_deepcopy_history_is_clean(model):
    clone = model.deepcopy()
    assert len(clone.history.undo_able_commands) == 0
    clone.reactions["PFK"].bounds = (0, 0)
    clone.reset()
    assert len(clone.reactions) == len(model.reactions)
    assert len(clone.interactions) == len(model.interactions)
    assert clone.reactions["PFK"].bounds == model.reactions["PFK"].bounds


def test_f6_target_interaction_setter_accepts_none():
    from mewpy.germ.variables import Target

    target = Target("t")
    target.interaction = None
    assert target.interaction is None


# --------------------------------------------------------------------------------------------------
# F7 - pickled variables stay attached to the model
# --------------------------------------------------------------------------------------------------
def test_f7_unpickled_variables_belong_to_unpickled_model(model):
    clone = pickle.loads(pickle.dumps(model))
    pfk = clone.reactions["PFK"]
    assert pfk.model is clone
    assert clone.genes["b3916"].model is clone
    old = pfk.bounds
    with clone:
        pfk.ko()
    assert pfk.bounds == old


# --------------------------------------------------------------------------------------------------
# F8 - RegulatoryExtension.get
# --------------------------------------------------------------------------------------------------
def test_f8_regulatory_extension_get(ext):
    assert ext.get("G_b0351") is ext.get_target("G_b0351")
    assert ext.get("surplusFDP") is ext.get_regulator("surplusFDP")
    assert ext.get("b0080_interaction") is ext.get_interaction("b0080_interaction")
    assert ext.get("missing", "default") == "default"


def test_f8_analysis_runs_on_regulatory_extension(ext):
    from mewpy.germ.analysis import find_conflicts
    from mewpy.germ.analysis.metabolic_analysis import single_gene_deletion

    legacy = single_gene_deletion(_legacy())["growth"].astype(float)
    growth = single_gene_deletion(ext)["growth"].astype(float)
    essential = sorted(gene[2:] for gene in growth[growth < 1e-6].index)
    assert essential == sorted(legacy[legacy < 1e-6].index)

    gene_conflicts, reaction_conflicts = find_conflicts(ext)
    assert gene_conflicts.empty and reaction_conflicts.empty


# --------------------------------------------------------------------------------------------------
# F9 - RegulatoryExtension save/load
# --------------------------------------------------------------------------------------------------
def test_f9_save_load_roundtrip(ext, tmp_path):
    from mewpy.germ.analysis import SRFBA
    from mewpy.germ.models import RegulatoryExtension

    path = tmp_path / "ext.json"
    ext.save(str(path))
    loaded = RegulatoryExtension.load(str(path), ext.simulator)

    interaction = loaded.get_interaction("b0080_interaction")
    assert all(isinstance(k, (int, float)) for k in interaction.regulatory_events)
    regulator_id = next(iter(interaction.regulators))
    assert interaction.regulators[regulator_id] is loaded.get_regulator(regulator_id)
    assert interaction.target is loaded.get_target(interaction.target.id)
    expected_types = {"regulator", "target"}
    assert loaded.get_regulator("b0113").types == ext.get_regulator("b0113").types == expected_types

    expected = SRFBA(ext).build().optimize().objective_value
    assert SRFBA(loaded).build().optimize().objective_value == pytest.approx(expected)


# --------------------------------------------------------------------------------------------------
# F10 - SimulatorBasedMetabolicModel
# --------------------------------------------------------------------------------------------------
@pytest.fixture
def sbm():
    from mewpy.germ.models.unified_factory import load_reframed_model

    return load_reframed_model(ECOLI)


def test_f10_simulator_model_pickle_and_copy(sbm):
    from mewpy.simulation import get_simulator

    clone = pickle.loads(pickle.dumps(sbm))
    assert set(clone.reactions) == set(sbm.reactions)
    assert set(sbm.copy().reactions) == set(sbm.reactions)
    assert set(copy.deepcopy(sbm).reactions) == set(sbm.reactions)
    pickle.loads(pickle.dumps(get_simulator(sbm)))


def test_f10_simulator_model_add_remove_raise(sbm):
    pfk = sbm.reactions["R_PFK"]
    with pytest.raises(NotImplementedError):
        sbm.remove(pfk)
    with pytest.raises(NotImplementedError):
        sbm.add(pfk)


def test_f10_simulator_model_bounds_reach_simulator(sbm):
    from mewpy.simulation import get_simulator

    pfk = sbm.reactions["R_PFK"]
    pfk.bounds = (0, 0)
    assert sbm.simulator.get_reaction_bounds("R_PFK") == (0, 0)
    assert get_simulator(sbm).simulate().objective_value == pytest.approx(0.704, abs=1e-3)


def test_f10_gene_extraction_uses_parsed_gpr():
    from mewpy.germ.models.simulator_model import SimulatorBasedMetabolicModel

    extract = SimulatorBasedMetabolicModel._extract_genes_from_gpr
    assert extract(None, "YAL001C-A and b0002") == {"YAL001C-A", "b0002"}
    assert extract(None, "10005_AT1 or 2593_AT1") == {"10005_AT1", "2593_AT1"}


# --------------------------------------------------------------------------------------------------
# F11 - shared regulator across regulatory events
# --------------------------------------------------------------------------------------------------
def test_f11_shared_regulator_events():
    from mewpy.germ.variables import Interaction, Regulator, Target

    a, b = Regulator("A"), Regulator("B")
    interaction = Interaction(
        "i",
        target=Target("T"),
        regulatory_events={1.0: _expr("A & B", {"A": a, "B": b}), 0.0: _expr("~A", {"A": a})},
    )
    interaction.remove_regulatory_event(1.0)
    assert "i" in a.interactions
    assert "i" not in b.interactions
    interaction.remove_regulatory_event(0.0)
    assert "i" not in a.interactions

    events = {1.0: _expr("A & B", {"A": a, "B": b}), 0.0: _expr("~A", {"A": a})}
    interaction.regulatory_events = events
    interaction.regulatory_events = {1.0: _expr("B", {"B": b})}
    assert "i" not in a.interactions and "i" in b.interactions


# --------------------------------------------------------------------------------------------------
# F12 - constructors and factories
# --------------------------------------------------------------------------------------------------
def test_f12_load_integrated_model():
    from mewpy.germ.models import load_integrated_model

    integrated = load_integrated_model(
        ECOLI, ECOLI_TRN, backend="reframed", regulatory_kwargs=TRN_KW
    )
    assert len(integrated.interactions) == 159
    assert "G_b0351" in integrated.targets


def test_f12_regulatory_model_constructor_keeps_regulators_and_targets():
    from mewpy.germ.models import RegulatoryModel
    from mewpy.germ.variables import Interaction, Regulator, Target

    a, b = Regulator("A"), Regulator("B")
    events = {1.0: _expr("A & B", {"A": a, "B": b})}
    interaction = Interaction("i", target=Target("T"), regulatory_events=events)
    network = RegulatoryModel("rm", interactions={"i": interaction})
    assert sorted(network.regulators) == ["A", "B"]
    assert sorted(network.targets) == ["T"]


def test_f12_model_constructor_history_is_clean():
    from mewpy.germ.models import Model
    from mewpy.germ.variables import Metabolite, Reaction

    met = Metabolite("a", compartment="c")
    rxn = Reaction("r", stoichiometry={met: -1})
    metabolic = Model.from_types(["metabolic"], identifier="x", reactions={"r": rxn})
    assert len(metabolic.history.undo_able_commands) == 0
    metabolic.reset()
    assert list(metabolic.reactions) == ["r"]


def test_f12_parsing_with_model_creates_unknown_variables():
    from mewpy.germ.models import RegulatoryModel
    from mewpy.germ.variables import Interaction, Reaction, Target

    network = RegulatoryModel("rm")
    interaction = Interaction.from_string("j", rule="A & NEW", target=Target("T2"), model=network)
    assert set(interaction.regulators) == {"A", "NEW"}
    assert {"A", "NEW"} <= set(network.regulators)
    rxn = Reaction.from_gpr_string("R", rule="g1 and g2", model=network)
    assert set(rxn.genes) == {"g1", "g2"}


# --------------------------------------------------------------------------------------------------
# Minor items
# --------------------------------------------------------------------------------------------------
def test_minor_solution_str():
    from mewpy.germ.solution import DynamicSolution, KOSolution, MultiSolution

    class Sol:
        def __init__(self, method):
            self.method = method

    assert str(MultiSolution(Sol("fba"), Sol("pfba"))) == "MultiSolution: fba,pfba"
    assert str(DynamicSolution(Sol("a"), Sol("b"))) == "DynamicSolution: t_0,t_1"
    assert str(KOSolution({"g1": Sol("a")})) == "KOSolution: ko_g1"


def test_minor_variable_clean_context(model):
    rxn = model.reactions["PFK"]
    with rxn:
        rxn.bounds = (0, 0)
        rxn.clean_context()
    assert rxn.bounds == (0, 0)


def test_minor_compartment_names_survive_copies(model):
    assert model.compartments["c"] == "cytosol"
    assert pickle.loads(pickle.dumps(model)).compartments["c"] == "cytosol"
    assert model.copy().compartments["c"] == "cytosol"
    assert model.deepcopy().compartments["c"] == "cytosol"


def test_minor_compartment_names_survive_json_roundtrip(model, tmp_path):
    from mewpy.io import Engines, Reader, Writer, read_model, write_model

    path = str(tmp_path / "model.json")
    write_model(Writer(Engines.JSON, path, model), warnings=False)
    loaded = read_model(Reader(Engines.JSON, path), warnings=False)
    assert loaded.compartments["c"] == "cytosol"
    assert loaded.compartments == model.compartments
    assert len(loaded.reactions) == len(model.reactions)


def test_minor_reaction_deepcopy_keeps_metabolite_attributes(model):
    pfk = model.reactions["PFK"]
    clone = pfk.deepcopy()
    assert clone.compartments == pfk.compartments == {"c"}
    formulas = {m.id: m.formula for m in pfk.yield_metabolites()}
    assert {m.id: m.formula for m in clone.yield_metabolites()} == formulas
