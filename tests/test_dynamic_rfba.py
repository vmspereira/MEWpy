"""Dynamic regulatory FBA: E. coli core with its regulatory network (Orth et al. 2010) in batch culture"""

from pathlib import Path

import pytest

GERM = Path(__file__).parent.parent / "examples" / "models" / "germ"
MODEL, TRN = str(GERM / "e_coli_core.xml"), str(GERM / "e_coli_core_trn.csv")
TRN_CSV = dict(regulatory_format="csv", sep=",", id_col=0, rule_col=2, aliases_cols=[1], header=0)
GLC, AC = "R_EX_glc__D_e", "R_EX_ac_e"
O2_LIMITED = {"EX_o2_e": (-15, 1000)}
FADR, ICLR, CRP = "b1187", "b4018", "b3357"


def _model(regulatory=True):
    from mewpy.germ.models import RegulatoryExtension

    if not regulatory:
        return RegulatoryExtension.from_sbml(MODEL, flavor="reframed")
    return RegulatoryExtension.from_sbml(MODEL, TRN, flavor="reframed", **TRN_CSV)


def _kinetics():
    from mewpy.simulation import MichaelisMenten

    return {"EX_glc__D_e": MichaelisMenten(10, 0.015), "EX_ac_e": MichaelisMenten(3, 0.5)}


def _diauxie(constraints=O2_LIMITED, t_end=12, **kwargs):
    """Glucose, then the acetate secreted while glucose lasts"""
    from mewpy.germ.analysis import DynamicRFBA

    dfba = DynamicRFBA(
        _model(),
        0.01,
        {"EX_glc__D_e": 10.0, "EX_ac_e": 0.0},
        uptake=_kinetics(),
        constraints=constraints,
        **kwargs,
    )
    return dfba, dfba.simulate(t_end=t_end, dt=0.1)


def _glucose_out(result, threshold=1e-3):
    """The first step whose glucose concentration is at most the threshold"""
    return next(i for i, c in enumerate(result.concentrations[GLC]) if c <= threshold)


def test_without_a_regulatory_network_it_is_dynamic_fba():
    from mewpy.germ.analysis import DynamicRFBA
    from mewpy.simulation import DynamicFBA

    conc = {"EX_glc__D_e": 10.0, "EX_ac_e": 0.0}
    plain = DynamicFBA(_model(False).simulator, 0.01, conc, uptake=_kinetics(), constraints=O2_LIMITED).simulate(12)
    regulated = DynamicRFBA(_model(False), 0.01, conc, uptake=_kinetics(), constraints=O2_LIMITED).simulate(12)
    assert regulated.biomass == plain.biomass
    assert regulated.concentrations == plain.concentrations
    assert all(not off for off in regulated.inactive)


def test_the_glyoxylate_shunt_opens_when_glucose_runs_out():
    """FadR (glucose present, or no acetate) activates IclR, which represses aceA (isocitrate lyase, ICL). Once
    glucose is gone, the shunt opens and the culture grows on the acetate it secreted"""
    _, result = _diauxie()
    out = _glucose_out(result)
    for i, (state, off) in enumerate(zip(result.regulatory_states, result.inactive)):
        glucose = i < out
        assert state[FADR] == state[ICLR] == (1 if glucose else 0), i
        assert ("R_ICL" in off) == glucose, i
    assert result.biomass[-1] > 1.03 * result.biomass[out]
    assert result.concentrations[AC][-1] < 0.5 * result.concentrations[AC][out]


def test_a_constitutive_fadr_blocks_growth_on_acetate():
    """FadR fixed on keeps the shunt closed: after glucose, acetate covers maintenance but not growth"""
    _, wild_type = _diauxie()
    _, fadr = _diauxie(initial_state={FADR: 1})
    out = _glucose_out(fadr)
    assert all("R_ICL" in off for off in fadr.inactive)
    assert all(mu == 0 for mu in fadr.growth_rate[out + 1 :])
    assert fadr.biomass[-1] < 1.001 * fadr.biomass[out] < wild_type.biomass[-1]


def test_a_gene_knockout_overrides_the_regulation():
    _, aceA = _diauxie(initial_state={"b4015": 0})
    assert all("R_ICL" in off for off in aceA.inactive)


def test_regulators_follow_the_environment_from_the_first_step():
    """The regulators that follow rules settle in every step: Crp, whose rule depends on glucose through CRPnoGLC,
    is off from the first step while glucose is present, and on once it runs out"""
    _, result = _diauxie()
    out = _glucose_out(result)
    assert result.regulatory_states[0][CRP] == 0
    assert result.regulatory_states[out][CRP] == 1


def test_reaction_regulators_take_the_flux_of_the_previous_step():
    _, result = _diauxie(t_end=1)
    biomass = "R_Biomass_Ecoli_core"
    states = [s[biomass] for s in result.regulatory_states]
    assert states[0] == 1  # the regulator's maximum coefficient
    assert states[1:] == pytest.approx(result.growth_rate[:-1], rel=1e-12)


def test_microaerobic_regulation_lowers_growth():
    """With some oxygen, the network keeps pyruvate formate-lyase (PFL) off, which the unregulated optimum uses: the
    first step is the FBA without the reactions the regulation switches off"""
    from mewpy.simulation import DynamicFBA

    microaerobic = {"EX_o2_e": (-5, 1000)}
    _, regulated = _diauxie(constraints=microaerobic, t_end=0.1)
    off = regulated.inactive[0]
    assert "R_PFL" in off
    conc = {"EX_glc__D_e": 10.0, "EX_ac_e": 0.0}

    def first_growth(constraints):
        return (
            DynamicFBA(_model(False).simulator, 0.01, conc, uptake=_kinetics(), constraints=constraints)
            .simulate(0.1)
            .growth_rate[0]
        )

    unregulated = first_growth(microaerobic)
    knocked_out = first_growth({**microaerobic, **{r: 0 for r in off}})
    assert regulated.growth_rate[0] == pytest.approx(knocked_out, rel=1e-9)
    assert regulated.growth_rate[0] < 0.97 * unregulated


@pytest.mark.parametrize("threshold", [1e-3, 1.0])
def test_a_metabolite_regulator_is_present_above_the_threshold(threshold):
    _, result = _diauxie(threshold=threshold)
    glucose = result.concentrations[GLC]
    for i, state in enumerate(result.regulatory_states):
        assert (state["M_glc__D_e"] > 0) == (glucose[i] > threshold), i


def test_input_errors():
    from mewpy.germ.analysis import DynamicRFBA

    with pytest.raises(ValueError, match="threshold"):
        DynamicRFBA(_model(), 0.01, {"EX_glc__D_e": 10.0}, threshold=-1)
