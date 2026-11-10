"""Dynamic FBA (static optimisation): integration, limits and backend independence."""

import math
from pathlib import Path

import pytest

EC_CORE = str(Path(__file__).parent.joinpath("data", "e_coli_core.xml"))
O2 = {"EX_o2_e": (-15, 1000)}


def _simulator(backend, envcond=O2):
    from mewpy.simulation import get_simulator

    if backend == "cobra":
        cobra = pytest.importorskip("cobra")
        model = cobra.io.read_sbml_model(EC_CORE)
    elif backend == "reframed":
        from reframed import load_cbmodel

        model = load_cbmodel(EC_CORE)
    else:
        from mewpy.io import read_sbml

        model = read_sbml(EC_CORE, regulatory=False, warnings=False)
    return get_simulator(model, envcond=envcond)


def _diauxie(backend, t_end=12):
    from mewpy.simulation.dfba import DynamicFBA, MichaelisMenten

    sim = _simulator(backend)
    dfba = DynamicFBA(
        sim,
        initial_biomass=0.01,
        initial_concentrations={"EX_glc__D_e": 10.0, "EX_ac_e": 0.0},
        uptake={"EX_glc__D_e": MichaelisMenten(vmax=10, km=0.015), "EX_ac_e": MichaelisMenten(vmax=3, km=0.5)},
    )
    return dfba, dfba.simulate(t_end=t_end, dt=0.1)


@pytest.mark.parametrize("backend", ["cobra", "reframed"])
def test_exponential_phase_is_exact(backend):
    """With the substrate in excess, uptake stays at vmax: biomass grows as X0 exp(mu t), with mu the FBA growth at
    that uptake, and glucose falls by q X0 (exp(mu t) - 1) / mu"""
    from mewpy.simulation.dfba import DynamicFBA, MichaelisMenten

    sim = _simulator(backend)
    glc = sim.resolve("EX_glc__D_e", "reaction")
    kinetics = MichaelisMenten(vmax=10, km=1e-9)
    dfba = DynamicFBA(sim, 0.01, {glc: 1e4}, uptake={glc: kinetics})
    reference = sim.simulate(method="pFBA", constraints={glc: (-10, 1000)})
    mu, q = reference.fluxes[dfba.biomass], -reference.fluxes[glc]
    result = dfba.simulate(t_end=2, dt=0.1)
    assert result.status == "completed"
    assert result.biomass[-1] == pytest.approx(0.01 * math.exp(mu * 2), rel=1e-6)
    consumed = q * 0.01 * (math.exp(mu * 2) - 1) / mu
    assert result.concentrations[glc][-1] == pytest.approx(1e4 - consumed, rel=1e-9)


@pytest.mark.parametrize("backend", ["cobra", "reframed"])
def test_diauxie_on_glucose_then_acetate(backend):
    """E. coli core under O2 limitation secretes acetate while glucose lasts, then grows on it; concentrations stay
    non-negative"""
    dfba, result = _diauxie(backend)
    df = result.dataframe
    glc = df[dfba.simulator.resolve("EX_glc__D_e", "reaction")]
    ac = df[dfba.simulator.resolve("EX_ac_e", "reaction")]
    assert (df.drop(columns="time") >= 0).all().all()
    t_glucose_out = df.time[glc.le(1e-9).idxmax()]
    peak = ac.idxmax()
    assert df.time[peak] == pytest.approx(t_glucose_out, abs=0.11)
    assert ac.iloc[-1] < 0.5 * ac[peak]
    assert df.biomass.iloc[-1] > df.biomass[peak]  # growth on acetate after glucose runs out


@pytest.mark.parametrize("backend", ["cobra", "reframed"])
def test_consumption_matches_the_medium(backend):
    """A step changes a concentration by v X0 (exp(mu dt) - 1) / mu; summed over the steps, this is the change in
    the medium, also in the step where glucose runs out (no uptake beyond what is left)"""
    dfba, result = _diauxie(backend)
    dt = 0.1
    exposure = [
        x * (math.exp(mu * dt) - 1) / mu if mu > 1e-9 else x * dt for x, mu in zip(result.biomass, result.growth_rate)
    ]
    for r_id, conc in result.concentrations.items():
        exchanged = sum(v * e for v, e in zip(result.exchange_fluxes[r_id], exposure))
        assert conc[-1] == pytest.approx(conc[0] + exchanged, abs=1e-9), r_id


def test_same_trajectory_on_every_backend():
    """The same call gives the same culture on COBRApy, REFRAMED and GERM (ids without the SBML prefix)"""
    trajectories = {}
    for backend in ("cobra", "reframed", "germ"):
        dfba, result = _diauxie(backend, t_end=8)
        sim = dfba.simulator
        trajectories[backend] = (
            result.time,
            result.biomass,
            {sim.to_canonical(r, "reaction"): c for r, c in result.concentrations.items()},
        )
    time, biomass, conc = trajectories["cobra"]
    for backend in ("reframed", "germ"):
        other_time, other_biomass, other_conc = trajectories[backend]
        assert other_time == time
        assert other_biomass == pytest.approx(biomass, rel=1e-5), backend
        for r_id in conc:
            assert other_conc[r_id] == pytest.approx(conc[r_id], rel=1e-5, abs=1e-7), (backend, r_id)


def test_growth_stops_when_the_substrate_runs_out():
    """Without a second substrate, the culture stops once glucose cannot sustain maintenance; the result says when"""
    from mewpy.simulation.dfba import DynamicFBA, MichaelisMenten

    sim = _simulator("cobra")
    dfba = DynamicFBA(sim, 0.01, {"EX_glc__D_e": 5.0}, uptake={"EX_glc__D_e": MichaelisMenten(10, 0.015)})
    result = dfba.simulate(t_end=20, dt=0.1)
    assert result.status.startswith("stopped at t = ")
    # what is left (< 0.02% of the 5 mM) cannot be taken up fast enough to cover maintenance
    assert result.concentrations["EX_glc__D_e"][-1] < 1e-3
    assert len(result.time) == len(result.biomass) == len(result.growth_rate) + 1


def test_ids_with_the_prefix_on_cobra():
    from mewpy.simulation.dfba import DynamicFBA, MichaelisMenten

    sim = _simulator("cobra")
    dfba = DynamicFBA(sim, 0.01, {"R_EX_glc__D_e": 10.0}, uptake={"R_EX_glc__D_e": MichaelisMenten(10, 0.015)})
    assert list(dfba.initial_concentrations) == ["EX_glc__D_e"]
    assert dfba.simulate(t_end=1, dt=0.1).status == "completed"


def test_input_errors():
    from mewpy.simulation.dfba import DynamicFBA, MichaelisMenten

    sim = _simulator("cobra")
    with pytest.raises(ValueError, match="Not exchange reactions"):
        DynamicFBA(sim, 0.01, {"PGI": 1.0})
    with pytest.raises(ValueError, match="without an initial concentration"):
        DynamicFBA(sim, 0.01, {"EX_glc__D_e": 1.0}, uptake={"EX_ac_e": MichaelisMenten(3, 0.5)})
    with pytest.raises(ValueError, match="initial_biomass"):
        DynamicFBA(sim, 0, {"EX_glc__D_e": 1.0})
