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


def test_a_product_that_must_be_secreted_keeps_its_bounds():
    """A tracked compound without kinetics whose lower bound is positive (forced secretion) keeps that bound"""
    from mewpy.simulation.dfba import DynamicFBA, MichaelisMenten

    sim = _simulator("cobra", envcond={})  # aerobic: no acetate overflow unless forced
    dfba = DynamicFBA(
        sim,
        0.01,
        {"EX_glc__D_e": 10.0, "EX_ac_e": 0.0},
        uptake={"EX_glc__D_e": MichaelisMenten(10, 0.015)},
        constraints={"EX_ac_e": (0.5, 1000)},
    )
    result = dfba.simulate(t_end=1, dt=0.1)
    assert min(result.exchange_fluxes["EX_ac_e"]) >= 0.5 - 1e-9


# ---------------------------------------------------------------------------------------------------------------------
# Community dynamic FBA
# ---------------------------------------------------------------------------------------------------------------------
GLC_KINETICS = {"EX_glc__D_e": (10, 0.015)}
AC_KINETICS = {"EX_ac_e": (3, 0.5)}


def _kinetics(*kinetics):
    from mewpy.simulation.dfba import MichaelisMenten

    return {r: MichaelisMenten(*p) for k in kinetics for r, p in k.items()}


def _cross_feeding(backends=("cobra", "cobra"), t_end=12):
    """A producer that takes up glucose and secretes acetate, and a specialist (no glucose transport) that grows on
    that acetate"""
    from mewpy.simulation.dfba import CommunityDynamicFBA

    community = CommunityDynamicFBA(
        {"producer": _simulator(backends[0]), "specialist": _simulator(backends[1])},
        initial_biomass={"producer": 0.01, "specialist": 0.01},
        initial_concentrations={"EX_glc__D_e": 10.0, "EX_ac_e": 0.0},
        uptake={"producer": _kinetics(GLC_KINETICS), "specialist": _kinetics(AC_KINETICS)},
        constraints={"producer": {"EX_ac_e": (0, 1000)}, "specialist": {"GLCpts": 0}},
    )
    return community, community.simulate(t_end=t_end, dt=0.1)


def _halves():
    """Two identical members with half the biomass each, competing for glucose and acetate"""
    from mewpy.simulation.dfba import CommunityDynamicFBA

    members = {"a": _simulator("cobra"), "b": _simulator("cobra")}
    community = CommunityDynamicFBA(
        members,
        initial_biomass={"a": 0.005, "b": 0.005},
        initial_concentrations={"EX_glc__D_e": 10.0, "EX_ac_e": 0.0},
        uptake={m: _kinetics(GLC_KINETICS, AC_KINETICS) for m in members},
    )
    return community, community.simulate(t_end=12, dt=0.1)


def test_a_one_member_community_is_a_culture():
    """On REFRAMED (prefixed ids), a community of one member gives the DynamicFBA trajectory; the community keys
    compounds by the unprefixed id"""
    from mewpy.simulation.dfba import CommunityDynamicFBA

    dfba, single = _diauxie("reframed")
    community = CommunityDynamicFBA(
        {"ec": _simulator("reframed")},
        initial_biomass={"ec": 0.01},
        initial_concentrations={"EX_glc__D_e": 10.0, "EX_ac_e": 0.0},
        uptake={"ec": _kinetics(GLC_KINETICS, AC_KINETICS)},
    ).simulate(t_end=12, dt=0.1)
    assert community.biomass["ec"] == single.biomass
    assert community.concentrations == {
        dfba.simulator.to_canonical(r, "reaction"): c for r, c in single.concentrations.items()
    }
    assert community.status.startswith(single.status.split(":")[0])


def test_cross_feeding():
    """The specialist does not grow until the producer has secreted acetate; then it grows on it, and the producer
    stops growing when glucose runs out"""
    _, result = _cross_feeding(t_end=24)
    specialist, producer = result.growth_rate["specialist"], result.growth_rate["producer"]
    assert specialist[0] == 0  # no acetate yet: no feasible solution (maintenance), no growth
    assert result.biomass["specialist"][-1] > 1.2 * result.biomass["specialist"][0]
    assert all(v >= 0 for v in result.exchange_fluxes["producer"]["EX_ac_e"])  # the producer only secretes acetate
    assert all(v <= 0 for v in result.exchange_fluxes["specialist"]["EX_ac_e"])
    glc = result.concentrations["EX_glc__D_e"]
    assert glc[-1] < 1e-3 and producer[-1] == 0
    assert max(result.concentrations["EX_ac_e"]) > result.concentrations["EX_ac_e"][-1]


@pytest.mark.parametrize("scenario", [_halves, _cross_feeding])
def test_community_consumption_matches_the_medium(scenario):
    """Summed over members and steps, v X0 (exp(mu dt) - 1) / mu is the change in the medium: when the members
    together would take more of a compound than what is left, they share it out instead of overdrawing it"""
    _, result = scenario()
    dt = 0.1
    for c, conc in result.concentrations.items():
        exchanged = 0.0
        for m, fluxes in result.exchange_fluxes.items():
            for x, mu, v in zip(result.biomass[m], result.growth_rate[m], fluxes.get(c, [])):
                exchanged += v * (x * (math.exp(mu * dt) - 1) / mu if mu > 1e-9 else x * dt)
        assert conc[-1] == pytest.approx(conc[0] + exchanged, abs=1e-9), c
        assert min(conc) >= 0


def test_two_halves_grow_like_one_culture():
    """Two identical members with half the biomass each grow as one culture; at depletion the shares leave a
    different residue than a single culture does, within 0.1% of the biomass"""
    _, single = _diauxie("cobra")
    _, halves = _halves()
    total = [a + b for a, b in zip(halves.biomass["a"], halves.biomass["b"])]
    assert halves.biomass["a"] == pytest.approx(halves.biomass["b"], rel=1e-9)
    assert total == pytest.approx(single.biomass, rel=1e-3)


def test_members_on_different_backends():
    """A COBRApy producer and a REFRAMED specialist give the trajectory of the all-COBRApy community"""
    _, cobra_only = _cross_feeding(("cobra", "cobra"))
    _, mixed = _cross_feeding(("cobra", "reframed"))
    assert mixed.concentrations.keys() == cobra_only.concentrations.keys()
    for m in cobra_only.biomass:
        assert mixed.biomass[m] == pytest.approx(cobra_only.biomass[m], rel=1e-5), m
    for c in cobra_only.concentrations:
        assert mixed.concentrations[c] == pytest.approx(cobra_only.concentrations[c], rel=1e-5, abs=1e-7), c


def test_members_from_a_community_model():
    from mewpy.model import CommunityModel
    from mewpy.simulation.dfba import CommunityDynamicFBA

    producer, specialist = _simulator("cobra"), _simulator("cobra")
    producer.id, specialist.id = "producer", "specialist"
    community = CommunityModel([producer, specialist], merge_biomasses=False, flavor="cobra", verbose=False)
    dfba = CommunityDynamicFBA(
        community,
        initial_biomass={"producer": 0.01, "specialist": 0.01},
        initial_concentrations={"EX_glc__D_e": 10.0, "EX_ac_e": 0.0},
        uptake={"producer": _kinetics(GLC_KINETICS), "specialist": _kinetics(AC_KINETICS)},
        constraints={"producer": {"EX_ac_e": (0, 1000)}, "specialist": {"GLCpts": 0}},
    )
    _, expected = _cross_feeding()
    result = dfba.simulate(t_end=12, dt=0.1)
    assert set(dfba.simulators) == {"producer", "specialist"}
    for m in expected.biomass:
        assert result.biomass[m] == pytest.approx(expected.biomass[m], rel=1e-9), m


def test_community_stops_when_no_member_can_grow():
    """Glucose only, and both members need it: once it runs out no member has a feasible solution"""
    from mewpy.simulation.dfba import CommunityDynamicFBA

    members = {"a": _simulator("cobra"), "b": _simulator("cobra")}
    result = CommunityDynamicFBA(
        members,
        initial_biomass={"a": 0.01, "b": 0.02},
        initial_concentrations={"EX_glc__D_e": 5.0},
        uptake={m: _kinetics(GLC_KINETICS) for m in members},
    ).simulate(t_end=20, dt=0.1)
    assert result.status.startswith("stopped at t = ") and "no member has a feasible solution" in result.status
    # same kinetics and model: the biomass ratio stays 1:2
    assert result.biomass["b"][-1] / result.biomass["a"][-1] == pytest.approx(2, rel=1e-6)


def test_community_input_errors():
    from mewpy.simulation.dfba import CommunityDynamicFBA

    members = {"a": _simulator("cobra"), "b": _simulator("cobra")}
    conc = {"EX_glc__D_e": 1.0}
    with pytest.raises(ValueError, match="Unknown members in initial_biomass"):
        CommunityDynamicFBA(members, {"a": 1, "b": 1, "c": 1}, conc)
    with pytest.raises(ValueError, match="no value for the members"):
        CommunityDynamicFBA(members, {"a": 1}, conc)
    with pytest.raises(ValueError, match="must be positive"):
        CommunityDynamicFBA(members, {"a": 1, "b": 0}, conc)
    with pytest.raises(ValueError, match="Not an exchange reaction of any member"):
        CommunityDynamicFBA(members, {"a": 1, "b": 1}, {"PGI": 1.0})
    with pytest.raises(ValueError, match="Uptake kinetics of b"):
        CommunityDynamicFBA(members, {"a": 1, "b": 1}, conc, uptake={"b": _kinetics(AC_KINETICS)})
    with pytest.raises(ValueError, match="Unknown members in constraints"):
        CommunityDynamicFBA(members, {"a": 1, "b": 1}, conc, constraints={"z": {}})


def test_plots():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    _, single = _diauxie("cobra", t_end=1)
    _, community = _cross_feeding(t_end=1)
    for result, n_biomass in ((single, 1), (community, 2)):
        ax = result.plot()
        right = [a for a in ax.figure.axes if a is not ax][0]
        assert len(ax.lines) == n_biomass and len(right.lines) == len(result.concentrations)
