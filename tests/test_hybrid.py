"""Hybrid kinetic/constraint-based simulation: unit conversion and the constraints passed to the CB model."""

from pathlib import Path
from types import SimpleNamespace

import pytest

DATA = Path(__file__).parent.joinpath("data")
GDW = 564.0
TO_CB = 3600 / GDW  # mM/s -> mmol/gDW/h


class _FakeKineticSimulation:
    """Returns fixed kinetic rates; records the factors it was given."""

    fluxes = {}
    calls = []

    def __init__(self, model=None, t_points=None, timeout=None):
        pass

    def simulate(self, parameters=None, initcon=None, factors=None):
        _FakeKineticSimulation.calls.append(factors)
        return SimpleNamespace(fluxes=dict(_FakeKineticSimulation.fluxes))


class _RecordingSimulator:
    """Records the constraints of each simulate() call."""

    def __init__(self, reactions=(), objective=None):
        self.model = SimpleNamespace(reactions=list(reactions))
        self.objective = objective or {}
        self.constraints = []

    def get_reaction_bounds(self, r_id):
        return (-1000, 1000)

    def simulate(self, objective=None, method=None, constraints=None, **kwargs):
        self.constraints.append(dict(constraints or {}))
        return SimpleNamespace(constraints=dict(constraints or {}))


def _mapping():
    from mewpy.simulation.hybrid import Map, Mapper

    # P is shared by v1 and v3; Q belongs to v2 only
    return Map(
        {
            "v1": Mapper("vmax1", {"P": 3600.0}, {"P": 7200.0}, 1),
            "v2": Mapper("vmax2", {"P": 3600.0, "Q": 1800.0}, {"P": 3600.0, "Q": 1800.0}, 1),
            "v3": Mapper("vmax3", {"P": 3600.0}, {"P": 3600.0}, 1),
        }
    )


PARAMS = {"vmax1": 1.0, "vmax2": 0.0, "vmax3": 0.5}
RATES = {"v1": 0.2, "v2": 0.0, "v3": 0.1}


def _expected(vmax1=1.0, tol=0.05):
    # v2 has vmax 0 and rate 0 (backward kcats); P is drawn once for v1, v2 and v3, so its bounds add up
    p_lb = (1 - tol) * TO_CB * (0.2 / 3600 + 0.0 / 3600 + 0.1 / 3600)
    p_ub = TO_CB * (vmax1 / 3600 + 0.0 / 3600 + 0.5 / 3600)
    return {"draw_prot_P": (p_lb, p_ub), "draw_prot_Q": (0.0, 0.0)}


def _assert_bounds(found, expected):
    assert set(found) >= set(expected)
    for r_id, (lb, ub) in expected.items():
        assert found[r_id] == pytest.approx((lb, ub), abs=1e-12), r_id


def test_kinetic_rates_convert_to_mmol_per_gdw_per_hour():
    """0.2004 mM/s (the Chassagnole PTS rate) is 0.2004 x 3600 / 564 = 1.279 mmol/gDW/h"""
    from mewpy.io.sbml import load_ODEModel
    from mewpy.simulation.hybrid import HybridSimulation

    cobra = pytest.importorskip("cobra")
    kmodel = load_ODEModel(str(DATA.joinpath("chassagnole2002.xml")))
    cbmodel = cobra.io.read_sbml_model(str(DATA.joinpath("e_coli_core.xml")))
    hybrid = HybridSimulation(kmodel, cbmodel, mapping={"vPTS": ("GLCpts", 1)})
    assert hybrid.unit_conv(0.2004) == pytest.approx(0.2004 * TO_CB)
    assert hybrid.mapping_conversion({"vPTS": 0.2004}) == {"GLCpts": pytest.approx(1.2791, abs=1e-4)}
    with pytest.warns(DeprecationWarning):
        HybridSimulation(kmodel, cbmodel, mapping={"vPTS": ("GLCpts", 1)}, D=0.278e-4)


def test_chemostat_glucose_uptake_supports_growth():
    """Chassagnole's steady state is a glucose-limited chemostat at D = 0.1 1/h: its PTS rate, as the glucose
    uptake of E. coli core, gives a growth rate of that order (0.0743 1/h; ATPM maintenance takes the rest)"""
    from mewpy.io.sbml import load_ODEModel
    from mewpy.simulation import get_simulator
    from mewpy.simulation.hybrid import HybridSimulation

    cobra = pytest.importorskip("cobra")
    kmodel = load_ODEModel(str(DATA.joinpath("chassagnole2002.xml")))
    sim = get_simulator(cobra.io.read_sbml_model(str(DATA.joinpath("e_coli_core.xml"))))
    result = HybridSimulation(kmodel, sim, mapping={"vPTS": ("EX_glc__D_e", -1)}).simulate(method="pFBA")
    assert result.fluxes["EX_glc__D_e"] == pytest.approx(-1.279, abs=1e-3)
    assert 0.05 < result.fluxes["Biomass_Ecoli_core"] < 0.12


def test_gecko_simulation_enzyme_bounds(monkeypatch):
    """Usage bounds: at most vmax/kcat, at least (1 - tol) |v|/kcat; a vmax of 0 blocks the enzyme;
    an enzyme shared by several kinetic reactions gets the sums"""
    from mewpy.simulation import hybrid

    monkeypatch.setattr(hybrid, "KineticSimulation", _FakeKineticSimulation)
    _FakeKineticSimulation.fluxes = RATES
    sim = object.__new__(hybrid.HybridGeckoSimulation)
    sim.__dict__.update(
        kmodel=SimpleNamespace(merge_constants=lambda: dict(PARAMS)),
        sim=_RecordingSimulator(),
        t_points=[0, 1e9],
        timeout=0,
        gDW=GDW,
        enzyme_mapping=_mapping(),
        protein_prefix="draw_prot_",
    )
    user = {"EX_glc__D_e": (-10, 0)}
    sim.simulate(constraints=user)
    _assert_bounds(sim.sim.constraints[-1], _expected())
    assert user == {"EX_glc__D_e": (-10, 0)}  # the caller's dict is not changed


def test_gecko_problem_scales_vmax_by_the_candidate_factor(monkeypatch):
    """A candidate's vmax factor (0 deletes the kinetic enzyme) also scales the enzyme's maximal usage"""
    from mewpy.problems import hybrid as phybrid

    monkeypatch.setattr(phybrid.GeckoOUProblem, "solution_to_constraints", lambda self, candidate: {})
    _FakeKineticSimulation.fluxes = dict(RATES, v1=0.0)
    _FakeKineticSimulation.calls = []
    problem = object.__new__(phybrid.HybridGeckoOUProblem)
    problem.__dict__.update(
        vmaxs=["vmax1", "vmax2", "vmax3"],
        prot_prefix="draw_prot_",
        ksim=_FakeKineticSimulation(),
        initcond=None,
        kmodel=SimpleNamespace(merge_constants=lambda: dict(PARAMS)),
        enzyme_mapping=_mapping(),
        gDW=GDW,
        apply_lb=True,
        lb_tolerance=0.05,
    )
    constraints = problem.solution_to_constraints({"vmax1": 0.0})
    assert _FakeKineticSimulation.calls[-1] == {"vmax1": 0.0}
    # v1 deleted: no rate and no capacity from v1, and the zero rate selects v1's backward kcat
    p_lb = 0.95 * TO_CB * (0.1 / 3600)
    p_ub = TO_CB * (0.5 / 3600)
    _assert_bounds(constraints, {"draw_prot_P": (p_lb, p_ub), "draw_prot_Q": (0.0, 0.0)})


def test_hybrid_simulation_projection_uses_the_constraints(monkeypatch):
    """The lMOMA projection of the kinetic fluxes sees the user's constraints (e.g. a knockout)"""
    from mewpy.simulation import hybrid

    seen = {}

    def fake_partial_lmoma(model, reactions, biomass, constraints=None):
        seen["constraints"] = dict(constraints) if constraints else None
        return SimpleNamespace(values={r: 0.5 * v for r, v in reactions.items()})

    monkeypatch.setattr(hybrid, "KineticSimulation", _FakeKineticSimulation)
    monkeypatch.setattr(hybrid, "_partial_lMOMA", fake_partial_lmoma)
    _FakeKineticSimulation.fluxes = {"vPTS": 0.2004}
    sim = object.__new__(hybrid.HybridSimulation)
    sim.__dict__.update(
        kmodel=SimpleNamespace(ratelaws={"vPTS": None}),
        sim=_RecordingSimulator(reactions=["GLCpts", "PGI"], objective={"Biomass": 1}),
        mapping={"vPTS": ("GLCpts", 1)},
        t_points=[0, 1e9],
        timeout=0,
        gDW=GDW,
    )
    user = {"PGI": 0}
    sim.simulate(constraints=user)
    assert seen["constraints"] == {"PGI": 0}
    projected = 0.5 * 0.2004 * TO_CB
    assert sim.sim.constraints[-1] == {"PGI": 0, "GLCpts": pytest.approx((projected, projected))}
    assert user == {"PGI": 0}


def _range_hybrid(monkeypatch, rate):
    from mewpy.simulation import hybrid

    monkeypatch.setattr(hybrid, "KineticSimulation", _FakeKineticSimulation)
    _FakeKineticSimulation.fluxes = {"vPGI": rate}
    sim = object.__new__(hybrid.HybridSimulation)
    sim.__dict__.update(
        kmodel=SimpleNamespace(ratelaws={"vPGI": None}),
        sim=_RecordingSimulator(reactions=["PGI"], objective={"Biomass": 1}),
        mapping={"vPGI": ("PGI", 1)},
        t_points=[0, 1e9],
        timeout=0,
        gDW=GDW,
    )
    return sim


def _bounds(value):
    return tuple(value) if isinstance(value, (tuple, list)) else (value, value)


@pytest.mark.parametrize("user, expected", [(0, (0, 0)), ((0, 0.5), (0.1 * TO_CB - 0.25, 0.5))])
def test_amplitude_ranges_respect_the_constraints(monkeypatch, user, expected):
    """A kinetic range is intersected with the user's bound on the same reaction; a knockout (0) has no overlap
    with it and holds"""
    sim = _range_hybrid(monkeypatch, 0.1)  # 0.638 mmol/gDW/h, so the range is 0.388 to 0.888
    sim.simulate(constraints={"PGI": user}, amplitude=0.5)
    assert _bounds(sim.sim.constraints[-1]["PGI"]) == pytest.approx(expected)


@pytest.mark.parametrize("user, expected", [(0, (0, 0)), ((0, 0.7), (0.1 * TO_CB, 0.7))])
def test_sampled_ranges_respect_the_constraints(monkeypatch, user, expected):
    import pandas as pd

    sim = _range_hybrid(monkeypatch, 0.1)
    samples = pd.DataFrame({"vPGI": [0.1, 0.1, 0.12, 0.14]})  # 10th-90th percentiles: 0.638 to 0.855 mmol/gDW/h
    sim.simulate_distribution(samples, constraints={"PGI": user})
    assert _bounds(sim.sim.constraints[-1]["PGI"]) == pytest.approx(expected)


def test_kinetic_ranges_respect_the_simulators_own_constraints(monkeypatch):
    """A knockout given to the simulator (not to simulate) also holds against a kinetic range"""
    sim = _range_hybrid(monkeypatch, 0.1)
    sim.sim._constraints = {"PGI": 0}
    sim.simulate(amplitude=0.5)
    assert _bounds(sim.sim.constraints[-1]["PGI"]) == (0, 0)
