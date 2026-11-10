# Copyright (C) 2026 Vitor Pereira

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
Dynamic flux balance analysis (dFBA), static optimisation approach, of a single organism (DynamicFBA) or of a
community whose members share the medium (CommunityDynamicFBA).

At each time step, uptake kinetics bound the exchange reactions, an FBA (by default pFBA) of each organism gives its
growth rate and exchange fluxes, and biomass and medium concentrations are integrated over the step for those
constant growth rates. Uptake is capped at what is left in the medium; when the members together would take more of
a compound than what is left, each member gets a share proportional to its demand. The models run through MEWpy's
Simulator, so the same code runs on COBRApy, REFRAMED and GERM models, the members of a community may use different
backends, and reaction ids can be given with or without the SBML prefix.

References:
    Varma A, Palsson BO (1994). Stoichiometric flux balance models quantitatively predict growth and metabolic
    by-product secretion in wild-type Escherichia coli W3110. Appl Environ Microbiol 60:3724-3731.
    doi:10.1128/aem.60.10.3724-3731.1994
    Mahadevan R, Edwards JS, Doyle FJ (2002). Dynamic flux balance analysis of diauxic growth in Escherichia coli.
    Biophys J 83:1331-1340. doi:10.1016/S0006-3495(02)73903-9
    Zhuang K, et al. (2011). Genome-scale dynamic modeling of the competition between Rhodoferax and Geobacter in
    anoxic subsurface environments. ISME J 5:305-316. doi:10.1038/ismej.2010.117
    Harcombe WR, et al. (2014). Metabolic resource allocation in individual microbes determines ecosystem
    interactions and spatial dynamics. Cell Rep 7:1104-1115. doi:10.1016/j.celrep.2014.03.070
"""

import itertools
import logging
import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Mapping, Optional, Union

import pandas as pd

from .simulation import SimulationMethod, Simulator, SStatus, current_bounds

logger = logging.getLogger(__name__)

# growth rates below this value (1/h) are treated as no growth when integrating a step
GROWTH_TOLERANCE = 1e-9
# a step may take up to this much (mM) more of a compound than what is left: solver round-off
CONCENTRATION_TOLERANCE = 1e-12


@dataclass
class MichaelisMenten:
    """Maximal uptake rate vmax * S / (km + S), in mmol/gDW/h, at a substrate concentration S in mM.

    :param vmax: the maximal uptake rate (mmol/gDW/h).
    :param km: the half-saturation constant (mM).
    """

    vmax: float
    km: float

    def __call__(self, concentration: float) -> float:
        if concentration <= 0:
            return 0.0
        return self.vmax * concentration / (self.km + concentration)


def _bounds(value):
    return tuple(value) if isinstance(value, (tuple, list)) else (value, value)


def _exposure(biomass, mu, dt):
    """The integral of biomass over a step at a constant growth rate mu (gDW h/L): a flux v (mmol/gDW/h) changes a
    concentration by v times this value"""
    if mu > GROWTH_TOLERANCE:
        return biomass * (math.exp(mu * dt) - 1) / mu
    return biomass * dt


def _solved(result):
    return result.status in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL) and result.fluxes is not None


def _biomass_reaction(simulator, biomass):
    if biomass:
        return simulator.resolve(biomass, "reaction")
    objective = [r for r, v in simulator.objective.items() if v]
    if len(objective) != 1:
        raise ValueError(f"The objective has {len(objective)} reactions; give the biomass reaction (biomass=)")
    return objective[0]


class _Member:
    """An organism of a dFBA: its simulator, the native exchange id of each tracked compound, its uptake kinetics
    {compound: callable}, biomass reaction and constraints.

    regulation, if given, is called at the start of each step with the member's tracked concentrations
    {compound: mM} and the fluxes of its previous step (None in the first step, {} after a step without a feasible
    solution); it returns {reaction: bounds} for the step, applied after the uptake caps (e.g. regulatory knockouts).
    """

    def __init__(self, simulator, exchanges, uptake, biomass, constraints, method, regulation=None):
        self.simulator = simulator
        self.exchanges = exchanges
        self.uptake = uptake
        self.biomass = biomass
        self.constraints = constraints
        self.method = method
        self.regulation = regulation
        self.step_constraints = {}
        self.last_fluxes = None

    def _current_bounds(self, compound):
        return _bounds(current_bounds(self.simulator, self.exchanges[compound], self.constraints))

    def max_uptake(self, compound, concentration):
        """The uptake rate (mmol/gDW/h) that the kinetics allow, or else the uptake bound of the model"""
        if compound in self.uptake:
            return max(self.uptake[compound](concentration), 0.0)
        lb, _ = self._current_bounds(compound)
        return -lb if lb < 0 else 0.0

    def solve(self, caps):
        """Simulates the member with the uptake of each tracked exchange capped at caps {compound: mmol/gDW/h}.

        :returns: ((growth rate, {compound: exchange flux}, all fluxes) or None if there is no feasible solution,
            the status)
        """
        bounds = dict(self.constraints)
        for compound, cap in caps.items():
            lb, ub = self._current_bounds(compound)
            # kinetics replace the uptake bound; a compound without kinetics that cannot be taken up keeps its bounds
            lower = -cap if (compound in self.uptake or lb < 0) else lb
            bounds[self.exchanges[compound]] = (lower, max(ub, lower))
        bounds.update(self.step_constraints)
        result = self.simulator.simulate(method=self.method, constraints=bounds)
        if not _solved(result):
            return None, result.status
        mu = result.fluxes[self.biomass]
        fluxes = {compound: result.fluxes[r_id] for compound, r_id in self.exchanges.items()}
        return (mu if mu > GROWTH_TOLERANCE else 0.0, fluxes, result.fluxes), result.status


def _solve_step(members, biomass, concentrations, dt):
    """Solves one step for every member, sharing out the compounds that the members together would overdraw.

    The first caps assume each member alone and its biomass at the start of the step. If the members' consumption
    over the step (growth included) exceeds what is left of a compound, each consumer's cap on it is cut to its
    share, proportional to its demand, and the members whose caps changed are solved again. Caps only decrease.

    :returns: ({member: (growth rate, exchange fluxes, all fluxes) or None}, {member: status})
    """
    caps = {
        m: {
            c: min(member.max_uptake(c, concentrations[c]), concentrations[c] / (biomass[m] * dt))
            for c in member.exchanges
        }
        for m, member in members.items()
    }
    solutions, statuses = {}, {}
    pending = list(members)
    # at most one pass per tracked compound, plus one to check the last caps
    for _ in range(len(concentrations) + 1):
        for m in pending:
            solutions[m], statuses[m] = members[m].solve(caps[m])
        exposure = {m: _exposure(biomass[m], s[0], dt) for m, s in solutions.items() if s}
        changed = set()
        for c, available in concentrations.items():
            demand = {m: -s[1][c] * exposure[m] for m, s in solutions.items() if s and s[1].get(c, 0.0) < 0}
            total = sum(demand.values())
            if total > available + CONCENTRATION_TOLERANCE:
                for m, d in demand.items():
                    caps[m][c] = min(caps[m][c], available * d / total / exposure[m])
                    changed.add(m)
        if not changed:
            return solutions, statuses
        pending = [m for m in members if m in changed]
    logger.warning("dFBA: the members still take more than what is left after sharing; the excess is dropped")
    return solutions, statuses


@dataclass
class _Run:
    time: list
    biomass: dict
    concentrations: dict
    growth_rate: dict
    exchange_fluxes: dict
    stopped: Optional[dict]  # the members' statuses at the step where none had a feasible solution


def _integrate(members, initial_biomass, initial_concentrations, t_end, dt):
    """Integrates the members on a shared medium from 0 to t_end. A member without a feasible solution in a step does
    not grow or exchange in that step; the integration stops at the first step where no member has one."""
    if dt <= 0 or t_end <= 0:
        raise ValueError("t_end and dt must be positive")
    biomass = dict(initial_biomass)
    conc = dict(initial_concentrations)
    run = _Run(
        time=[0.0],
        biomass={m: [x] for m, x in biomass.items()},
        concentrations={c: [v] for c, v in conc.items()},
        growth_rate={m: [] for m in members},
        exchange_fluxes={m: {c: [] for c in member.exchanges} for m, member in members.items()},
        stopped=None,
    )
    for step in range(int(round(t_end / dt))):
        for member in members.values():
            if member.regulation is not None:
                member_conc = {c: conc[c] for c in member.exchanges}
                member.step_constraints = member.regulation(member_conc, member.last_fluxes)
        solutions, statuses = _solve_step(members, biomass, conc, dt)
        if all(s is None for s in solutions.values()):
            run.stopped = statuses
            break
        change = dict.fromkeys(conc, 0.0)
        for m, member in members.items():
            mu, fluxes, member.last_fluxes = solutions[m] or (0.0, dict.fromkeys(member.exchanges, 0.0), {})
            # exact integration over the step for a constant growth rate (Varma & Palsson 1994)
            exposure = _exposure(biomass[m], mu, dt)
            for c, v in fluxes.items():
                change[c] += v * exposure
                run.exchange_fluxes[m][c].append(v)
            biomass[m] *= math.exp(mu * dt)
            run.biomass[m].append(biomass[m])
            run.growth_rate[m].append(mu)
        for c in conc:
            conc[c] = max(0.0, conc[c] + change[c])  # round-off, or an excess left after sharing (logged)
            run.concentrations[c].append(conc[c])
        run.time.append((step + 1) * dt)
    return run


def _plot(time, biomass, concentrations, ids, ax):
    """Biomass {label: series} on the left axis (solid lines) and concentrations on the right axis (dashed lines)"""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    colors = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])
    for label, series in biomass.items():
        ax.plot(time, series, color=next(colors), label=label)
    ax.set_xlabel("time (h)")
    ax.set_ylabel("biomass (gDW/L)")
    right = ax.twinx()
    for c in ids or list(concentrations):
        right.plot(time, concentrations[c], color=next(colors), linestyle="--", label=c)
    right.set_ylabel("concentration (mM)")
    lines = ax.get_legend_handles_labels()
    more = right.get_legend_handles_labels()
    ax.legend(lines[0] + more[0], lines[1] + more[1], loc="best")
    return ax


class DFBAResult:
    """The trajectory of a dynamic FBA.

    :ivar time: the time points (h), from 0.
    :ivar biomass: the biomass concentration (gDW/L) at each time point.
    :ivar concentrations: {exchange id: concentrations (mM) at each time point}.
    :ivar growth_rate: the growth rate (1/h) of each step; the step from time[i] to time[i + 1] has growth_rate[i].
    :ivar exchange_fluxes: {exchange id: flux (mmol/gDW/h) of each step}.
    :ivar status: "completed", or why the simulation stopped.
    """

    def __init__(self, time, biomass, concentrations, growth_rate, exchange_fluxes, status):
        self.time = time
        self.biomass = biomass
        self.concentrations = concentrations
        self.growth_rate = growth_rate
        self.exchange_fluxes = exchange_fluxes
        self.status = status

    @property
    def dataframe(self) -> pd.DataFrame:
        """Time, biomass and concentrations, one row per time point"""
        data = {"time": self.time, "biomass": self.biomass}
        data.update(self.concentrations)
        return pd.DataFrame(data)

    def __repr__(self):
        return (
            f"{type(self).__name__}(status={self.status!r}, time 0-{self.time[-1]:g} h, steps={len(self.growth_rate)}, "
            f"biomass {self.biomass[0]:g} -> {self.biomass[-1]:g} gDW/L)"
        )

    def plot(self, ids: Optional[List[str]] = None, ax=None):
        """Biomass (left axis) and medium concentrations (right axis, dashed) over time.

        :param ids: the exchange ids to plot (default: all tracked).
        :param ax: a matplotlib Axes (default: a new figure).
        :returns: the biomass Axes.
        """
        return _plot(self.time, {"biomass": self.biomass}, self.concentrations, ids, ax)


class CommunityDFBAResult:
    """The trajectory of a community dynamic FBA.

    :ivar time: the time points (h), from 0.
    :ivar biomass: {member: biomass concentration (gDW/L) at each time point}.
    :ivar concentrations: {exchange id without the SBML prefix: concentrations (mM) at each time point}.
    :ivar growth_rate: {member: growth rate (1/h) of each step}; the step from time[i] to time[i + 1] has
        growth_rate[member][i], 0 when the member had no feasible solution in that step.
    :ivar exchange_fluxes: {member: {exchange id: flux (mmol/gDW/h) of each step}}, for the tracked exchanges the
        member has.
    :ivar status: "completed", or why the simulation stopped.
    """

    def __init__(self, time, biomass, concentrations, growth_rate, exchange_fluxes, status):
        self.time = time
        self.biomass = biomass
        self.concentrations = concentrations
        self.growth_rate = growth_rate
        self.exchange_fluxes = exchange_fluxes
        self.status = status

    @property
    def dataframe(self) -> pd.DataFrame:
        """Time, the biomass of each member (one column per member id) and concentrations, one row per time point"""
        data = {"time": self.time}
        data.update(self.biomass)
        data.update(self.concentrations)
        return pd.DataFrame(data)

    def __repr__(self):
        biomass = ", ".join(f"{m} {x[0]:g} -> {x[-1]:g}" for m, x in self.biomass.items())
        return (
            f"CommunityDFBAResult(status={self.status!r}, time 0-{self.time[-1]:g} h, steps={len(self.time) - 1}, "
            f"biomass {biomass} gDW/L)"
        )

    def plot(self, ids: Optional[List[str]] = None, ax=None):
        """The biomass of each member (left axis) and medium concentrations (right axis, dashed) over time.

        :param ids: the exchange ids to plot (default: all tracked).
        :param ax: a matplotlib Axes (default: a new figure).
        :returns: the biomass Axes.
        """
        return _plot(self.time, self.biomass, self.concentrations, ids, ax)


class DynamicFBA:
    """Dynamic flux balance analysis of a batch culture (static optimisation approach).

    :param model: a COBRApy, REFRAMED or GERM model, or a Simulator. The simulator's medium and constraints apply
        at every step; the bounds of the tracked exchanges are set from the medium concentrations.
    :param initial_biomass: the initial biomass concentration (gDW/L).
    :param initial_concentrations: {exchange id: initial concentration (mM)} of the tracked medium compounds,
        substrates and products (e.g. 0 for a secreted product).
    :param uptake: {exchange id: maximal uptake rate}: a MichaelisMenten, or any callable that takes the
        compound's concentration (mM) and returns its maximal uptake rate (mmol/gDW/h). A tracked compound without
        kinetics keeps the uptake bound of the model. In both cases uptake is limited to what is left in the medium.
    :param biomass: the biomass reaction (default: the single reaction of the simulator's objective).
    :param method: the simulation method of each step (default pFBA, which also fixes the exchange fluxes among
        alternative optima).
    :param constraints: additional constraints applied at every step.
    """

    _ORGANISM = "organism"

    def __init__(
        self,
        model,
        initial_biomass: float,
        initial_concentrations: Mapping[str, float],
        uptake: Mapping[str, Callable[[float], float]] = None,
        biomass: str = None,
        method: Union[str, SimulationMethod] = SimulationMethod.pFBA,
        constraints: Dict = None,
    ):
        from . import get_simulator

        self.simulator = model if isinstance(model, Simulator) else get_simulator(model)
        sim = self.simulator
        if initial_biomass <= 0:
            raise ValueError("initial_biomass must be positive")
        self.initial_biomass = initial_biomass
        self.initial_concentrations = sim.resolve_keys(dict(initial_concentrations), "reaction")
        self.uptake = sim.resolve_keys(dict(uptake or {}), "reaction")
        self.constraints = sim.resolve_keys(dict(constraints or {}), "reaction")
        self.method = method

        exchanges = set(sim.get_exchange_reactions())
        not_exchanges = sorted(r for r in list(self.initial_concentrations) + list(self.uptake) if r not in exchanges)
        if not_exchanges:
            raise ValueError(f"Not exchange reactions: {not_exchanges}")
        untracked = sorted(r for r in self.uptake if r not in self.initial_concentrations)
        if untracked:
            raise ValueError(f"Uptake kinetics without an initial concentration: {untracked}")
        negative = sorted(r for r, c in self.initial_concentrations.items() if c < 0)
        if negative:
            raise ValueError(f"Negative initial concentrations: {negative}")
        self.biomass = _biomass_reaction(sim, biomass)

    def simulate(self, t_end: float, dt: float = 0.1) -> DFBAResult:
        """Integrates the culture from 0 to t_end.

        The simulation stops early when a step has no feasible solution (e.g. the substrates cannot sustain the
        maintenance demand); the result's status says when and why.

        :param t_end: the final time (h).
        :param dt: the time step (h).
        :returns: a DFBAResult.
        """
        one = self._ORGANISM
        run = _integrate({one: self._member()}, {one: self.initial_biomass}, self.initial_concentrations, t_end, dt)
        status = "completed"
        if run.stopped:
            status = f"stopped at t = {run.time[-1]:g} h: {run.stopped[one].name}"
            logger.info("dFBA %s", status)
        return DFBAResult(
            run.time, run.biomass[one], run.concentrations, run.growth_rate[one], run.exchange_fluxes[one], status
        )

    def _member(self, regulation=None):
        exchanges = {r_id: r_id for r_id in self.initial_concentrations}
        return _Member(self.simulator, exchanges, self.uptake, self.biomass, self.constraints, self.method, regulation)


class CommunityDynamicFBA:
    """Dynamic flux balance analysis of a community in batch culture: each member has its own biomass, model and
    uptake kinetics, and all members share the medium (static optimisation approach).

    In each step every member is simulated on its own, with its uptake capped by its kinetics and by what is left in
    the medium; when the members together would take more of a compound than what is left, each one gets a share
    proportional to its demand. A member whose step has no feasible solution (e.g. the medium cannot cover its
    maintenance) neither grows nor exchanges in that step. Members interact only through the medium: what one member
    secretes is available to the others from the next step on.

    :param community: a CommunityModel, or {member id: a COBRApy, REFRAMED or GERM model, or a Simulator}. Members
        may use different backends. Each member's medium and constraints apply at every step.
    :param initial_biomass: {member id: initial biomass concentration (gDW/L)}, for every member.
    :param initial_concentrations: {exchange id: initial concentration (mM)} of the tracked medium compounds. Ids
        may be given with or without the SBML prefix; a compound is exchanged by the members whose models have that
        exchange reaction.
    :param uptake: {member id: {exchange id: maximal uptake rate}}: a MichaelisMenten, or any callable that takes
        the compound's concentration (mM) and returns its maximal uptake rate (mmol/gDW/h). A tracked compound
        without kinetics keeps the member's uptake bound.
    :param biomass: {member id: biomass reaction} (default: the single reaction of each member's objective).
    :param method: the simulation method of each member and step (default pFBA).
    :param constraints: {member id: constraints applied at every step}.
    """

    def __init__(
        self,
        community,
        initial_biomass: Mapping[str, float],
        initial_concentrations: Mapping[str, float],
        uptake: Mapping[str, Mapping[str, Callable[[float], float]]] = None,
        biomass: Mapping[str, str] = None,
        method: Union[str, SimulationMethod] = SimulationMethod.pFBA,
        constraints: Mapping[str, Dict] = None,
    ):
        from . import get_simulator

        members = getattr(community, "organisms", community)
        self.simulators = {m: s if isinstance(s, Simulator) else get_simulator(s) for m, s in members.items()}
        if not self.simulators:
            raise ValueError("The community has no members")
        uptake, biomass, constraints = uptake or {}, biomass or {}, constraints or {}
        for name, mapping in (
            ("initial_biomass", initial_biomass),
            ("uptake", uptake),
            ("biomass", biomass),
            ("constraints", constraints),
        ):
            unknown = sorted(set(mapping) - set(self.simulators))
            if unknown:
                raise ValueError(f"Unknown members in {name}: {unknown}")
        missing = sorted(set(self.simulators) - set(initial_biomass))
        if missing:
            raise ValueError(f"initial_biomass has no value for the members {missing}")
        not_positive = sorted(m for m, x in initial_biomass.items() if x <= 0)
        if not_positive:
            raise ValueError(f"initial_biomass must be positive: {not_positive}")
        self.initial_biomass = {m: initial_biomass[m] for m in self.simulators}

        # medium compounds by exchange id without the SBML prefix; the native id of each in every member that has it
        self.exchanges = {m: {} for m in self.simulators}
        self.initial_concentrations = {}
        not_exchanges = []
        for r_id, value in initial_concentrations.items():
            compound = None
            for m, sim in self.simulators.items():
                native = sim.resolve(r_id, "reaction")
                if native in set(sim.get_exchange_reactions()):
                    compound = compound or sim.to_canonical(native, "reaction")
                    self.exchanges[m][compound] = native
            if compound is None:
                not_exchanges.append(r_id)
            else:
                self.initial_concentrations[compound] = value
        if not_exchanges:
            raise ValueError(f"Not an exchange reaction of any member: {sorted(not_exchanges)}")
        negative = sorted(c for c, v in self.initial_concentrations.items() if v < 0)
        if negative:
            raise ValueError(f"Negative initial concentrations: {negative}")

        self.uptake = {}
        for m, sim in self.simulators.items():
            kinetics = {
                sim.to_canonical(sim.resolve(r, "reaction"), "reaction"): f for r, f in uptake.get(m, {}).items()
            }
            untracked = sorted(c for c in kinetics if c not in self.exchanges[m])
            if untracked:
                raise ValueError(
                    f"Uptake kinetics of {m} for exchanges it does not have or that are not tracked: {untracked}"
                )
            self.uptake[m] = kinetics
        self.biomass = {m: _biomass_reaction(sim, biomass.get(m)) for m, sim in self.simulators.items()}
        self.constraints = {
            m: sim.resolve_keys(dict(constraints.get(m, {})), "reaction") for m, sim in self.simulators.items()
        }
        self.method = method

    def simulate(self, t_end: float, dt: float = 0.1) -> CommunityDFBAResult:
        """Integrates the community from 0 to t_end.

        The simulation stops early when no member has a feasible solution in a step; the result's status says when.

        :param t_end: the final time (h).
        :param dt: the time step (h).
        :returns: a CommunityDFBAResult.
        """
        members = {
            m: _Member(sim, self.exchanges[m], self.uptake[m], self.biomass[m], self.constraints[m], self.method)
            for m, sim in self.simulators.items()
        }
        run = _integrate(members, self.initial_biomass, self.initial_concentrations, t_end, dt)
        status = "completed"
        if run.stopped:
            statuses = ", ".join(f"{m} {s.name}" for m, s in run.stopped.items())
            status = f"stopped at t = {run.time[-1]:g} h: no member has a feasible solution ({statuses})"
            logger.info("dFBA %s", status)
        return CommunityDFBAResult(
            run.time, run.biomass, run.concentrations, run.growth_rate, run.exchange_fluxes, status
        )
