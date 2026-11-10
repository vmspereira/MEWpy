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
Dynamic flux balance analysis (dFBA), static optimisation approach.

At each time step, uptake kinetics bound the exchange reactions, an FBA (by default pFBA) gives the growth rate and
the exchange fluxes, and biomass and medium concentrations are integrated over the step for that constant growth
rate. The model runs through MEWpy's Simulator, so the same code runs on COBRApy, REFRAMED and GERM models, with
reaction ids given with or without the SBML prefix.

References:
    Varma A, Palsson BO (1994). Stoichiometric flux balance models quantitatively predict growth and metabolic
    by-product secretion in wild-type Escherichia coli W3110. Appl Environ Microbiol 60:3724-3731.
    doi:10.1128/aem.60.10.3724-3731.1994
    Mahadevan R, Edwards JS, Doyle FJ (2002). Dynamic flux balance analysis of diauxic growth in Escherichia coli.
    Biophys J 83:1331-1340. doi:10.1016/S0006-3495(02)73903-9
"""
import logging
import math
from dataclasses import dataclass
from typing import Callable, Dict, List, Mapping, Optional, Union

import pandas as pd

from .simulation import SimulationMethod, Simulator, SStatus, current_bounds

logger = logging.getLogger(__name__)

# growth rates below this value (1/h) are treated as no growth when integrating a step
GROWTH_TOLERANCE = 1e-9


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
            f"DFBAResult(status={self.status!r}, time 0-{self.time[-1]:g} h, steps={len(self.growth_rate)}, "
            f"biomass {self.biomass[0]:g} -> {self.biomass[-1]:g} gDW/L)"
        )

    def plot(self, ids: Optional[List[str]] = None, ax=None):
        """Biomass (left axis) and medium concentrations (right axis) over time.

        :param ids: the exchange ids to plot (default: all tracked).
        :param ax: a matplotlib Axes (default: a new figure).
        :returns: the biomass Axes.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()
        ax.plot(self.time, self.biomass, color="k", label="biomass")
        ax.set_xlabel("time (h)")
        ax.set_ylabel("biomass (gDW/L)")
        right = ax.twinx()
        for r_id in ids or list(self.concentrations):
            right.plot(self.time, self.concentrations[r_id], label=r_id)
        right.set_ylabel("concentration (mM)")
        lines = ax.get_legend_handles_labels()
        more = right.get_legend_handles_labels()
        ax.legend(lines[0] + more[0], lines[1] + more[1], loc="best")
        return ax


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

        if biomass:
            self.biomass = sim.resolve(biomass, "reaction")
        else:
            objective = [r for r, v in sim.objective.items() if v]
            if len(objective) != 1:
                raise ValueError(f"The objective has {len(objective)} reactions; give the biomass reaction (biomass=)")
            self.biomass = objective[0]

    def _step_bounds(self, concentrations, exposure):
        """The bounds of the tracked exchanges for a step: uptake kinetics, and no more uptake than what is left.

        :param exposure: the integral of biomass over the step (gDW h/L).
        """
        bounds = dict(self.constraints)
        for r_id, conc in concentrations.items():
            lb, ub = _bounds(current_bounds(self.simulator, r_id, self.constraints))
            available = conc / exposure
            if r_id in self.uptake:
                max_uptake = min(self.uptake[r_id](conc), available)
            else:
                max_uptake = min(-lb, available) if lb < 0 else 0.0
            lower = -max(max_uptake, 0.0)
            bounds[r_id] = (lower, max(ub, lower))
        return bounds

    def _solve(self, concentrations, exposure):
        return self.simulator.simulate(method=self.method, constraints=self._step_bounds(concentrations, exposure))

    def simulate(self, t_end: float, dt: float = 0.1) -> DFBAResult:
        """Integrates the culture from 0 to t_end.

        The simulation stops early when a step has no feasible solution (e.g. the substrates cannot sustain the
        maintenance demand); the result's status says when and why.

        :param t_end: the final time (h).
        :param dt: the time step (h).
        :returns: a DFBAResult.
        """
        if dt <= 0 or t_end <= 0:
            raise ValueError("t_end and dt must be positive")
        biomass = self.initial_biomass
        conc = dict(self.initial_concentrations)
        time, biomass_t = [0.0], [biomass]
        conc_t = {r: [c] for r, c in conc.items()}
        growth_t, flux_t = [], {r: [] for r in conc}
        status = "completed"
        n_steps = int(round(t_end / dt))
        for step in range(n_steps):
            # the uptake cap assumes the biomass of the start of the step; growth during the step raises the
            # consumption, so a step that would take more than what is left is solved again with the cap set from
            # the biomass integral at that growth rate
            result = self._solve(conc, biomass * dt)
            if _solved(result):
                exposure = _exposure(biomass, result.fluxes[self.biomass], dt)
                if any(c + result.fluxes[r] * exposure < 0 for r, c in conc.items()):
                    result = self._solve(conc, exposure)
            if not _solved(result):
                status = f"stopped at t = {time[-1]:g} h: {result.status.name}"
                logger.info("dFBA %s", status)
                break
            mu = result.fluxes[self.biomass]
            if mu <= GROWTH_TOLERANCE:
                mu = 0.0
            # exact integration over the step for a constant growth rate (Varma & Palsson 1994)
            exposure = _exposure(biomass, mu, dt)
            for r_id in conc:
                v = result.fluxes[r_id]
                conc[r_id] = max(0.0, conc[r_id] + v * exposure)  # max() only absorbs round-off
                conc_t[r_id].append(conc[r_id])
                flux_t[r_id].append(v)
            biomass *= math.exp(mu * dt)
            biomass_t.append(biomass)
            growth_t.append(mu)
            time.append((step + 1) * dt)
        return DFBAResult(time, biomass_t, conc_t, growth_t, flux_t, status)
