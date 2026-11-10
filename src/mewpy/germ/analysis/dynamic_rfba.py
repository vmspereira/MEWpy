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
Dynamic regulatory FBA: a batch culture simulated by dynamic FBA, in which a transcriptional regulatory network
switches reactions off at each time step.

At the start of each step the regulators are set from the culture. A metabolite regulator whose exchange is tracked
is present while its concentration exceeds a threshold; one whose exchange is not tracked is present when the
medium allows its uptake. A reaction regulator takes the flux of the previous step. The other regulators follow their
rules, evaluated until they no longer change. The target genes are then decoded from the regulators, and the
reactions whose GPR is false cannot carry flux in the step.

References:
    Covert MW, Schilling CH, Palsson BO (2001). Regulation of gene expression in flux balance models of metabolism.
    J Theor Biol 213:73-88. doi:10.1006/jtbi.2001.2405
    Covert MW, et al. (2004). Integrating high-throughput and computational data elucidates bacterial networks.
    Nature 429:92-96. doi:10.1038/nature02456
"""

from typing import Callable, Dict, Mapping, Union

from mewpy.germ.analysis.analysis_utils import regulator_exchanges
from mewpy.germ.analysis.rfba import RFBA
from mewpy.simulation import SimulationMethod
from mewpy.simulation.dfba import DFBAResult, DynamicFBA, _bounds
from mewpy.simulation.simulation import current_bounds
from mewpy.util.constants import ModelConstants


class DynamicRFBAResult(DFBAResult):
    """The trajectory of a dynamic regulatory FBA (see DFBAResult).

    :ivar regulatory_states: the regulators' states of each step, {regulator id: state}.
    :ivar inactive: the reactions that the regulation switched off in each step, sorted.
    """

    def __init__(self, result: DFBAResult, regulatory_states, inactive):
        super().__init__(
            result.time,
            result.biomass,
            result.concentrations,
            result.growth_rate,
            result.exchange_fluxes,
            result.status,
        )
        self.regulatory_states = regulatory_states
        self.inactive = inactive


class DynamicRFBA(DynamicFBA):
    """Dynamic regulatory flux balance analysis of a batch culture.

    :param model: a RegulatoryExtension (a metabolic model with a regulatory network). Its simulator's medium and
        constraints apply at every step.
    :param initial_biomass: the initial biomass concentration (gDW/L).
    :param initial_concentrations: {exchange id: initial concentration (mM)} of the tracked medium compounds.
    :param uptake: {exchange id: maximal uptake rate}, e.g. a MichaelisMenten (see DynamicFBA).
    :param biomass: the biomass reaction (default: the single reaction of the objective).
    :param method: the simulation method of each step (default pFBA).
    :param constraints: additional constraints applied at every step.
    :param initial_state: states of regulators, fixed for the whole simulation, and of genes, which replace the
        decoded gene states (e.g. a knockout: {gene: 0}). Ids may be TRN or SBML ids. A reaction regulator not given
        here takes its maximum coefficient in the first step.
    :param threshold: the concentration (mM) above which a tracked metabolite regulator is present.
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
        initial_state: Mapping[str, float] = None,
        threshold: float = 1e-3,
    ):
        super().__init__(
            model.simulator,
            initial_biomass,
            initial_concentrations,
            uptake=uptake,
            biomass=biomass,
            method=method,
            constraints=constraints,
        )
        if threshold < 0:
            raise ValueError("threshold must not be negative")
        self.model = model
        self.rfba = RFBA(model)
        self.fixed_state = self.rfba._ids.keys(dict(initial_state or {}))
        self.threshold = threshold

        sim = self.simulator
        regulators = [r_id for r_id, _ in self.rfba._get_regulators()]
        exchanges = regulator_exchanges(model)
        self._metabolite_regulators = {r: exchanges[r][0] for r in regulators if r in exchanges}
        self._reaction_regulators = [r for r in regulators if r in sim.reactions and r not in exchanges]

    def _regulation(self, states, inactive):
        """The step's regulatory knockouts, from the tracked concentrations and the previous step's fluxes"""
        state = self.rfba.initial_state(self.fixed_state)
        fixed = self.fixed_state

        def regulate(concentrations, fluxes):
            for regulator, r_id in self._metabolite_regulators.items():
                if regulator in fixed:
                    continue
                if r_id in concentrations:
                    c = concentrations[r_id]
                    state[regulator] = c if c > self.threshold else 0.0
                else:
                    lb, _ = _bounds(current_bounds(self.simulator, r_id, self.constraints))
                    state[regulator] = -lb if lb < 0 else 0.0
            if fluxes is not None:
                for regulator in self._reaction_regulators:
                    if regulator not in fixed:
                        v = fluxes.get(regulator, 0.0)
                        state[regulator] = v if abs(v) > ModelConstants.TOLERANCE else 0.0
            # the regulators that follow rules settle, as in RFBA
            state.update(self.rfba.settle_regulatory_state(state, fixed=fixed))
            genes = {**self.rfba.decode_metabolic_state(state), **fixed}
            knockouts = self.rfba.decode_constraints(genes)
            states.append(dict(state))
            inactive.append(sorted(knockouts))
            return knockouts

        return regulate

    def simulate(self, t_end: float, dt: float = 0.1) -> DynamicRFBAResult:
        """Integrates the culture from 0 to t_end (see DynamicFBA.simulate).

        :param t_end: the final time (h).
        :param dt: the time step (h).
        :returns: a DynamicRFBAResult.
        """
        states, inactive = [], []
        self._regulated = self._regulation(states, inactive)
        result = super().simulate(t_end, dt)
        steps = len(result.growth_rate)
        # a step without a feasible solution ends the simulation: its regulation is not part of the trajectory
        return DynamicRFBAResult(result, states[:steps], inactive[:steps])

    def _member(self, regulation=None):
        return super()._member(regulation or self._regulated)
