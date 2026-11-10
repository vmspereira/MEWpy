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
Strain design objectives from a batch culture of each candidate, simulated by dynamic FBA: product titer, yield
and volumetric productivity, as defined by DySScO.

Reference:
    Zhuang K, et al. (2013). Dynamic strain scanning optimization: an efficient strain design strategy for balanced
    yield, titer, and productivity. DySScO strategy for strain design. BMC Biotechnol 13:8.
    doi:10.1186/1472-6750-13-8
"""

from typing import Callable, Dict, Mapping, Union

from mewpy.simulation import SimulationMethod, SStatus
from mewpy.simulation.dfba import DFBAResult, DynamicFBA

from .evaluator import PhenotypeEvaluationFunction
from .phenotype import _SimulatorCopyCache

__all__ = ["BatchCulture", "BatchTiter", "BatchYield", "BatchProductivity"]


def _freeze(mapping):
    return frozenset((k, tuple(v) if isinstance(v, list) else v) for k, v in (mapping or {}).items())


class BatchCulture(_SimulatorCopyCache):
    """A batch culture simulated by dynamic FBA (see mewpy.simulation.DynamicFBA).

    The batch objectives (BatchTiter, BatchYield, BatchProductivity) simulate the culture of each candidate, with
    the candidate's modifications as constraints; objectives that share a BatchCulture simulate each candidate once.

    :param initial_biomass: the initial biomass concentration (gDW/L).
    :param initial_concentrations: {exchange id: initial concentration (mM)} of the tracked compounds; the product
        and the substrate of the objectives must be among them.
    :param uptake: {exchange id: maximal uptake rate}, e.g. a MichaelisMenten.
    :param t_end: the batch time (h).
    :param dt: the time step (h).
    :param biomass: the biomass reaction (default: the single reaction of the model's objective).
    :param method: the simulation method of each step (default pFBA).
    :param constraints: constraints applied in every batch, before the candidate's.
    """

    def __init__(
        self,
        initial_biomass: float,
        initial_concentrations: Mapping[str, float],
        uptake: Mapping[str, Callable[[float], float]] = None,
        t_end: float = 24.0,
        dt: float = 0.1,
        biomass: str = None,
        method: Union[str, SimulationMethod] = SimulationMethod.pFBA,
        constraints: Dict = None,
    ):
        if t_end <= 0 or dt <= 0:
            raise ValueError("t_end and dt must be positive")
        self.initial_biomass = initial_biomass
        self.initial_concentrations = dict(initial_concentrations)
        self.uptake = dict(uptake or {})
        self.t_end = t_end
        self.dt = dt
        self.biomass = biomass
        self.method = method
        self.constraints = dict(constraints or {})

    def simulate(self, simulator, constraints: Dict = None) -> DFBAResult:
        """The culture of a simulator, with additional constraints (e.g. a candidate's modifications).

        :param simulator: a Simulator.
        :param constraints: constraints applied after the culture's own.
        :returns: a DFBAResult.
        """
        merged = dict(self.constraints)
        merged.update(constraints or {})
        dfba = DynamicFBA(
            simulator,
            self.initial_biomass,
            self.initial_concentrations,
            uptake=self.uptake,
            biomass=self.biomass,
            method=self.method,
            constraints=merged,
        )
        return dfba.simulate(t_end=self.t_end, dt=self.dt)

    def candidate(self, result) -> DFBAResult:
        """The culture of the candidate behind a simulation result: its simulator, medium and constraints.

        The last candidate's culture is kept, so objectives that share this BatchCulture simulate it once.
        """
        key = (_freeze(result.envcond), _freeze(result.model_constraints), _freeze(result.simulation_constraints))
        last = self.__dict__.get("_last")
        if last is None or last[0] is not result.model or last[1] != key:
            simulator = self._simulator_copy(result)
            last = (result.model, key, self.simulate(simulator, result.simulation_constraints))
            self._last = last
        return last[2]

    def __getstate__(self):
        # caches hold simulators: rebuilt where the culture is unpickled
        state = self.__dict__.copy()
        state.pop("_sim_cache", None)
        state.pop("_last", None)
        return state


class _BatchObjective(PhenotypeEvaluationFunction):
    """An objective computed from the batch culture of each candidate (product concentrations in mM)"""

    def __init__(self, culture: BatchCulture, product: str, maximize: bool = True, worst_fitness: float = None):
        super().__init__(maximize=maximize, worst_fitness=worst_fitness)
        self.culture = culture
        self.product = product

    def _concentrations(self, result, r_id):
        """The concentrations of a tracked compound along the candidate's culture"""
        culture = self.culture.candidate(result)
        native = self._native_id(result, r_id)
        if native not in culture.concentrations:
            raise ValueError(f"{r_id} is not a tracked compound of the batch culture (initial_concentrations)")
        return culture.concentrations[native]

    def get_fitness(self, simul_results, candidate, **kwargs):
        result = simul_results.get(self.culture.method)
        if result is None or result.status not in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL):
            return self.no_solution
        return self._from_culture(result)

    def _from_culture(self, result):
        raise NotImplementedError

    def required_simulations(self):
        # the simulation gives the candidate's simulator and constraints, and screens out infeasible candidates
        return [self.culture.method]


class BatchTiter(_BatchObjective):
    """The product titer: its concentration (mM) at the end of the batch culture.

    :param culture: the BatchCulture.
    :param product: the product's exchange reaction, a tracked compound of the culture.
    """

    def _from_culture(self, result):
        return self._concentrations(result, self.product)[-1]

    def short_str(self):
        return "BatchTiter"

    def method_str(self):
        return f"BatchTiter {self.product} at {self.culture.t_end:g} h"


class BatchYield(_BatchObjective):
    """The product yield of the batch culture: product formed per substrate consumed (mol/mol).

    :param culture: the BatchCulture.
    :param product: the product's exchange reaction, a tracked compound of the culture.
    :param substrate: the substrate's exchange reaction, a tracked compound of the culture.
    """

    def __init__(self, culture: BatchCulture, product: str, substrate: str, maximize: bool = True, worst_fitness=None):
        super().__init__(culture, product, maximize=maximize, worst_fitness=worst_fitness)
        self.substrate = substrate

    def _from_culture(self, result):
        product = self._concentrations(result, self.product)
        substrate = self._concentrations(result, self.substrate)
        consumed = substrate[0] - substrate[-1]
        if consumed <= 0:
            return self.no_solution
        return (product[-1] - product[0]) / consumed

    def short_str(self):
        return "BatchYield"

    def method_str(self):
        return f"BatchYield {self.product} on {self.substrate} at {self.culture.t_end:g} h"


class BatchProductivity(_BatchObjective):
    """The volumetric productivity of the batch culture: product formed per batch time (mM/h), the batch time being
    the culture's t_end, also when growth stops earlier.

    :param culture: the BatchCulture.
    :param product: the product's exchange reaction, a tracked compound of the culture.
    """

    def _from_culture(self, result):
        product = self._concentrations(result, self.product)
        return (product[-1] - product[0]) / self.culture.t_end

    def short_str(self):
        return "BatchProductivity"

    def method_str(self):
        return f"BatchProductivity {self.product} over {self.culture.t_end:g} h"
