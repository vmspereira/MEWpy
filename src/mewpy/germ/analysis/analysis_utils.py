from collections import namedtuple
from dataclasses import dataclass
from math import exp, log1p
from typing import TYPE_CHECKING, Callable, Dict, Tuple, Union

from mewpy.germ.algebra import And, Or
from mewpy.solvers.solution import Status
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.germ.models import MetabolicModel, Model, RegulatoryModel
    from mewpy.solvers import Solution


def decode_solver_solution(solution: "Solution", minimize: bool = False) -> Tuple[float, str]:
    """
    It decodes the solution of the solver and returns the objective value.

    :param solution: the solution of the solver
    :param minimize: whether the objective was minimized. An unbounded minimization is decoded as the
    lower bound, an unbounded maximization as the upper bound
    :return: the objective value and the status of the solution
    """
    sol_status = solution.status.value

    if solution.status == Status.OPTIMAL:
        sol_f_obj = solution.fobj

    elif solution.status == Status.UNBOUNDED or solution.status == Status.INF_OR_UNB:

        sol_f_obj = ModelConstants.REACTION_LOWER_BOUND if minimize else ModelConstants.REACTION_UPPER_BOUND

    else:
        sol_f_obj = 0.0

    return sol_f_obj, sol_status


def run_method_and_decode(
    method,
    objective: Union[str, Dict[str, float]] = None,
    constraints: Dict[str, Tuple[float, float]] = None,
    **kwargs,
) -> Tuple[float, str]:
    """
    It runs a method and decodes the objective value and status returned by the solver.
    :param method: the method to be run
    :param objective: an alternative temporary objective function
    :param constraints: alternative temporary constraints
    :param kwargs: additional arguments to be passed to the method
    :return: the objective value and the status of the solution
    """
    solver_kwargs = {"get_values": False}

    if objective:
        if hasattr(objective, "keys"):
            solver_kwargs["linear"] = objective.copy()
        else:
            solver_kwargs["linear"] = {str(objective): 1.0}

    if constraints:
        solver_kwargs["constraints"] = constraints

    if "minimize" in kwargs:
        solver_kwargs["minimize"] = kwargs["minimize"]

    solution = method.optimize(to_solver=True, solver_kwargs=solver_kwargs, **kwargs)
    minimize = bool(kwargs.get("minimize", False))
    objective_value, status = decode_solver_solution(solution=solution, minimize=minimize)
    return objective_value, status


# ---------------------------------
# CoRegFlux utils
# ---------------------------------
CoRegMetabolite = namedtuple("CoRegMetabolite", ("id", "concentration", "exchange"))
CoRegBiomass = namedtuple("CoRegBiomass", ("id", "biomass_yield"))


@dataclass
class CoRegResult:
    objective_value: float = None
    values: Dict[str, float] = None
    metabolites: Dict[str, CoRegMetabolite] = None
    biomass: CoRegBiomass = None
    state: Dict[str, float] = None
    constraints: Dict[str, Tuple[float, float]] = None
    status: Status = Status.OPTIMAL


def build_metabolites(
    model: Union["Model", "MetabolicModel", "RegulatoryModel"], metabolites: Dict[str, float]
) -> Dict[str, CoRegMetabolite]:
    res = {}

    # Build map from metabolite to exchange reaction
    # Handle both APIs: RegulatoryExtension has get_exchange_reactions(), legacy has exchanges dict
    if hasattr(model, "get_exchange_reactions"):
        exchange_reactions = model.get_exchange_reactions()
    elif hasattr(model, "exchanges"):
        exchange_reactions = model.exchanges
    else:
        exchange_reactions = []

    met_to_exchange = {}
    for ex_rxn_id in exchange_reactions:
        # Handle both APIs for getting reactions
        if hasattr(model, "get_reaction"):
            rxn = model.get_reaction(ex_rxn_id)
        else:
            rxn = model.reactions[ex_rxn_id]

        for met in rxn.stoichiometry.keys():
            # legacy reactions key the stoichiometry by Metabolite objects
            met_to_exchange[getattr(met, "id", met)] = ex_rxn_id

    for metabolite, concentration in metabolites.items():
        # Find exchange reaction for this metabolite
        exchange = met_to_exchange.get(metabolite)
        if exchange is None:
            # Skip metabolites without exchange reactions
            continue

        res[metabolite] = CoRegMetabolite(id=metabolite, concentration=concentration, exchange=exchange)
    return res


def build_biomass(model: Union["Model", "MetabolicModel", "RegulatoryModel"], biomass: float) -> CoRegBiomass:
    # model.objective is a dict with reaction IDs as keys (might be Variable objects or strings)
    variable = next(iter(model.objective))

    # Extract just the reaction ID string
    if hasattr(variable, "id"):
        # Variable object - use its id attribute
        variable_id = variable.id
    else:
        # String - might be full representation like "BIOMASS || equation"
        # Extract just the ID part before " || " if present
        variable_id = str(variable).split(" || ")[0].strip()

    return CoRegBiomass(id=variable_id, biomass_yield=biomass)


def concentration_to_lb(concentration, biomass, time_step):
    return concentration / (biomass * time_step)


def _soft_plus(x):
    # log(1 + exp(x)) without overflowing exp for large x (e.g. scaled states of 1000)
    if x > 0:
        return x + log1p(exp(-x))
    return log1p(exp(x))


def lb_soft_plus(soft_plus, coefficient):
    return -1 * _soft_plus(soft_plus + abs(coefficient))


def ub_soft_plus(soft_plus, coefficient):
    return _soft_plus(soft_plus + coefficient)


def euler_step_biomass(old_biomass_yield, growth_rate, time_step):
    return old_biomass_yield * exp(growth_rate * time_step)


def euler_step_metabolites(metabolite_concentration, metabolite_rate, old_biomass_yield, growth_rate, time_step):
    next_concentration = metabolite_rate / growth_rate * old_biomass_yield * (1 - exp(growth_rate * time_step))
    return metabolite_concentration - next_concentration


def biomass_yield_to_rate(biomass):
    return 1 * biomass


def metabolites_constraints(
    constraints: Dict[str, Tuple[float, float]],
    metabolites: Dict[str, CoRegMetabolite],
    biomass: CoRegBiomass,
    time_step: float,
):
    constraints = constraints.copy()
    for metabolite in metabolites.values():

        next_lb = concentration_to_lb(
            concentration=metabolite.concentration, biomass=biomass.biomass_yield, time_step=time_step
        )

        rxn = metabolite.exchange

        if next_lb < ModelConstants.TOLERANCE:
            next_lb = 0

        prev_lb, prev_ub = constraints[rxn]

        if prev_lb == 0:
            next_lb = -next_lb

        else:

            if abs(prev_lb) < abs(next_lb):
                next_lb = prev_lb

            else:
                next_lb = -abs(next_lb)

        constraints[rxn] = (next_lb, prev_ub)

    return constraints


class IdResolver:
    """
    Maps user ids to the model's native ids: an id in the other SBML form (with or without the R_/M_/G_
    prefix, see Simulator.resolve) or an original TRN id renamed by RegulatoryExtension (id_map).
    Ids that are already native, or that match nothing, are returned unchanged.
    """

    def __init__(self, model):
        self.model = model
        simulator = getattr(model, "simulator", None)
        if simulator is None:
            from mewpy.simulation import get_simulator

            simulator = get_simulator(model)
        self.simulator = simulator
        self.id_map = dict(getattr(model, "id_map", None) or {})
        known = set()
        for attribute in ("reactions", "genes", "metabolites", "regulators", "targets"):
            known.update(getattr(model, attribute, None) or ())
        self.known = known

    def id(self, identifier):
        if not isinstance(identifier, str) or identifier in self.known:
            return identifier
        if identifier in self.id_map:
            return self.id_map[identifier]
        for kind in ("reaction", "gene", "metabolite"):
            native = self.simulator.resolve(identifier, kind)
            if native != identifier:
                return native
        return identifier

    def ids(self, identifiers):
        return [self.id(identifier) for identifier in identifiers]

    def keys(self, mapping):
        """A copy of mapping with resolved keys (tuple keys, e.g. PROM's (target, regulator), item by item)."""
        if not mapping:
            return mapping
        return {
            (tuple(self.id(item) for item in key) if isinstance(key, tuple) else self.id(key)): value
            for key, value in mapping.items()
        }

    def solver_kwargs(self, solver_kwargs):
        """A copy of solver_kwargs with the keys of its constraints and linear objective resolved."""
        if not solver_kwargs:
            return solver_kwargs
        resolved = dict(solver_kwargs)
        for key in ("constraints", "linear"):
            if isinstance(resolved.get(key), dict):
                resolved[key] = self.keys(resolved[key])
            elif isinstance(resolved.get(key), str):
                resolved[key] = self.id(resolved[key])
        return resolved


def _metabolite_exchanges(model) -> Dict[str, Tuple[str, float]]:
    """Metabolite id -> (exchange reaction id, lower bound), from the simulator of a RegulatoryExtension."""
    if isinstance(model.reactions, dict) or not hasattr(model, "get_exchange_reactions"):
        # legacy GERM models link metabolite regulators to their exchange reaction directly
        return {}
    exchanges = {}
    for ex_rxn_id in model.get_exchange_reactions():
        rxn = model.get_reaction(ex_rxn_id)
        for met_id in rxn["stoichiometry"]:
            exchanges[met_id] = (ex_rxn_id, rxn["lb"])
    return exchanges


def regulator_exchanges(model: Union["Model", "MetabolicModel", "RegulatoryModel"]) -> Dict[str, Tuple[str, float]]:
    """
    Metabolite regulators and their exchange reaction: regulator id -> (exchange reaction id, lower bound).
    Legacy GERM metabolites know their exchange reaction; RegulatoryExtension regulators are matched to the
    simulator metabolites (the TRN ids are reconciled with the simulator ids when the network is loaded).
    """
    exchanges = _metabolite_exchanges(model)
    result = {}
    for regulator in model.yield_regulators():
        if isinstance(regulator, tuple):
            regulator = regulator[1]
        if regulator.is_metabolite() and getattr(regulator, "exchange_reaction", None):
            result[regulator.id] = (regulator.exchange_reaction.id, regulator.exchange_reaction.lower_bound)
        elif regulator.id in exchanges:
            result[regulator.id] = exchanges[regulator.id]
    return result


def regulators_initial_state(
    model: Union["Model", "MetabolicModel", "RegulatoryModel"], state: Dict[str, float] = None
) -> Dict[str, float]:
    """
    Initial state of every regulator: the value given in state; for a metabolite with an exchange reaction,
    the uptake bound |lower bound| of the exchange (or the value given in state for the exchange), i.e. whether
    the metabolite is in the medium; otherwise the maximum coefficient of the regulator.

    :param model: the model to be simulated
    :param state: values given explicitly for regulators or exchange reactions
    :return: dict of regulator ids and their initial values
    """
    if not state:
        state = {}

    exchanges = regulator_exchanges(model)
    initial_state = {}
    for regulator in model.yield_regulators():
        # RegulatoryExtension yields (id, regulator) tuples, legacy models yield regulators
        if isinstance(regulator, tuple):
            regulator = regulator[1]

        if regulator.id in state:
            initial_state[regulator.id] = state[regulator.id]
            continue

        exchange = exchanges.get(regulator.id)
        if exchange is not None:
            exchange_id, lower_bound = exchange
            initial_state[regulator.id] = state[exchange_id] if exchange_id in state else abs(lower_bound)
        else:
            initial_state[regulator.id] = max(regulator.coefficients)

    return initial_state


def gene_defaults(model: Union["Model", "MetabolicModel", "RegulatoryModel"]) -> Dict[str, float]:
    """
    Default activity of every metabolic gene: the maximum coefficient for legacy GERM genes
    (what a GPR symbol falls back to when the gene is missing from the state), 1.0 otherwise.
    """
    genes = model.genes
    if isinstance(genes, dict):
        return {gene_id: max(getattr(gene, "coefficients", (1.0,))) for gene_id, gene in genes.items()}
    return {getattr(gene, "id", gene): 1.0 for gene in genes}


def reaction_gpr(model: Union["Model", "MetabolicModel", "RegulatoryModel"], rxn_id: str):
    """Parsed GPR of a reaction for both RegulatoryExtension (get_parsed_gpr) and legacy GERM models."""
    if hasattr(model, "get_parsed_gpr"):
        return model.get_parsed_gpr(rxn_id)
    return model.reactions[rxn_id].gpr


def gpr_gene_ids(gpr) -> list:
    """Names of all gene symbols in a GPR (an Expression or a Symbolic), at any nesting depth."""
    symbolic = getattr(gpr, "symbolic", gpr)
    return [atom.name for atom in symbolic.atoms(symbols_only=True)]


def continuous_gpr(
    model: Union["Model", "MetabolicModel", "RegulatoryModel"], state: Dict[str, float], scale: bool = False
):
    operators = {And: min, Or: max}

    states = {}
    # Iterate over reaction IDs and get parsed GPR for each
    for rxn_id in model.reactions:
        gpr = reaction_gpr(model, rxn_id)

        if gpr.is_none:
            continue

        # gene symbols at any depth (gpr.variables only holds the root operator's children)
        if not set(gpr_gene_ids(gpr)).issubset(state):
            continue

        states[rxn_id] = gpr.evaluate(values=state, operators=operators)

    if scale and states:
        _max_state = max(states.values())

        if _max_state > 0:
            return {key: round((val / _max_state), 6) * 1000 for key, val in states.items()}

    return states


def gene_state_constraints(
    model: Union["Model", "MetabolicModel", "RegulatoryModel"],
    constraints: Dict[str, Tuple[float, float]],
    state: Dict[str, float],
    soft_plus: float = 0,
    tolerance: float = ModelConstants.TOLERANCE,
    scale: bool = False,
):
    constraints = constraints.copy()

    reactions_state = continuous_gpr(model=model, state=state, scale=scale)
    # find the reactions which lower bound was changed by the rules
    # find the reactions which upper bound was changed by the rules
    # evaluate the soft plus over the bounds computed using the continuous version of the gpr rules
    # gene_state_bounds does not contain reactions without gpr rules
    for reaction, coefficient in reactions_state.items():
        old_lb, old_ub = constraints[reaction]

        if coefficient > tolerance:

            if abs(old_lb) > tolerance:
                new_lb = lb_soft_plus(soft_plus=soft_plus, coefficient=coefficient)
            else:
                new_lb = old_lb

            if old_ub > tolerance:
                new_ub = ub_soft_plus(soft_plus=soft_plus, coefficient=coefficient)
            else:
                new_ub = old_ub

        else:
            new_lb = old_lb
            new_ub = old_ub

        constraints[reaction] = (new_lb, new_ub)

    return constraints


def system_state_update(
    model: Union["Model", "MetabolicModel", "RegulatoryModel"],
    flux_state: Dict[str, float],
    metabolites: Dict[str, CoRegMetabolite],
    biomass: CoRegBiomass,
    time_step: float,
    biomass_fn: Callable = None,
) -> Tuple[CoRegBiomass, Dict[str, CoRegMetabolite]]:

    growth_rate = flux_state[biomass.id]
    old_biomass_yield = biomass.biomass_yield

    if biomass_fn is not None:
        growth_rate = biomass_fn(growth_rate)

    if growth_rate == 0:
        growth_rate = ModelConstants.TOLERANCE

    biomass_yield = euler_step_biomass(
        old_biomass_yield=old_biomass_yield, growth_rate=growth_rate, time_step=time_step
    )

    next_biomass = build_biomass(model, biomass_yield)

    next_metabolites = {}

    for met_id, met in metabolites.items():

        concentration = euler_step_metabolites(
            metabolite_concentration=met.concentration,
            metabolite_rate=flux_state[met.exchange],
            old_biomass_yield=old_biomass_yield,
            growth_rate=growth_rate,
            time_step=time_step,
        )

        if concentration < 0:
            concentration = 0

        next_metabolites[met_id] = CoRegMetabolite(id=met_id, concentration=concentration, exchange=met.exchange)

    return next_biomass, next_metabolites
