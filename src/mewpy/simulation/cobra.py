# Copyright (C) 2019-2026 Vitor Pereira and Centre of Biological Engineering,
#     University of Minho, Portugal

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
##############################################################################
Simulation for COBRApy models. Provide a wrapper with a common interface
to COBRApy.

Author: Vítor Pereira
##############################################################################
"""

import logging
from collections import OrderedDict
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from cobra.core.model import Model
from cobra.core.solution import Solution
from cobra.exceptions import OptimizationError, SolverNotFound
from cobra.flux_analysis import moma, pfba, room
from tqdm import tqdm

from mewpy.util.constants import ModelConstants
from mewpy.util.utilities import AttrDict

from . import SimulationMethod, SStatus, get_default_solver
from .simulation import ModelContainer, SimulationResult, Simulator


def _is_partial(model, reference):
    """Whether a reference flux distribution leaves out some of the model's reactions"""
    return any(rxn.id not in reference for rxn in model.reactions)


def _add_partial_moma(model, reference, linear):
    """COBRApy's add_moma over the reactions of a partial reference only: the distance to the reference is minimised
    over those reactions, as REFRAMED's lMOMA does"""
    import cobra.util.solver as sutil
    from optlang.symbolics import Zero, add

    if not linear and sutil.interface_to_str(model.problem) not in sutil.qp_solvers:
        model.solver = sutil.choose_solver(model, qp=True)
    prob = model.problem
    v = prob.Variable("moma_old_objective")
    c = prob.Constraint(model.solver.objective.expression - v, lb=0.0, ub=0.0, name="moma_old_objective_constraint")
    to_add = [v, c]
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)
    obj_vars = []
    for r_id, flux in reference.items():
        r = model.reactions.get_by_id(r_id)
        if linear:
            components = sutil.add_absolute_expression(
                model, r.flux_expression, name="moma_dist_" + r.id, difference=flux, add=False
            )
            to_add.extend(components)
            obj_vars.append(components.variable)
        else:
            dist = prob.Variable("moma_dist_" + r.id)
            const = prob.Constraint(r.flux_expression - dist, lb=flux, ub=flux, name="moma_constraint_" + r.id)
            to_add.extend([dist, const])
            obj_vars.append(dist**2)
    model.add_cons_vars(to_add)
    if linear:
        model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})
    else:
        model.objective = prob.Objective(add(obj_vars), direction="min", sloppy=True)


def _add_partial_room(model, reference, delta=0.03, epsilon=1e-03):
    """COBRApy's add_room over the reactions of a partial reference only: only those reactions count as changed or
    not, as in REFRAMED's ROOM"""
    from optlang.symbolics import Zero

    prob = model.problem
    variable = prob.Variable("room_old_objective")
    constraint = prob.Constraint(
        model.solver.objective.expression - variable, ub=0.0, lb=0.0, name="room_old_objective_constraint"
    )
    model.objective = prob.Objective(Zero, direction="min", sloppy=True)
    vars_and_cons = [variable, constraint]
    obj_vars = []
    for r_id, flux in reference.items():
        rxn = model.reactions.get_by_id(r_id)
        y = prob.Variable("y_" + rxn.id, type="binary")
        w_u = flux + (delta * abs(flux)) + epsilon
        upper_const = prob.Constraint(
            rxn.flux_expression - y * (rxn.upper_bound - w_u), ub=w_u, name="room_constraint_upper_" + rxn.id
        )
        w_l = flux - (delta * abs(flux)) - epsilon
        lower_const = prob.Constraint(
            rxn.flux_expression - y * (rxn.lower_bound - w_l), lb=w_l, name="room_constraint_lower_" + rxn.id
        )
        vars_and_cons.extend([y, upper_const, lower_const])
        obj_vars.append(y)
    model.add_cons_vars(vars_and_cons)
    model.objective.set_linear_coefficients({v: 1.0 for v in obj_vars})


if TYPE_CHECKING:
    from pandas import DataFrame

LOGGER = logging.getLogger(__name__)


class CobraModelContainer(ModelContainer):
    """A basic container for COBRApy models.

    :param model: A metabolic model.

    """

    def __init__(self, model: Model = None):
        self.model = model
        self._gene_to_reaction = None

    @property
    def id(self):
        """The model identifier."""
        return self.model.id

    @id.setter
    def id(self, sid: str):
        self.model.id = sid

    @property
    def reactions(self) -> List[str]:
        """The list of reaction identifiers"""
        return [rxn.id for rxn in self.model.reactions]

    def get_reaction(self, r_id: str) -> Dict[str, dict]:
        """Returns a dictionary of a reaction properties.

        :param r_id: The reaction identifier.
        :type r_id: str
        :return: A dictionary of properties
        :rtype: AttrDict
        """
        r_id = self.resolve(r_id, "reaction")
        rxn = self.model.reactions.get_by_id(r_id)
        stoichiometry = {met.id: val for met, val in rxn.metabolites.items()}
        res = {
            "id": r_id,
            "name": rxn.name,
            "lb": rxn.lower_bound,
            "ub": rxn.upper_bound,
            "stoichiometry": stoichiometry,
        }
        res["gpr"] = rxn.gene_reaction_rule
        res["annotations"] = rxn.annotation
        return AttrDict(res)

    @property
    def genes(self) -> List[str]:
        """Lists the model gene identifiers

        :return: A list of gene identifiers
        :rtype: list
        """
        return [gene.id for gene in self.model.genes]

    def get_gene(self, g_id):
        g_id = self.resolve(g_id, "gene")
        g = self.model.genes.get_by_id(g_id)
        gr = self.get_gene_reactions()
        r = gr[g_id]
        res = {"id": g_id, "name": g.name, "reactions": r}
        return AttrDict(res)

    @property
    def metabolites(self) -> List[str]:
        """Lists the model metabolite identifiers

        :return: A list of metabolite identifiers
        :rtype: list
        """
        return [met.id for met in self.model.metabolites]

    def get_metabolite(self, m_id) -> Dict[str, Any]:
        """Returns a metabolite propeties.

        :param m_id: The metabolite identifier
        :type m_id: str
        :return: A dictionary of properties
        :rtype: Dict[str,Any]
        """
        m_id = self.resolve(m_id, "metabolite")
        met = self.model.metabolites.get_by_id(m_id)
        res = {"id": m_id, "name": met.name, "compartment": met.compartment, "formula": met.formula}
        return AttrDict(res)

    @property
    def medium(self) -> Dict[str, float]:
        """Returns the model default medium.

        :return: The constraints on the model exchanges.
        :rtype: dict
        """
        return self.model.medium

    @property
    def compartments(self) -> dict:
        """Get the dictionary of current compartment descriptions.

        :return: the dictionary of current compartment
        :rtype: dict
        """
        return self.model.compartments

    def get_compartment(self, c_id: str) -> Dict[str, Any]:
        """Get a dictionary of a compartment descriptions.

        :param c_id: The compartment identifier
        :type c_id: str
        :return: a dictionary of a compartment descriptions.
        :rtype: dict
        """
        from cobra.medium import find_external_compartment

        c = self.model.compartments[c_id]
        e = find_external_compartment(self.model)
        res = {"id": c_id, "name": c, "external": (e == c_id)}
        return AttrDict(res)

    def get_gpr(self, reaction_id: str) -> str:
        """Returns the gpr rule (str) for a given reaction ID.

        :param str reaction_id: The reaction identifier.
        :returns: A string representation of the GPR rule.

        """
        reaction_id = self.resolve(reaction_id, "reaction")
        if reaction_id not in self._native_ids("reaction"):
            raise ValueError(f"Reactions {reaction_id} does not exist")
        reaction = self.model.reactions.get_by_id(reaction_id)
        if reaction.gene_reaction_rule:
            return str(reaction.gene_reaction_rule)
        else:
            return None

    def get_exchange_reactions(self) -> List[str]:
        """Get the list of exchange reactions

        :return: The list of exchange reactions
        :rtype: List[str]
        """
        rxns = [r.id for r in self.model.exchanges]
        return rxns

    def get_gene_reactions(self) -> Dict[str, List[str]]:
        """Get a map of genes to reactions.

        :return: A map of genes to reactions.
        :rtype: Dict[str,List[str]]
        """
        if self._gene_to_reaction is None:
            gr = dict()
            for rxn_id in self.reactions:
                rxn = self.model.reactions.get_by_id(rxn_id)
                genes = rxn.genes
                for g in genes:
                    if g.id in gr.keys():
                        gr[g.id].append(rxn_id)
                    else:
                        gr[g.id] = [rxn_id]
            self._gene_to_reaction = gr
        return self._gene_to_reaction


class Simulation(CobraModelContainer, Simulator):

    def __init__(
        self,
        model: Model,
        envcond: Union[dict, None] = None,
        constraints: Union[dict, None] = None,
        solver=None,
        reference: Union[dict, None] = None,
        reset_solver: bool = ModelConstants.RESET_SOLVER,
    ):
        """
        Generic Simulation class for cobra Model.
        Defines the simulation conditions, and makes available a set of methods.

        :param model: A constraint-based model instance.
        :type model: Model
        :param envcond: Dictionary of environmental conditions, defaults to None
        :type envcond: Union[dict, None], optional
        :param constraints: A dictionary of reaction constraints., defaults to None
        :type constraints: Union[dict, None], optional
        :param solver: An instance of the LP solver., defaults to None
        :type solver: _type_, optional
        :param reference: A dictionary of the wild type flux values, defaults to None
        :type reference: Union[dict, None], optional
        :param reset_solver: If the solver is reset prior to each simulation, defaults to ModelConstants.RESET_SOLVER
        :type reset_solver: bool, optional
        :raises ValueError: Model is incompatible.
        """

        if not isinstance(model, Model):
            raise ValueError("Model is incompatible or inexistent")

        self.model = model
        solver_name = get_default_solver()
        try:
            self.model.solver = solver_name
        except SolverNotFound:
            # mewpy names such as 'scip' or 'optlang' have no optlang interface
            interface = self.model.solver.interface.__name__
            LOGGER.debug("No optlang interface for %s; keeping %s", solver_name, interface)
        self._environmental_conditions = OrderedDict() if envcond is None else self._resolve_mapping(envcond)
        constraints = self._resolve_mapping(constraints)
        self._constraints = (
            dict()
            if constraints is None
            else {k: v for k, v in constraints.items() if k not in list(self._environmental_conditions.keys())}
        )

        self.solver = solver
        self._gene_to_reaction = None
        self._reference = self._resolve_mapping(reference)
        self.__status_mapping = {
            "optimal": SStatus.OPTIMAL,
            "unbounded": SStatus.UNBOUNDED,
            "infeasible": SStatus.INFEASIBLE,
            "infeasible_or_unbounded": SStatus.INF_OR_UNB,
            "suboptimal": SStatus.SUBOPTIMAL,
            "unknown": SStatus.UNKNOWN,
        }
        self._reset_solver = reset_solver
        self.reverse_syntax = []
        self._m_r_lookup = None

        self._MAX_STR = "maximize"
        self._MIN_STR = "minimize"

        # apply the env. cond. and additional constraints to the model
        for r_id, bounds in self._environmental_conditions.items():
            self._set_model_reaction_bounds(r_id, bounds)
        for r_id, bounds in self._constraints.items():
            self._set_model_reaction_bounds(r_id, bounds)

        # if modifications on the envirenment are permited
        # during simulations
        self._allow_env_changes = False
        self.biomass_reaction = None
        try:
            self.biomass_reaction = list(self.objective.keys())[0]
        except (IndexError, KeyError, AttributeError):
            # Objective may be empty or not properly defined
            pass

    @property
    def environmental_conditions(self):
        return self._environmental_conditions.copy()

    @environmental_conditions.setter
    def environmental_conditions(self, a):
        raise ValueError("Can not change environmental conditions. Use set_reaction_bounds instead")

    def _set_model_reaction_bounds(self, r_id, bounds):
        if isinstance(bounds, tuple):
            lb = bounds[0]
            ub = bounds[1]
        elif isinstance(bounds, (int, float)):
            lb = bounds
            ub = bounds
        else:
            raise ValueError(f"Invalid bounds definition {bounds}")
        rxn = self.model.reactions.get_by_id(r_id)
        rxn.bounds = (lb, ub)

    @property
    def objective(self):
        from cobra.util.solver import linear_reaction_coefficients

        d = dict(linear_reaction_coefficients(self.model))
        return {k.id: v for k, v in d.items()}

    @objective.setter
    def objective(self, objective):
        objective = self._resolve_objective(objective)
        if isinstance(objective, str):
            self.model.objective = objective
        elif isinstance(objective, dict):
            from cobra.util.solver import set_objective

            linear_coef = {self.model.reactions.get_by_id(r_id): v for r_id, v in objective.items()}
            set_objective(self.model, linear_coef)
        else:
            raise ValueError(
                "The objective must be a reaction identifier or a dictionary of \
                reaction identifier with respective coeficients."
            )

    def add_compartment(self, comp_id, name=None, external=False):
        """Adds a compartment

        :param str comp_id: Compartment ID
        :param str name: Compartment name, default None
        :param bool external: If the compartment is external, default False.
        """
        self.model.compartments = {comp_id: name}

    def add_metabolite(self, id, formula=None, name=None, compartment=None):
        from cobra import Metabolite

        meta = Metabolite(id, formula=formula, name=name, compartment=compartment)
        self.model.add_metabolites([meta])
        self._invalidate_id_cache()

    def add_gene(self, id, name):
        pass

    def add_reaction(
        self,
        rxn_id,
        name=None,
        stoichiometry=None,
        reversible=True,
        lb=ModelConstants.REACTION_LOWER_BOUND,
        ub=ModelConstants.REACTION_UPPER_BOUND,
        gpr=None,
        objective=0,
        replace=True,
        annotations={},
        reaction_type=None,
    ):
        """Adds a reaction to the model

        :param rxn_id: The reaction identifier
        :type rxn_id: str
        :param name: The name of the reaction, defaults to None
        :type name: str, optional
        :param stoichiometry: The stoichiometry of the reaction, defaults to None
        :type stoichiometry: dict[str,float], optional
        :param reversible: If the reaction is reversible or not, defaults to True
        :type reversible: bool, optional
        :param lb: The reaction lower bound, defaults to ModelConstants.REACTION_LOWER_BOUND
        :type lb: float, optional
        :param ub: The reaction upper bound, defaults to ModelConstants.REACTION_UPPER_BOUND
        :type ub: float, optional
        :param gpr: The gene-protein-reactio rule, defaults to None
        :type gpr: str, optional
        :param objective: the objective coeficient, defaults to 0
        :type objective: float, optional
        :param replace: If replaces the reaction with a same identifier, defaults to True
        :type replace: bool, optional
        """
        from cobra.core.reaction import Reaction, set_objective

        reaction = Reaction(rxn_id)

        reaction.name = name if name else rxn_id
        stoich = {self.model.metabolites.get_by_id(self.resolve(k, "metabolite")): v for k, v in stoichiometry.items()}
        reaction.add_metabolites(stoich)

        # in cobrapy reversibility is ignored for compatibility reasons
        # instead the lb is evaluated to infer reversability
        reaction.lower_bound = lb
        reaction.upper_bound = ub
        if gpr and isinstance(gpr, str):
            reaction.gene_reaction_rule = gpr
        if annotations:
            reaction.annotation = annotations
        if replace and rxn_id in self.reactions:
            self.remove_reaction(rxn_id)

        self.model.add_reactions([reaction])
        self._invalidate_id_cache()
        if objective != 0:
            set_objective(self.model, {reaction: objective})

    def remove_reaction(self, r_id):
        """Removes a reaction from the model.

        Args:
            r_id (str): The reaction identifier.
        """
        self.model.remove_reactions([self.resolve(r_id, "reaction")])
        self._invalidate_id_cache()

    def remove_reactions(self, rxn_ids: List[str]):
        """_summary_

        Args:
            rxn_ids (List[str]): _description_
        """
        self.model.remove_reactions(self.resolve_many(rxn_ids, "reaction"))
        self._invalidate_id_cache()

    def update_stoichiometry(self, rxn_id, stoichiometry):
        """Updates the stoichiometry of a reaction by creating a
        new reaction with a same id and the new stoichiometry

        :param rxn_id: Reaction identifier
        :type rxn_id: str
        :param stoichiometry: the new stoichiometry
        :type stoichiometry: dict[str:float]
        """
        rxn_id = self.resolve(rxn_id, "reaction")
        rxn = self.model.reactions.get_by_id(rxn_id)
        objective = self.objective.get(rxn_id, 0)
        self.add_reaction(
            rxn_id,
            name=rxn.name,
            stoichiometry=stoichiometry,
            lb=rxn.lower_bound,
            ub=rxn.upper_bound,
            gpr=rxn.gene_reaction_rule,
            annotations=rxn.annotation,
            objective=objective,
            replace=True,
        )

    def get_uptake_reactions(self):
        """
        :returns: The list of uptake reactions.
        """
        drains = self.get_exchange_reactions()
        rxns = [
            r
            for r in drains
            if self.model.reactions.get_by_id(r).reversibility
            or (
                (
                    self.model.reactions.get_by_id(r).lower_bound is None
                    or self.model.reactions.get_by_id(r).lower_bound < 0
                )
                and len(self.model.reactions.get_by_id(r).reactants) > 0
            )
            or (
                (
                    self.model.reactions.get_by_id(r).upper_bound is None
                    or self.model.reactions.get_by_id(r).upper_bound > 0
                )
                and len(self.model.reactions.get_by_id(r).products) > 0
            )
        ]
        return rxns

    def reverse_reaction(self, reaction_id):
        """
        Identify if a reaction is reversible and returns the reverse reaction if it is the case.

        :param reaction_id: A reaction identifier.
        :returns: A reaction identifier or None.

        """
        reaction_id = self.resolve(reaction_id, "reaction")
        rxn = self.model.reactions.get_by_id(reaction_id)
        reactions = self.reactions
        if rxn.lower_bound < 0:
            return reaction_id
        else:
            rev = rxn.reverse_id
            if rev in reactions:
                return rev
        return None

    def metabolite_reaction_lookup(self, force_recalculate=False):
        """Return the network topology as a nested map from metabolite to reaction to coefficient.
        :return: a dictionary lookup table
        """

        if not self._m_r_lookup or force_recalculate:
            self._m_r_lookup = OrderedDict([(m_id, OrderedDict()) for m_id in self.metabolites])
            for _, reaction in enumerate(self.model.reactions):
                for m, coeff in reaction.metabolites.items():
                    self._m_r_lookup[m.id][reaction.id] = coeff

        return self._m_r_lookup

    def metabolite_elements(self, metabolite_id):
        return self.model.metabolites.get_by_id(self.resolve(metabolite_id, "metabolite")).elements

    def get_reaction_bounds(self, reaction_id):
        """
        Returns the bounds for a given reaction.

        :param reaction_id: str, reaction ID
        :return: lb(s), ub(s), tuple

        """
        lb, ub = self.model.reactions.get_by_id(self.resolve(reaction_id, "reaction")).bounds
        # return lb if lb > -np.inf else ModelConstants.REACTION_LOWER_BOUND,\
        #    ub if ub < np.inf else ModelConstants.REACTION_UPPER_BOUND
        return lb, ub

    def set_reaction_bounds(self, reaction_id, lb, ub, track=True):
        """
        Sets the bounds for a given reaction.
        :param reaction_id: str, reaction ID
        :param float lb: lower bound
        :param float ub: upper bound
        :param bool track: if the changes are to be logged. Default True
        """
        reaction_id = self.resolve(reaction_id, "reaction")
        if track:
            if reaction_id in self.get_uptake_reactions():
                self._environmental_conditions[reaction_id] = (lb, ub)
            else:
                self._constraints[reaction_id] = (lb, ub)
        rxn = self.model.reactions.get_by_id(reaction_id)
        rxn.bounds = (lb, ub)

    def find_bounds(self):
        """
        Return the median upper and lower bound of the metabolic model.
        Bounds can vary from model to model. Cobrapy defaults to (-1000, 1000).
        """
        lower_bounds = np.asarray([rxn.lower_bound for rxn in self.model.reactions], dtype=float)
        upper_bounds = np.asarray([rxn.upper_bound for rxn in self.model.reactions], dtype=float)
        lower_bound = np.nanmedian(lower_bounds[lower_bounds != 0.0])
        upper_bound = np.nanmedian(upper_bounds[upper_bounds != 0.0])
        if np.isnan(lower_bound):
            LOGGER.warning("Could not identify a median lower bound.")
            lower_bound = ModelConstants.REACTION_LOWER_BOUND
        if np.isnan(upper_bound):
            LOGGER.warning("Could not identify a median upper bound.")
            upper_bound = ModelConstants.REACTION_UPPER_BOUND
        return lower_bound, upper_bound

    def find_unconstrained_reactions(self):
        """Return list of reactions that are not constrained at all."""
        lower_bound, upper_bound = self.find_bounds()
        return [
            rxn.id for rxn in self.model.reactions if rxn.lower_bound <= lower_bound and rxn.upper_bound >= upper_bound
        ]

    # The simulator
    def simulate(
        self,
        objective: Dict[str, float] = None,
        method: Union[SimulationMethod, str, Callable] = SimulationMethod.FBA,
        maximize: bool = True,
        constraints: Dict[str, Union[float, Tuple[float, float]]] = None,
        reference: Dict[str, float] = None,
        scalefactor: float = None,
        solver=None,
        slim: bool = False,
        shadow_prices: bool = False,
        **kwargs,
    ) -> SimulationResult:
        """
        Simulates a phenotype when applying a set constraints using the specified method.

        :param dic objective: The simulation objective. If none, the model objective is considered.
        :param method: The SimulationMethod (FBA, pFBA, lMOMA, etc ...)
        :param boolean maximize: The optimization direction.
        :param dic constraints: A dictionary of constraints to be applied to the model.
        :param dic reference: A dictionary of reaction flux values.
        :param float scalefactor: A positive scaling factor for the solver. Default None.

        """
        # ids may be given with or without the SBML prefix
        objective = self._resolve_objective(objective)
        constraints = self.resolve_keys(constraints)
        reference = self._resolve_mapping(reference)

        if callable(method):
            return self._simulate_callable(
                method, objective=objective, maximize=maximize, constraints=constraints, **kwargs
            )

        if not objective:
            objective = self.model.objective

        # MOMA, lMOMA and ROOM compare against the wild type, so the reference is computed
        # before the simulation constraints are applied (cobra would otherwise use the mutant's pFBA).
        ref_solution, partial_reference = None, False
        if method in (SimulationMethod.lMOMA, SimulationMethod.MOMA, SimulationMethod.ROOM):
            if reference is None:
                reference = self.reference
            # reference ids may be given with or without the SBML prefix, as for constraints
            reference = self._resolve_mapping(reference)
            unknown = [r_id for r_id in reference if r_id not in self.model.reactions]
            if unknown:
                raise ValueError(f"Reference reactions not in the model: {unknown[:10]}")
            partial_reference = _is_partial(self.model, reference)
            ref_solution = Solution(None, "optimal", fluxes=pd.Series(reference, dtype=float))

        simul_constraints = {}
        if constraints:
            if not self._allow_env_changes:
                simul_constraints.update(
                    {k: v for k, v in constraints.items() if k not in list(self._environmental_conditions.keys())}
                )
            else:
                simul_constraints.update(constraints)

        with self.model as model:
            if isinstance(objective, dict):
                model.objective = {model.reactions.get_by_id(r_id): coef for r_id, coef in objective.items()}
            else:
                model.objective = objective
            for rxn in list(simul_constraints.keys()):
                reac = model.reactions.get_by_id(rxn)

                # constraints defined as a tuple (lower_bound, upper_bound) or
                # as a single float value
                if isinstance(simul_constraints.get(rxn), tuple):
                    reac.bounds = (simul_constraints.get(rxn)[0], simul_constraints.get(rxn)[1])
                else:
                    reac.bounds = (simul_constraints.get(rxn), simul_constraints.get(rxn))

            # If working directly over optlang use 'max' and 'min'
            # such is the case with pytfa.core.Model.
            objective_sense = self._MAX_STR if maximize else self._MIN_STR

            try:
                # FBA
                if method == SimulationMethod.FBA:
                    if slim:
                        # slim_optimize uses the objective direction; set it from maximize (reverted on exit)
                        model.objective_direction = "max" if maximize else "min"
                        solution = model.slim_optimize()
                        if model.solver.status != "optimal":
                            solution = None  # the infeasible-result convention (Simulator.simulate)
                    else:
                        solution = model.optimize(objective_sense=objective_sense)
                # pFBA
                elif method == SimulationMethod.pFBA:
                    # the fraction of the optimum kept is shared by every backend
                    solution = pfba(model, fraction_of_optimum=ModelConstants.PFBA_OBJ_FRAC)

                # MOMA/lMOMA; a reference that leaves out reactions counts only the reactions it has (as REFRAMED)
                elif method == SimulationMethod.lMOMA or method == SimulationMethod.MOMA:
                    linear = True if method == SimulationMethod.lMOMA else False
                    if partial_reference:
                        _add_partial_moma(model, reference, linear)
                        solution = model.optimize()
                    else:
                        solution = moma(model, solution=ref_solution, linear=linear)

                # ROOM
                elif method == SimulationMethod.ROOM:
                    if partial_reference:
                        _add_partial_room(model, reference)
                        solution = model.optimize()
                    else:
                        solution = room(model, solution=ref_solution)

                # Special case in which only the simulation context is required without any simulation result
                elif method == SimulationMethod.NONE:
                    solution = Solution(None, "unknown", None)

                else:
                    raise Exception("Unknown method to perform the simulation.")
            except OptimizationError:
                # pfba, moma and room raise on infeasible problems; report the status like FBA does
                solution = Solution(None, model.solver.status, None)

        if slim and method == SimulationMethod.FBA:
            return solution

        status = self.__status_mapping[solution.status]
        solved = status in (SStatus.OPTIMAL, SStatus.SUBOPTIMAL)
        if slim:
            return solution.objective_value if solved else None

        fluxes = solution.fluxes.to_dict(into=OrderedDict) if solved and solution.fluxes is not None else None
        prices = solution.shadow_prices if solved else None
        result = SimulationResult(
            self,
            solution.objective_value if solved else None,
            fluxes=fluxes,
            status=status,
            envcond=self.environmental_conditions,
            model_constraints=self._constraints.copy(),
            simul_constraints=dict(constraints) if constraints else None,
            maximize=self._optimised_sense(method, maximize),
            method=method,
            shadow_prices=prices.to_dict(into=OrderedDict) if prices is not None else None,
        )
        return result

    def FVA(
        self,
        reactions: Union[List[str], None] = None,
        obj_frac: float = 0.9,
        constraints: Dict[str, Union[float, Tuple[float, float]]] = None,
        loopless: bool = False,
        solver=None,
        format: str = "dict",
        processes: int = 1,
    ) -> Union[dict, "DataFrame"]:
        """ Flux Variability Analysis (FVA).

        :param model: An instance of a constraint-based model.
        :param int processes: Worker processes for cobra's FVA (default 1, no pool). With the 'spawn' start method
            (the macOS/Windows default) a value > 1 needs the calling script to use an ``if __name__ == "__main__":``
            guard, otherwise it hangs.
        :param float obj_frac: The minimum fraction of the maximum growth rate (default 0.9). Requires that the \
            objective value is at least the fraction times maximum objective value. A value of 0.85 for instance \
            means that the objective has to be at least at 85% percent of its maximum.
        :param list reactions: List of reactions to analyze (default: all).
        :param dic constraints: Additional constraints (optional).
        :param boolean loopless: Run looplessFBA internally (very slow) (default: false).
        :param solver: A pre-instantiated solver instance (optional).
        :param format: The return format: 'dict' returns a dictionary, 'df' returns a pandas DataFrame.
        :returns: Flux variation ranges. Returns dict[str, list[float, float]] if format='dict',
                  or pandas.DataFrame with columns ['Reaction ID', 'Minimum', 'Maximum'] if format='df'.

        """
        from cobra.flux_analysis.variability import flux_variability_analysis

        constraints = self.resolve_keys(constraints)
        simul_constraints = {}
        if constraints:
            simul_constraints.update(
                {k: v for k, v in constraints.items() if k not in list(self._environmental_conditions.keys())}
            )

        if reactions is None:
            _reactions = self.reactions
        elif isinstance(reactions, str):
            _reactions = [self.resolve(reactions, "reaction")]
        elif isinstance(reactions, list):
            _reactions = self.resolve_many(reactions, "reaction")
        else:
            raise ValueError("Invalid reactions.")

        with self.model as model:

            if simul_constraints:
                for rxn in list(simul_constraints.keys()):
                    reac = model.reactions.get_by_id(rxn)
                    if isinstance(simul_constraints.get(rxn), tuple):
                        reac.bounds = (simul_constraints.get(rxn)[0], simul_constraints.get(rxn)[1])
                    else:
                        reac.bounds = (simul_constraints.get(rxn), simul_constraints.get(rxn))

            df = flux_variability_analysis(
                model, reaction_list=_reactions, loopless=loopless, fraction_of_optimum=obj_frac, processes=processes
            )

        variability = {}
        for r_id in _reactions:
            variability[r_id] = [float(df.loc[r_id].iloc[0]), float(df.loc[r_id].iloc[1])]

        if format == "df":
            import pandas as pd

            e = variability.items()
            f = [[a, b, c] for a, [b, c] in e]
            # the same layout on every backend: 'Reaction ID' is a column, with a default index
            return pd.DataFrame(f, columns=["Reaction ID", "Minimum", "Maximum"])
        else:
            return variability

    def set_objective(self, reaction_id: str):
        self.model.objective = self.resolve(reaction_id, "reaction")

    def create_empty_model(self, model_id: str):
        return Simulation(Model(model_id))


class GeckoSimulation(Simulation):
    """
    Simulator for geckopy.gecko.GeckoModel
    """

    def __init__(
        self,
        model,
        envcond=None,
        constraints=None,
        solver=None,
        reference=None,
        reset_solver=ModelConstants.RESET_SOLVER,
        protein_prefix=None,
    ):
        try:
            from geckopy.gecko import GeckoModel

            if not isinstance(model, GeckoModel):
                raise ValueError("The model is not an instance of geckopy.gecko.GeckoModel")
        except ImportError:
            raise RuntimeError("The geckopy package is not installed.")

        super(GeckoSimulation, self).__init__(model, envcond, constraints, solver, reference, reset_solver)
        self.protein_prefix = protein_prefix if protein_prefix else self.to_native("draw_prot_", "reaction")
        self._essential_proteins = None
        self._protein_rev_reactions = None
        self._prot_react = None

    @property
    def proteins(self):
        return list(self.model.proteins)

    @property
    def protein_pool_exchange(self):
        return self.model.protein_pool_exchange

    def essential_proteins(self, min_growth=0.01):
        if self._essential_proteins is not None:
            return self._essential_proteins
        wt_solution = self.simulate()
        wt_growth = wt_solution.objective_value
        self._essential_proteins = []
        proteins = self.model.proteins
        for p in tqdm(proteins):
            rxn = "{}{}".format(self.protein_prefix, p)
            res = self.simulate(constraints={rxn: 0})
            if res:
                if (
                    res.status == SStatus.OPTIMAL and res.objective_value < wt_growth * min_growth
                ) or res.status == SStatus.INFEASIBLE:
                    self._essential_proteins.append(p)
        return self._essential_proteins

    def map_prot_react(self):
        if self._prot_react is None:
            self._prot_react = dict()
            for p in self.proteins:
                self._prot_react[p] = []

            for r_id in self.reactions:
                rxn = self.model.reactions.get_by_id(r_id)
                lsub = rxn.reactants
                for m in lsub:
                    if "prot_" in m.id:
                        p = m.id[5:-2]
                        try:
                            reaction_list = self._prot_react[p]
                            reaction_list.append(r_id)

                        except KeyError:
                            # Protein not in mapping
                            pass
        return self._prot_react

    def protein_reactions(self, protein):
        """
        Returns the list of reactions associated to a protein.
        """
        mapper = self.map_prot_react()
        return mapper[protein]

    def get_protein(self, p_id):
        res = {"Protein": p_id, "reactions": self.protein_reactions(p_id)}
        return AttrDict(res)

    def find_proteins(self, pattern=None, sort=False):
        """A user friendly method to find proteins in the model.

        :param pattern: The pattern which can be a regular expression, defaults to None in which case all entries are listed.
        :type pattern: str, optional
        :param sort: if the search results should be sorted, defaults to False
        :type sort: bool, optional

        :return: the search results
        :rtype: pandas dataframe
        """
        values = self.proteins
        if pattern:
            import re

            if isinstance(pattern, list):
                patt = "|".join(pattern)
                re_expr = re.compile(patt)
            else:
                re_expr = re.compile(pattern)
            values = [x for x in values if re_expr.search(x) is not None]
        if sort:
            values.sort()

        import pandas as pd

        data = [self.get_protein(x) for x in values]
        df = pd.DataFrame(data)
        df = df.set_index(df.columns[0])
        return df

    def reverse_reaction(self, reaction_id):
        """
        Returns the list of reactions associated with a protein.

        :params (str) protein: The protein identifier.
        :returns: A list of reaction identifiers
        """
        f, d = zip(*self.protein_rev_reactions.values())
        if reaction_id in f:
            return d[f.index(reaction_id)]
        elif reaction_id in d:
            return f[d.index(reaction_id)]
        else:
            return None

    @property
    def protein_rev_reactions(self):
        """
        Identify if a reaction is reversible and returns the
        reverse reaction if it is the case.

        :param (str) reaction_id: The reaction identifier.
        :returns: A reaction identifier or None.
        """
        if not self._protein_rev_reactions:
            proteins = self.model.proteins
            reactions = self.reactions
            in_sub = {}
            for p in proteins:
                sub = []
                for r_id in reactions:
                    rxn = self.model.reactions.get_by_id(r_id)
                    lsub = rxn.reactants
                    for m in lsub:
                        if p in m.id:
                            sub.append(r_id)
                in_sub[p] = sub
            pairs = {}
            for k, s in in_sub.items():
                revs = [r for r in s if "_REV" in r]
                if len(revs) > 0:
                    for r in revs:
                        la = [a for a in s if r.replace("_REV", "") == a]
                        la.append(r)
                        if len(la) == 2:
                            if k in pairs.keys():
                                pairs[k].append((la[0], la[1]))
                            else:
                                pairs[k] = [(la[0], la[1])]
            self._protein_rev_reactions = pairs
        return self._protein_rev_reactions

    def get_Kcats(self, protein: str):
        """
        Returns a dictionary of reactions and respective Kcat for a
        specific protein/enzyme·

        :params (str) protein: the protein identifier.

        :returns: A dictionary of reactions and respective Kcat values.
        """
        import re

        re_expr = re.compile(f"{protein}_")
        values = [x for x in self.metabolites if re_expr.search(x) is not None]
        if len(values) == 1:
            m_r = self.metabolite_reaction_lookup()
            r_d = m_r[values[0]]
            return {k: -1 / v for k, v in r_d.items() if self.protein_prefix not in k}
        elif len(values) > 1:
            raise ValueError(f"More than one protein match {values}")
        else:
            raise ValueError(f"Protein {protein} not founded.")

    def set_Kcat(self, protein, reaction, kcat):
        """Alters an enzyme kcat for a given reaction.

        :param protein: The protein identifier
        :type protein: str
        :param reaction: The reaction identifier
        :type reaction: str
        :param invkcat: The inverse kcat value, a real positive value.
        :type invkcat: float
        """
        if kcat <= 0:
            raise ValueError("kcat value needs to be positive.")

        rxn = self.model.reactions.get_by_id(self.resolve(reaction, "reaction"))
        m = None
        for met in rxn.metabolites:
            if protein in met.id:
                m = met
                break
        if m is not None:
            rxn.subtract_metabolites({m: -1 / kcat})
        else:
            LOGGER.warn(f"Could not identify {protein} " f"protein specie in reaction {reaction}")
