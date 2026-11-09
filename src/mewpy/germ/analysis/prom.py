"""
Probabilistic Regulation of Metabolism (PROM) - Clean Implementation

This module implements PROM using RegulatoryExtension only.
No backwards compatibility with legacy GERM models.
"""

from typing import TYPE_CHECKING, Any, Dict, Sequence, Tuple, Union

import pandas as pd

from mewpy.germ.analysis.fba import _RegulatoryAnalysisBase
from mewpy.germ.models.regulatory_extension import RegulatoryExtension
from mewpy.germ.solution import KOSolution
from mewpy.solvers.solution import Solution, Status
from mewpy.solvers.solver import Solver
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.germ.models import RegulatoryModel
    from mewpy.germ.variables import Gene, Regulator, Target


def _run_and_decode_solver(lp, additional_constraints: Dict[str, Tuple[float, float]] = None, **kwargs):
    if not additional_constraints:
        additional_constraints = {}

    # merge into a new dict: always apply the additional constraints, never mutate the caller's constraints
    kwargs["constraints"] = {**(kwargs.get("constraints") or {}), **additional_constraints}

    solution = lp.solver.solve(**kwargs)
    if solution.status == Status.OPTIMAL:
        return solution.fobj
    else:
        return


class PROM(_RegulatoryAnalysisBase):
    """
    Probabilistic Regulation of Metabolism (PROM) using RegulatoryExtension.

    PROM predicts the growth phenotype and flux response after transcriptional
    perturbation, given a metabolic and regulatory network. PROM introduces
    probabilities to represent gene states and gene-transcription factor interactions.

    For more details: https://doi.org/10.1073/pnas.1005139107
    """

    def __init__(
        self, model: RegulatoryExtension, solver: Union[str, Solver] = None, build: bool = False, attach: bool = False
    ):
        """
        Initialize PROM with a RegulatoryExtension model.

        :param model: A RegulatoryExtension instance wrapping a simulator
        :param solver: The solver to be used. If None, a new instance will be created from the default solver.
        :param build: If True, the linear problem will be built upon initialization.
        :param attach: If True, the linear problem will be attached to the model.
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)
        self.method = "PROM"

    def _build(self):
        """
        Build the PROM problem.

        It also builds a regular FBA problem to be used for the growth prediction.
        The linear problem is then loaded into the solver.
        """
        self._build_mass_constraints()
        self._linear_objective = dict(self.model.objective)
        self._minimize = False

    def _max_rates(self, solver_kwargs: Dict[str, Any], reactions: Sequence[str] = None):
        """Compute maximum rates using FVA, for the given reactions (default: all reactions)."""
        # Wild-type reference
        reference = self.solver.solve(
            **{**solver_kwargs, "linear": self._linear_objective, "minimize": self._minimize, "get_values": True}
        )
        if reference.status != Status.OPTIMAL:
            raise RuntimeError("The solver did not find an optimal solution for the wild-type conditions.")
        reference = reference.values.copy()
        reference_constraints = {key: (reference[key] * 0.99, reference[key]) for key in self._linear_objective}

        # FVA of the reaction at fraction of 0.99 (for wild-type growth rate)
        rates = {}
        for reaction in self.model.reactions if reactions is None else reactions:
            min_rxn = _run_and_decode_solver(
                self,
                additional_constraints=reference_constraints,
                **{**solver_kwargs, "get_values": False, "linear": {reaction: 1}, "minimize": True},
            )
            max_rxn = _run_and_decode_solver(
                self,
                additional_constraints=reference_constraints,
                **{**solver_kwargs, "get_values": False, "linear": {reaction: 1}, "minimize": False},
            )

            reference_rate = reference[reaction]

            # Handle None values from infeasible solutions
            if min_rxn is None:
                min_rxn = reference_rate
            if max_rxn is None:
                max_rxn = reference_rate

            if reference_rate < 0:
                value = min((min_rxn, max_rxn, reference_rate))
            elif reference_rate > 0:
                value = max((min_rxn, max_rxn, reference_rate))
            else:
                value = max((abs(min_rxn), abs(max_rxn), abs(reference_rate)))

            if abs(value) < ModelConstants.TOLERANCE:
                value = 0.0

            rates[reaction] = value

        return rates

    def _gene_reactions(self, gene_id: str) -> list:
        """Reaction ids of a metabolic gene (RegulatoryExtension or legacy GERM model)."""
        if hasattr(self.model, "get_gene"):
            return list(self.model.get_gene(gene_id).reactions)
        gene = self.model.genes[gene_id]
        return list(gene.reactions) if hasattr(gene, "reactions") else []

    def _optimize_ko(
        self,
        probabilities: Dict[Tuple[str, str], float],
        regulator: Union["Gene", "Regulator"],
        reference: Dict[str, float],
        max_rates: Dict[str, float],
        to_solver: bool = False,
        solver_kwargs: Dict[str, Any] = None,
    ):
        """Optimize with regulator knockout."""
        solver_constrains = solver_kwargs.get("constraints") or {}

        # Reaction bounds: model bounds overridden by the caller's constraints. Legacy GERM models yield
        # Reaction objects and RegulatoryExtension yields ids; constraints are always keyed by reaction id.
        prom_constraints = {}
        for rxn in self.model.yield_reactions():
            rxn_id = getattr(rxn, "id", rxn)
            rxn_data = self._get_reaction(rxn_id)
            prom_constraints[rxn_id] = (
                rxn_data.get("lb", ModelConstants.REACTION_LOWER_BOUND),
                rxn_data.get("ub", ModelConstants.REACTION_UPPER_BOUND),
            )
        for key, bounds in solver_constrains.items():
            prom_constraints[key] = tuple(bounds) if isinstance(bounds, (tuple, list)) else (bounds, bounds)
        base_bounds = dict(prom_constraints)

        genes = self.model.genes
        state = {gene: 1 for gene in genes}

        # If the regulator to be KO is a metabolic gene, the associated reactions are KO too
        # Check if regulator ID exists in model genes list
        if regulator.id in genes:
            # Handle both APIs: RegulatoryExtension has get_gene(), legacy has genes dict
            if hasattr(self.model, "get_gene"):
                gene_data = self.model.get_gene(regulator.id)
                reactions_list = gene_data.reactions
            else:
                gene_obj = genes[regulator.id]
                reactions_list = gene_obj.reactions if hasattr(gene_obj, "reactions") else []

            for rxn_id in reactions_list:
                prom_constraints[rxn_id] = (
                    -ModelConstants.TOLERANCE,
                    ModelConstants.TOLERANCE,
                )

        # Find the target genes of the deleted regulator
        target_reactions = {}
        for target in regulator.yield_targets():
            # Check if this target corresponds to a metabolic gene
            if target.id in genes:
                state[target.id] = 0

                # Handle both APIs
                if hasattr(self.model, "get_gene"):
                    gene_data = self.model.get_gene(target.id)
                    reactions_list = gene_data.reactions
                else:
                    gene_obj = genes[target.id]
                    reactions_list = gene_obj.reactions if hasattr(gene_obj, "reactions") else []

                for rxn_id in reactions_list:
                    target_reactions[rxn_id] = rxn_id  # Store ID, not object

        # GPR evaluation using changed gene state
        inactive_reactions = {}
        for rxn_id in target_reactions.keys():
            # Handle both APIs for accessing GPR
            if hasattr(self.model, "get_parsed_gpr"):
                # RegulatoryExtension: has get_parsed_gpr() method
                gpr = self.model.get_parsed_gpr(rxn_id)
            else:
                # Legacy model: access reaction's gpr attribute directly
                rxn = self.model.reactions[rxn_id]
                gpr = rxn.gpr

            if gpr.is_none:
                continue

            if gpr.evaluate(values=state):
                continue

            inactive_reactions[rxn_id] = rxn_id  # Store ID, not object

        # For each target regulated by the regulator
        for target in regulator.yield_targets():
            # Check if target is a metabolic gene
            if target.id not in genes:
                continue

            # Composed key for interactions_probabilities
            target_regulator = (target.id, regulator.id)

            if target_regulator not in probabilities:
                continue

            interaction_probability = probabilities[target_regulator]

            # Get reactions for this gene - handle both APIs
            if hasattr(self.model, "get_gene"):
                gene_data = self.model.get_gene(target.id)
                reactions_list = gene_data.reactions
            else:
                gene_obj = genes[target.id]
                reactions_list = gene_obj.reactions if hasattr(gene_obj, "reactions") else []

            # For each reaction associated with this single target
            for rxn_id in reactions_list:
                if rxn_id not in inactive_reactions:
                    continue

                if interaction_probability >= 1:
                    continue

                # Reaction old bounds
                rxn_lb, rxn_ub = tuple(prom_constraints[rxn_id])

                # Probability flux
                probability_flux = max_rates[rxn_id] * interaction_probability

                # Wild-type flux value
                wt_flux = reference[rxn_id]

                # Reaction bounds before the knockout (model bounds or the caller's constraints)
                reaction_lower_bound, reaction_upper_bound = base_bounds[rxn_id]

                # Update flux bounds according to probability flux
                if wt_flux < 0:
                    rxn_lb = max((reaction_lower_bound, probability_flux, rxn_lb))
                    rxn_lb = min((rxn_lb, -ModelConstants.TOLERANCE))
                elif wt_flux > 0:
                    rxn_ub = min((reaction_upper_bound, probability_flux, rxn_ub))
                    rxn_ub = max((rxn_ub, ModelConstants.TOLERANCE))
                else:
                    # If it is zero, the reaction is not changed
                    continue

                prom_constraints[rxn_id] = (rxn_lb, rxn_ub)

        solution = self.solver.solve(
            **{
                **solver_kwargs,
                "linear": self._linear_objective,
                "minimize": self._minimize,
                "get_values": True,
                "constraints": prom_constraints,
            }
        )

        if to_solver:
            return solution

        return Solution.from_solver(method=self.method, solution=solution, model=self.model, minimize=self._minimize)

    def _optimize(
        self,
        initial_state: Dict[Tuple[str, str], float] = None,
        regulators: Sequence[Union["Gene", "Regulator"]] = None,
        to_solver: bool = False,
        solver_kwargs: Dict[str, Any] = None,
    ) -> Union[Dict[str, Solution], Dict[str, Solution]]:
        """Internal optimization method."""
        # Wild-type reference
        solver_kwargs["get_values"] = True
        reference = self.solver.solve(linear=self._linear_objective, minimize=self._minimize, **solver_kwargs)
        if reference.status != Status.OPTIMAL:
            raise RuntimeError("The solver did not find an optimal solution for the wild-type conditions.")
        reference = reference.values.copy()

        # Max and min fluxes, only of the reactions the knockouts can restrict (those of the regulators'
        # metabolic targets): the other rates are never used, and FVA costs two solves per reaction
        genes = set(self.model.genes)
        regulated = {
            rxn_id
            for regulator in regulators
            for target in regulator.yield_targets()
            if target.id in genes
            for rxn_id in self._gene_reactions(target.id)
        }
        max_rates = self._max_rates(solver_kwargs=solver_kwargs, reactions=sorted(regulated))

        # Single regulator knockout
        if len(regulators) == 1:
            ko_solution = self._optimize_ko(
                probabilities=initial_state,
                regulator=regulators[0],
                reference=reference,
                max_rates=max_rates,
                to_solver=to_solver,
                solver_kwargs=solver_kwargs,
            )
            return {regulators[0].id: ko_solution}

        # Multiple regulator knockouts
        kos = {}
        for regulator in regulators:
            ko_solution = self._optimize_ko(
                probabilities=initial_state,
                regulator=regulator,
                reference=reference,
                max_rates=max_rates,
                to_solver=to_solver,
                solver_kwargs=solver_kwargs,
            )
            kos[regulator.id] = ko_solution
        return kos

    def optimize(
        self,
        initial_state: Dict[Tuple[str, str], float] = None,
        regulators: Union[str, Sequence["str"]] = None,
        to_solver: bool = False,
        solver_kwargs: Dict[str, Any] = None,
    ) -> Union[KOSolution, Dict[str, Solution]]:
        """
        Solve the PROM linear problem.

        :param initial_state: Dictionary with the probabilities of the interactions
                            between the regulators and the targets.
        :param regulators: List of regulators to be knocked out. If None, all regulators are knocked out.
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False.
        :param solver_kwargs: Solver parameters to be set temporarily.
        :return: A KOSolution instance or a list of SolverSolution instances if to_solver is True.
        """
        # Build solver if out of sync
        if not self.synchronized:
            self.build()

        if not initial_state:
            initial_state = {}

        # Ids may be given in either SBML form or as original TRN ids
        initial_state = self._ids.keys(initial_state)

        if not regulators:
            # RegulatoryExtension yields (id, regulator) tuples, legacy models yield regulators
            regulators = [regulator for _, regulator in self._get_regulators()]
        else:
            if isinstance(regulators, str):
                regulators = [regulators]
            regulators = self._ids.ids(regulators)
            # Get regulator objects - handle both APIs
            if hasattr(self.model, "get_regulator"):
                # RegulatoryExtension API
                regulators = [self.model.get_regulator(regulator) for regulator in regulators]
            else:
                # Legacy model API - regulators is a dict
                regulators = [self.model.regulators[regulator] for regulator in regulators]

        # a copy: the caller's dict is not mutated. 'linear'/'minimize' set a temporary objective
        solver_kwargs = dict(self._ids.solver_kwargs(solver_kwargs)) if solver_kwargs else {}
        linear, minimize = self._objective_from(solver_kwargs)
        objective, direction = self._linear_objective, self._minimize
        self._linear_objective, self._minimize = linear, minimize

        # Concrete optimize
        try:
            solutions = self._optimize(
                initial_state=initial_state, regulators=regulators, to_solver=to_solver, solver_kwargs=solver_kwargs
            )
        finally:
            self._linear_objective, self._minimize = objective, direction

        if to_solver:
            return solutions

        return KOSolution(solutions)


# ----------------------------------------------------------------------------------------------------------------------
# Probability of Target-Regulator interactions
# ----------------------------------------------------------------------------------------------------------------------
def target_regulator_interaction_probability(
    model: Union[RegulatoryExtension, "RegulatoryModel"],
    expression: pd.DataFrame,
    binary_expression: pd.DataFrame,
) -> Tuple[Dict[Tuple[str, str], float], Dict[Tuple[str, str], float]]:
    """
    Compute the conditional probability of a target gene being active when the regulator is inactive.

    Uses the formula:
        P(target = 1 | regulator = 0) = count(target = 1, regulator = 0) / # samples

    This probability is computed for each combination of target-regulator.
    This method is used in PROM analysis.

    :param model: A RegulatoryExtension or legacy GERM RegulatoryModel instance
    :param expression: Quantile preprocessed expression matrix
    :param binary_expression: Quantile preprocessed expression matrix binarized
    :return: Dictionary with the conditional probability of a target gene being active when the regulator is inactive,
            Dictionary with missed interactions
    """
    try:
        # noinspection PyPackageRequirements
        from scipy.stats import ks_2samp
    except ImportError:
        raise ImportError(
            "The package scipy is not installed. "
            "To compute the probability of target-regulator interactions, please install scipy "
            "(pip install scipy)."
        )

    missed_interactions = {}
    interactions_probabilities = {}

    # Handle both legacy GERM models (yield Interaction) and RegulatoryExtension (yield tuples)
    interactions_gen = model.yield_interactions()
    first_item = next(interactions_gen, None)
    if first_item is None:
        return interactions_probabilities, missed_interactions

    # Check if yielded items are tuples (RegulatoryExtension) or objects (legacy)
    if isinstance(first_item, tuple):
        # RegulatoryExtension: yields (id, Interaction) tuples
        def _get_interactions():
            yield first_item[1]  # Unwrap first item
            for _, interaction in interactions_gen:
                yield interaction

    else:
        # Legacy GERM: yields Interaction objects directly
        def _get_interactions():
            yield first_item
            yield from interactions_gen

    for interaction in _get_interactions():
        target = interaction.target

        if not interaction.regulators or target.id not in expression.index:
            missed_interactions[(target.id, target.id)] = 1
            interactions_probabilities[(target.id, target.id)] = 1
            continue

        target_expression = expression.loc[target.id]
        target_binary = binary_expression.loc[target.id]

        for regulator in interaction.yield_regulators():
            if regulator.id not in expression.index:
                missed_interactions[(target.id, regulator.id)] = 1
                interactions_probabilities[(target.id, regulator.id)] = 1
                continue

            regulator_binary = binary_expression.loc[regulator.id]

            target_expression_1_regulator = target_expression[regulator_binary == 1]
            target_expression_0_regulator = target_expression[regulator_binary == 0]

            if len(target_expression_1_regulator) == 0 and len(target_expression_0_regulator) == 0:
                missed_interactions[(target.id, regulator.id)] = 1
                interactions_probabilities[(target.id, regulator.id)] = 1
                continue

            _, p_val = ks_2samp(target_expression_1_regulator, target_expression_0_regulator)
            if p_val < 0.05:
                target_binary_0_regulator = target_binary[regulator_binary == 0]
                probability = sum(target_binary_0_regulator) / len(target_binary_0_regulator)
                interactions_probabilities[(target.id, regulator.id)] = probability
                missed_interactions[(target.id, regulator.id)] = 0
            else:
                missed_interactions[(target.id, regulator.id)] = 1
                interactions_probabilities[(target.id, regulator.id)] = 1

    return interactions_probabilities, missed_interactions
