"""
CoRegFlux - Clean Implementation

This module implements CoRegFlux using RegulatoryExtension only.
No backwards compatibility with legacy GERM models.
"""

from typing import TYPE_CHECKING, Dict, List, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from mewpy.germ.analysis.analysis_utils import (
    CoRegBiomass,
    CoRegMetabolite,
    CoRegResult,
    biomass_yield_to_rate,
    build_biomass,
    build_metabolites,
    gene_state_constraints,
    metabolites_constraints,
    system_state_update,
)
from mewpy.germ.analysis.fba import _RegulatoryAnalysisBase
from mewpy.germ.models.regulatory_extension import RegulatoryExtension
from mewpy.germ.solution import DynamicSolution
from mewpy.germ.variables import Gene, Target
from mewpy.solvers.solution import Solution, Status
from mewpy.solvers.solver import Solver
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.germ.models import MetabolicModel, Model, RegulatoryModel


def _run_and_decode(lp, additional_constraints=None, solver_kwargs=None):
    """Run solver and decode results."""
    if not solver_kwargs:
        solver_kwargs = {}

    if additional_constraints:
        solver_kwargs["constraints"] = {**solver_kwargs.get("constraints", {}), **additional_constraints}

    solution = lp.solver.solve(linear=lp._linear_objective, minimize=lp._minimize, **solver_kwargs)

    reactions = lp.model.reactions

    if not solution.values:
        return {rxn: 0 for rxn in reactions}, 0

    return solution.values, solution.fobj


def result_to_solution(result: CoRegResult, model: RegulatoryExtension, to_solver: bool = False) -> Solution:
    """
    Convert a CoRegResult object to a Solution object.

    :param result: the CoRegResult object
    :param model: the RegulatoryExtension model
    :param to_solver: if True, it returns a Solution object
    :return: the Solution object
    """
    if to_solver:
        return Solution(status=Status.OPTIMAL, fobj=result.objective_value, values=result.values)

    solution = Solution(
        objective_value=result.objective_value,
        values=result.values,
        status=Status.OPTIMAL,
        method="CoRegFlux",
        model=model,
    )

    solution.metabolites = {key: met.concentration for key, met in result.metabolites.items()}
    solution.biomass = result.biomass.biomass_yield
    solution.constraints = result.constraints
    return solution


class CoRegFlux(_RegulatoryAnalysisBase):
    """
    CoRegFlux - Integration of transcriptional regulatory networks and gene expression.

    CoRegFlux integrates reverse engineered transcriptional regulatory networks and
    gene expression into metabolic models to improve phenotype prediction. It builds
    a linear regression estimator to predict target gene expression as a function of
    regulator co-expression, using influence scores as input.

    Author: Pauline Trébulle, Daniel Trejo-Banos, Mohamed Elati
    For more details: https://dx.doi.org/10.1186%2Fs12918-017-0507-0
    """

    def __init__(
        self, model: RegulatoryExtension, solver: Union[str, Solver] = None, build: bool = False, attach: bool = False
    ):
        """
        Initialize CoRegFlux with a RegulatoryExtension model.

        :param model: A RegulatoryExtension instance wrapping a simulator
        :param solver: A Solver instance or solver name. If None, a new solver is instantiated.
        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)

    def next_state(
        self,
        solver_kwargs: Dict = None,
        state: Dict[str, float] = None,
        metabolites: Dict[str, CoRegMetabolite] = None,
        biomass: CoRegBiomass = None,
        time_step: float = None,
        soft_plus: float = 0,
        tolerance: float = ModelConstants.TOLERANCE,
        scale: bool = False,
    ) -> CoRegResult:
        """
        Compute the next state of the system given the current state and time step.

        :param solver_kwargs: solver arguments
        :param state: current state of the system
        :param metabolites: metabolites constraints
        :param biomass: biomass constraints
        :param time_step: time step
        :param soft_plus: soft plus parameter
        :param tolerance: tolerance
        :param scale: whether to scale the metabolites
        :return: next state of the system
        """
        result = CoRegResult()

        # Get reaction constraints from simulator
        # yield_reactions() returns reaction IDs (strings), not objects
        constraints = {}
        for rxn_id in self.model.yield_reactions():
            rxn_data = self._get_reaction(rxn_id)
            constraints[rxn_id] = (
                rxn_data.get("lb", ModelConstants.REACTION_LOWER_BOUND),
                rxn_data.get("ub", ModelConstants.REACTION_UPPER_BOUND),
            )

        if metabolites:
            # Update coregflux constraints using metabolite concentrations
            constraints = metabolites_constraints(
                constraints=constraints, metabolites=metabolites, biomass=biomass, time_step=time_step
            )

        if state:
            # Update coregflux bounds using gene state
            constraints = gene_state_constraints(
                model=self.model,
                constraints=constraints,
                state=state,
                soft_plus=soft_plus,
                tolerance=tolerance,
                scale=scale,
            )

        # Retrieve the FBA simulation from the inferred constraints
        values, objective_value = _run_and_decode(self, additional_constraints=constraints, solver_kwargs=solver_kwargs)
        result.values = values
        result.objective_value = objective_value

        # Update the system state by solving an Euler step for metabolites and biomass
        next_biomass, next_metabolites = system_state_update(
            model=self.model,
            flux_state=values,
            metabolites=metabolites,
            biomass=biomass,
            time_step=time_step,
            biomass_fn=biomass_yield_to_rate,
        )

        result.metabolites = next_metabolites
        result.biomass = next_biomass
        result.constraints = constraints
        return result

    def _dynamic_optimize(
        self,
        to_solver: bool = False,
        solver_kwargs: Dict = None,
        initial_state: Sequence[Dict[str, float]] = None,
        metabolites: Dict[str, CoRegMetabolite] = None,
        biomass: CoRegBiomass = None,
        time_steps: Sequence[float] = None,
        soft_plus: float = 0,
        tolerance: float = ModelConstants.TOLERANCE,
        scale: bool = False,
    ) -> Union[DynamicSolution, Dict[float, Solution]]:
        """Dynamic optimization over multiple time steps."""
        solutions = []

        previous_time_step = 0
        for i_initial_state, time_step in zip(initial_state, time_steps):
            time_step_diff = time_step - previous_time_step

            next_state = self.next_state(
                solver_kwargs=solver_kwargs,
                state=i_initial_state,
                metabolites=metabolites,
                biomass=biomass,
                time_step=time_step_diff,
                soft_plus=soft_plus,
                tolerance=tolerance,
                scale=scale,
            )

            solution = result_to_solution(result=next_state, model=self.model, to_solver=to_solver)
            solutions.append(solution)

            metabolites = next_state.metabolites
            biomass = next_state.biomass

            previous_time_step = time_step

        # DynamicSolution expects positional args, not keyword 'solutions'
        return DynamicSolution(*solutions, time=time_steps)

    def optimize(
        self,
        initial_state: Union[Dict[str, float], Sequence[Dict[str, float]]] = None,
        metabolites: Dict[str, Union[float, CoRegMetabolite]] = None,
        biomass: Union[float, CoRegBiomass] = None,
        time_steps: Union[float, Sequence[float]] = None,
        soft_plus: float = 0,
        tolerance: float = ModelConstants.TOLERANCE,
        scale: bool = False,
        to_solver: bool = False,
        solver_kwargs: Dict = None,
    ) -> Union[Solution, DynamicSolution]:
        """
        Solve the CoRegFlux problem.

        :param initial_state: Initial gene state or sequence of states
        :param metabolites: Initial metabolite concentrations
        :param biomass: Initial biomass
        :param time_steps: Time step(s) for simulation
        :param soft_plus: Soft plus parameter
        :param tolerance: Tolerance for constraints
        :param scale: Whether to scale metabolites
        :param to_solver: Whether to return raw solver solution
        :param solver_kwargs: Additional solver arguments
        :return: Solution or DynamicSolution
        """
        # Build solver if out of sync
        if not self.synchronized:
            self.build()

        if not solver_kwargs:
            solver_kwargs = {}

        # Build metabolites and biomass objects
        if metabolites is None:
            metabolites = {}
        metabolites = build_metabolites(model=self.model, metabolites=metabolites)

        if biomass is None:
            # Calculate initial growth rate if not provided
            _, growth_rate = _run_and_decode(self, solver_kwargs=solver_kwargs)
            biomass = growth_rate
        biomass = build_biomass(model=self.model, biomass=biomass)

        # Handle single vs dynamic simulation
        if isinstance(initial_state, dict):
            # Single time step
            if time_steps is None:
                time_steps = 0.1

            result = self.next_state(
                solver_kwargs=solver_kwargs,
                state=initial_state,
                metabolites=metabolites,
                biomass=biomass,
                time_step=time_steps,
                soft_plus=soft_plus,
                tolerance=tolerance,
                scale=scale,
            )

            return result_to_solution(result=result, model=self.model, to_solver=to_solver)

        else:
            # Dynamic simulation with multiple time steps
            if time_steps is None:
                time_steps = np.linspace(0, 1, len(initial_state))

            return self._dynamic_optimize(
                to_solver=to_solver,
                solver_kwargs=solver_kwargs,
                initial_state=initial_state,
                metabolites=metabolites,
                biomass=biomass,
                time_steps=time_steps,
                soft_plus=soft_plus,
                tolerance=tolerance,
                scale=scale,
            )


# ----------------------------------------------------------------------------------------------------------------------
# Gene Expression Prediction
# ----------------------------------------------------------------------------------------------------------------------
def _get_target_regulators(gene: Union["Gene", "Target"] = None) -> List[str]:
    """
    Return the list of regulators of a target gene.

    :param gene: Target gene
    :return: List of regulators of the target gene
    """
    if gene is None:
        return []

    if gene.is_target():
        return [regulator.id for regulator in gene.yield_regulators()]

    return []


def _filter_influence_and_expression(
    interactions: Dict[str, List[str]], influence: pd.DataFrame, expression: pd.DataFrame, experiments: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Filter influence, expression and experiments matrices to keep only targets and their regulators.

    :param interactions: Dictionary with the interactions between targets and regulators
    :param influence: Influence matrix
    :param expression: Expression matrix
    :param experiments: Experiments matrix
    :return: Filtered influence matrix, filtered expression matrix, filtered experiments matrix
    """
    targets = pd.Index(set(interactions.keys()))
    regulators = pd.Index(set([regulator for regulators in interactions.values() for regulator in regulators]))

    # Filter the expression matrix for the target genes only
    expression = expression.loc[expression.index.intersection(targets)].copy()
    influence = influence.loc[influence.index.intersection(regulators)].copy()
    experiments = experiments.loc[experiments.index.intersection(regulators)].copy()
    return influence, expression, experiments


def _predict_experiment(
    interactions: Dict[str, List[str]], influence: pd.DataFrame, expression: pd.DataFrame, experiment: pd.Series
) -> pd.Series:
    """Predict gene expression for a single experiment using linear regression."""
    try:
        # noinspection PyPackageRequirements
        from sklearn.linear_model import LinearRegression
    except ImportError:
        raise ImportError(
            "The package sklearn is not installed. "
            "To compute the probability of target-regulator interactions, please install sklearn "
            "(pip install sklearn)."
        )

    predictions = {}

    for target, regulators in interactions.items():
        if not regulators:
            predictions[target] = np.nan
            continue

        if target not in expression.index:
            predictions[target] = np.nan
            continue

        if not set(regulators).issubset(influence.index):
            predictions[target] = np.nan
            continue

        if not set(regulators).issubset(experiment.index):
            predictions[target] = np.nan
            continue

        # Train linear regression model:
        # y = expression of the target gene for all samples
        # x1 = influence score of regulator 1 in training dataset
        # x2 = influence score of regulator 2 in training dataset
        # etc.
        x = influence.loc[regulators].transpose().to_numpy()
        y = expression.loc[target].to_numpy()

        regressor = LinearRegression()
        regressor.fit(x, y)

        # Predict the expression of the target gene for the experiment
        x_pred = experiment.loc[regulators].to_frame().transpose().to_numpy()
        predictions[target] = regressor.predict(x_pred)[0]

    return pd.Series(predictions)


def predict_gene_expression(
    model: RegulatoryExtension, influence: pd.DataFrame, expression: pd.DataFrame, experiments: pd.DataFrame
) -> pd.DataFrame:
    """
    Predict gene expression in experiments using co-expression of regulators.

    Adapted from CoRegFlux documentation:
    - A GEM model containing GPRs
    - A TRN network containing target genes and co-activators/co-repressors
    - GEM genes and TRN targets must match
    - An influence dataset containing influence scores of regulators (similar to correlation)
      Format: (rows: regulators, columns: samples)
      Influence scores calculated using CoRegNet algorithm
    - A gene expression dataset (training dataset)
      Format: (rows: genes, columns: samples)
    - An experiments dataset containing influence scores for test conditions
      Format: (rows: regulators, columns: experiments/conditions)

    :param model: A RegulatoryExtension instance
    :param influence: Influence scores of regulators in training dataset
    :param expression: Expression of genes in training dataset
    :param experiments: Influence scores of regulators in test dataset
    :return: Predicted expression of genes in test dataset
    """
    # Filter only gene expression and influences data of metabolic genes in the model
    # yield_targets() returns (target_id, target_object) tuples
    # Handle both legacy GERM models (yield Target) and RegulatoryExtension (yield tuples)
    targets_gen = model.yield_targets()
    first_target = next(targets_gen, None)

    if first_target is None:
        interactions = {}
    elif isinstance(first_target, tuple):
        # RegulatoryExtension: yields (id, Target) tuples
        interactions = {first_target[1].id: _get_target_regulators(first_target[1])}
        interactions.update({target.id: _get_target_regulators(target) for _, target in targets_gen})
    else:
        # Legacy GERM: yields Target objects directly
        interactions = {first_target.id: _get_target_regulators(first_target)}
        interactions.update({target.id: _get_target_regulators(target) for target in targets_gen})
    influence, expression, experiments = _filter_influence_and_expression(
        interactions=interactions, influence=influence, expression=expression, experiments=experiments
    )

    predictions = []
    for column in experiments.columns:
        experiment_prediction = _predict_experiment(
            interactions=interactions, influence=influence, expression=expression, experiment=experiments[column]
        )
        predictions.append(experiment_prediction)

    predictions = pd.concat(predictions, axis=1)
    predictions.columns = experiments.columns
    return predictions.dropna()
