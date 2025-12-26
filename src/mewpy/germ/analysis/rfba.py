"""
Regulatory Flux Balance Analysis (RFBA) - Clean Implementation

This module implements RFBA using RegulatoryExtension only.
No backwards compatibility with legacy GERM models.
"""
from typing import Union, Dict, Tuple
from warnings import warn

from mewpy.germ.analysis.fba import _RegulatoryAnalysisBase
from mewpy.solvers import Solution
from mewpy.solvers.solver import Solver
from mewpy.solvers import solver_instance
from mewpy.germ.models.regulatory_extension import RegulatoryExtension
from mewpy.util.constants import ModelConstants
from mewpy.germ.solution import DynamicSolution


class RFBA(_RegulatoryAnalysisBase):
    """
    Regulatory Flux Balance Analysis (RFBA) using RegulatoryExtension.

    RFBA integrates transcriptional regulatory networks with metabolic models
    to predict cellular behavior under regulatory constraints.

    This implementation:
    - Works exclusively with RegulatoryExtension instances
    - Delegates all metabolic operations to external simulators (cobrapy/reframed)
    - Stores only regulatory network information
    - Falls back to FBA if no regulatory network present

    For more details: Covert et al. 2004, https://doi.org/10.1038/nature02456
    """

    def __init__(self,
                 model: RegulatoryExtension,
                 solver: Union[str, Solver, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Initialize RFBA with a RegulatoryExtension model.

        :param model: A RegulatoryExtension instance wrapping a simulator
        :param solver: Solver instance or name. If None, uses default solver.
        :param build: If True, builds the problem immediately. Default: False
        :param attach: If True, attaches problem to model. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)
        self.method = "RFBA"

    def build(self):
        """
        Build the RFBA linear problem.

        Creates the solver instance and sets up the objective function.
        For models without regulatory networks, this is equivalent to FBA.
        """
        # Get simulator from RegulatoryExtension
        simulator = self.model.simulator

        # Create solver from simulator
        self._solver = solver_instance(simulator)

        # Set up objective (always string keys from simulator)
        self._linear_objective = dict(self.model.objective)
        self._minimize = False

        # Mark as synchronized
        self._synchronized = True

        return self

    def decode_regulatory_state(self, state: Dict[str, float]) -> Dict[str, float]:
        """
        Solve the boolean regulatory network for a given state.

        Evaluates all regulatory interactions and returns the resulting
        target gene states.

        :param state: Dictionary mapping regulator IDs to their states (0.0 or 1.0)
        :return: Dictionary mapping target gene IDs to their resulting states
        """
        # If no regulatory network, return empty dict
        if not self.model.has_regulatory_network():
            return {}

        # Evaluate all interactions synchronously
        result = {}
        for interaction in self.model.yield_interactions():
            # An interaction can have multiple regulatory events
            # (e.g., coefficient 1.0 if condition A, 0.0 otherwise)
            for coefficient, event in interaction.regulatory_events.items():
                if event.is_none:
                    continue

                # Evaluate the regulatory event with current state
                eval_result = event.evaluate(values=state)
                if eval_result:
                    result[interaction.target.id] = coefficient
                else:
                    result[interaction.target.id] = 0.0

        return result

    def decode_metabolic_state(self, state: Dict[str, float]) -> Dict[str, float]:
        """
        Decode metabolic state from regulatory state.

        For most models, this is identical to decode_regulatory_state.

        :param state: Dictionary mapping regulator IDs to their states
        :return: Dictionary mapping metabolic gene IDs to their states
        """
        return self.decode_regulatory_state(state)

    def decode_constraints(self, state: Dict[str, float]) -> Dict[str, Tuple[float, float]]:
        """
        Decode metabolic constraints from gene states.

        Evaluates GPR rules for all reactions and returns constraints
        for reactions whose GPRs evaluate to False (genes knocked out).

        :param state: Dictionary mapping gene IDs to their states (0.0 or 1.0)
        :return: Dictionary mapping reaction IDs to bounds (0.0, 0.0) for inactive reactions
        """
        constraints = {}

        # Evaluate GPRs for all reactions
        for rxn_id in self.model.reactions:
            # Get cached parsed GPR expression from RegulatoryExtension
            gpr = self.model.get_parsed_gpr(rxn_id)

            if gpr.is_none:
                continue

            # Evaluate GPR with current gene states
            is_active = gpr.evaluate(values=state)

            if not is_active:
                # Reaction is knocked out - set bounds to zero
                constraints[rxn_id] = (0.0, 0.0)

        return constraints

    def initial_state(self, initial_state: Dict[str, float] = None) -> Dict[str, float]:
        """
        Get initial regulatory state for RFBA simulation.

        :param initial_state: Optional user-provided initial state
        :return: Complete initial state dictionary
        """
        if initial_state is None:
            initial_state = {}

        # Initialize all regulators to active (1.0) by default
        state = {}

        if self.model.has_regulatory_network():
            # Get all regulators
            for reg_id, regulator in self.model.yield_regulators():
                # Use user-provided state if available, otherwise default to active
                state[reg_id] = initial_state.get(reg_id, 1.0)

        # Override with any user-provided states
        state.update(initial_state)

        return state

    def optimize(self,
                 initial_state: Dict[str, float] = None,
                 dynamic: bool = False,
                 to_solver: bool = False,
                 solver_kwargs: Dict = None) -> Union[Solution, DynamicSolution]:
        """
        Solve the RFBA problem.

        :param initial_state: Initial regulatory state. If None, all regulators start active.
        :param dynamic: If True, performs dynamic RFBA (iterative). Default: False
        :param to_solver: If True, returns raw solver solution. Default: False
        :param solver_kwargs: Additional arguments for solver
        :return: Solution object (or DynamicSolution for dynamic=True)
        """
        # Build solver if not synchronized
        if not self.synchronized:
            self.build()

        if solver_kwargs is None:
            solver_kwargs = {}

        # Get initial state
        state = self.initial_state(initial_state)

        if not dynamic:
            # Steady-state RFBA
            return self._optimize_steady_state(state, to_solver, solver_kwargs)
        else:
            # Dynamic RFBA (iterative until convergence)
            return self._optimize_dynamic(state, to_solver, solver_kwargs)

    def _optimize_steady_state(self,
                                state: Dict[str, float],
                                to_solver: bool,
                                solver_kwargs: Dict) -> Solution:
        """
        Perform steady-state RFBA simulation.

        :param state: Initial regulatory state
        :param to_solver: Whether to return raw solver solution
        :param solver_kwargs: Solver arguments
        :return: Solution object
        """
        # Decode regulatory state to metabolic state
        metabolic_state = self.decode_metabolic_state(state)

        # Get constraints from metabolic state
        constraints = self.decode_constraints(metabolic_state)

        # Merge with user-provided constraints
        if 'constraints' in solver_kwargs:
            constraints.update(solver_kwargs['constraints'])

        # Solve
        solution = self.solver.solve(
            linear=self._linear_objective,
            minimize=self._minimize,
            constraints=constraints,
            get_values=True,
            **solver_kwargs
        )

        if to_solver:
            return solution

        return Solution.from_solver(
            method=self.method,
            solution=solution,
            model=self.model,
            minimize=self._minimize
        )

    def _optimize_dynamic(self,
                          state: Dict[str, float],
                          to_solver: bool,
                          solver_kwargs: Dict) -> DynamicSolution:
        """
        Perform dynamic RFBA simulation (iterative until convergence).

        Dynamic RFBA iteratively:
        1. Solves FBA with current regulatory state
        2. Updates regulatory state based on solution
        3. Repeats until steady state (no changes) or max iterations

        :param state: Initial regulatory state
        :param to_solver: Whether to return raw solver solutions
        :param solver_kwargs: Solver arguments
        :return: DynamicSolution containing all iterations
        """
        solutions = []
        max_iterations = 100  # Safety limit

        for iteration in range(max_iterations):
            # Solve with current state
            solution = self._optimize_steady_state(state, to_solver=True, solver_kwargs=solver_kwargs)
            solutions.append(solution)

            # Check if solution is optimal
            if solution.status.name != 'OPTIMAL':
                warn(f"Non-optimal solution at iteration {iteration}")
                break

            # Get new regulatory state from solution
            # Update state based on reaction fluxes and metabolite concentrations
            new_state = self._update_state_from_solution(state, solution)

            # Check for convergence (state hasn't changed)
            if self._states_equal(state, new_state):
                break

            # Update state for next iteration
            state = new_state

        return DynamicSolution(solutions=solutions, method=self.method)

    def _update_state_from_solution(self,
                                     current_state: Dict[str, float],
                                     solution) -> Dict[str, float]:
        """
        Update regulatory state based on FBA solution.

        This examines reaction fluxes and metabolite concentrations
        to determine new regulator states.

        :param current_state: Current regulatory state
        :param solution: FBA solution
        :return: Updated regulatory state
        """
        new_state = current_state.copy()

        # If no regulatory network, return unchanged
        if not self.model.has_regulatory_network():
            return new_state

        # Update regulator states based on solution
        # Note: This is model-specific logic that may need customization
        # For now, we use a simple approach based on flux values

        for reg_id in new_state.keys():
            # Check if regulator is a reaction or metabolite
            if reg_id in self.model.reactions:
                # Regulator is a reaction - update based on flux
                flux = solution.values.get(reg_id, 0.0)
                new_state[reg_id] = 1.0 if abs(flux) > ModelConstants.TOLERANCE else 0.0

            elif reg_id in self.model.metabolites:
                # Regulator is a metabolite - could use concentration if available
                # For now, keep current state
                pass

        return new_state

    def _states_equal(self, state1: Dict[str, float], state2: Dict[str, float]) -> bool:
        """
        Check if two regulatory states are equal.

        :param state1: First state
        :param state2: Second state
        :return: True if states are equal, False otherwise
        """
        if set(state1.keys()) != set(state2.keys()):
            return False

        for key in state1:
            if abs(state1[key] - state2[key]) > ModelConstants.TOLERANCE:
                return False

        return True
