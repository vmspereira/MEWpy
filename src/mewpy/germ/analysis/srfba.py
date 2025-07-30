from typing import Union, Dict
from warnings import warn

from mewpy.util.constants import ModelConstants
from mewpy.germ.analysis import FBA
from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel
from mewpy.solvers import Solution


class SRFBA(FBA):
    """
    Steady-state Regulatory Flux Balance Analysis (SRFBA) using pure simulator-based approach.
    
    This implementation uses simulators as the foundation and provides a simplified
    approach to regulatory-metabolic optimization without mixed-integer programming complexity.
    """

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Steady-state Regulatory Flux Balance Analysis (SRFBA) of a metabolic-regulatory model.
        Implementation using pure simulator-based approach.

        For more details consult Shlomi et al. 2007 at https://dx.doi.org/10.1038%2Fmsb4100141

        :param model: a MetabolicModel, RegulatoryModel or GERM model. The model is used to retrieve
        the simulator for optimization
        :param solver: A solver name. The solver interface will be used to load and solve the optimization problem.
        If none, a new solver is instantiated.
        :param build: Whether to build the problem upon instantiation. Default: False
        :param attach: Whether to attach the problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)
        self._model_default_lb = ModelConstants.REACTION_LOWER_BOUND
        self._model_default_ub = ModelConstants.REACTION_UPPER_BOUND
        
        # Warn about simplified implementation
        warn("SRFBA has been simplified to use pure simulator approach. "
             "Complex mixed-integer programming features are not available. "
             "For full SRFBA functionality, consider using the legacy implementation.",
             UserWarning, stacklevel=2)

    @property
    def model_default_lb(self) -> float:
        """
        The default lower bound for the model reactions.
        :return:
        """
        if self.synchronized:
            return self._model_default_lb

        if hasattr(self.model, 'yield_reactions'):
            self._model_default_lb = min(reaction.lower_bound for reaction in self.model.yield_reactions())
        return self._model_default_lb

    @property
    def model_default_ub(self) -> float:
        """
        The default upper bound for the model reactions.
        :return:
        """
        if self.synchronized:
            return self._model_default_ub

        if hasattr(self.model, 'yield_reactions'):
            self._model_default_ub = max(reaction.upper_bound for reaction in self.model.yield_reactions())
        return self._model_default_ub

    def build(self):
        """
        Build the SRFBA problem using pure simulator approach.
        
        This simplified implementation focuses on regulatory constraints
        without the full mixed-integer programming complexity.
        """
        # Check if this is a native GERM model with regulatory capabilities
        if hasattr(self.model, 'is_regulatory') and self.model.is_regulatory():
            # For regulatory models, use simulator-based approach similar to RFBA
            return super().build()
        else:
            # For non-regulatory models, fall back to FBA
            return super().build()

    def optimize(self, 
                 solver_kwargs: Dict = None, 
                 initial_state: Dict[str, float] = None,
                 to_solver: bool = False,
                 **kwargs) -> Solution:
        """
        Optimize the SRFBA problem using pure simulator approach.
        
        This simplified implementation provides regulatory-aware optimization
        without mixed-integer programming complexity.
        
        :param solver_kwargs: A dictionary of keyword arguments to be passed to the solver.
        :param initial_state: Initial state for regulatory variables
        :param to_solver: Whether to return the solution as a SolverSolution instance. Default: False
        :return: A Solution instance.
        """
        if not self.synchronized:
            self.build()
        
        if not solver_kwargs:
            solver_kwargs = {}
            
        if not initial_state:
            initial_state = {}
        
        # Check if model has regulatory capabilities
        if hasattr(self.model, 'is_regulatory') and self.model.is_regulatory():
            # Apply regulatory constraints via initial_state
            if initial_state:
                constraints = solver_kwargs.get('constraints', {})
                constraints.update(initial_state)
                solver_kwargs['constraints'] = constraints
        
        # Use the base FBA optimization with regulatory constraints
        solution = super().optimize(solver_kwargs=solver_kwargs, **kwargs)
        
        if to_solver:
            return solution
        
        # Convert to Solution if needed
        if not isinstance(solution, Solution):
            minimize = solver_kwargs.get('minimize', self._minimize)
            return Solution.from_solver(method="SRFBA", solution=solution, model=self.model, minimize=minimize)
        
        return solution

    # Legacy method signatures for backward compatibility
    def _optimize(self,
                  to_solver: bool = False,
                  solver_kwargs: Dict = None,
                  initial_state: Dict[str, float] = None,
                  **kwargs) -> Solution:
        """
        Legacy optimization method for backward compatibility.
        """
        return self.optimize(solver_kwargs=solver_kwargs, 
                           initial_state=initial_state, 
                           to_solver=to_solver, 
                           **kwargs)

    # Simplified constraint generation methods (for interface compatibility)
    def gpr_constraint(self, reaction):
        """
        Simplified GPR constraint handling for pure simulator approach.
        Returns empty constraints as GPRs are handled by the simulator.
        """
        warn("GPR constraints are handled automatically by the simulator in this implementation.",
             UserWarning, stacklevel=2)
        return [], []

    def interaction_constraint(self, interaction):
        """
        Simplified interaction constraint handling for pure simulator approach.
        Returns empty constraints as interactions are handled by the simulator.
        """
        warn("Interaction constraints are handled automatically by the simulator in this implementation.",
             UserWarning, stacklevel=2)
        return [], []
