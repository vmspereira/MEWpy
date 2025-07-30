from typing import Union, Dict

from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel
from mewpy.solvers.solution import Solution
from mewpy.solvers.solver import Solver
from mewpy.solvers import solver_instance


class FBA:
    """
    Flux Balance Analysis (FBA) of a metabolic model using pure simulator-based approach.
    
    This implementation uses simulators as the foundation for all models, providing
    a clean, unified architecture for metabolic analysis.
    """

    def __init__(self,
                 model: Union[Model, MetabolicModel, RegulatoryModel],
                 solver: Union[str, Solver, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Flux Balance Analysis (FBA) of a metabolic model. Pure simulator-based implementation.

        For more details consult: https://dx.doi.org/10.1038%2Fnbt.1614

        :param model: a MetabolicModel, RegulatoryModel or GERM model. The model is used to retrieve
        the simulator for optimization
        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve the optimization problem.
        If none, a new solver is instantiated.
        :param build: Whether to build the problem upon instantiation. Default: False
        :param attach: Whether to attach the problem to the model upon instantiation. Default: False
        """
        self.model = model
        self.solver_name = solver
        self._solver = None
        self._linear_objective = None
        self._minimize = False
        self._synchronized = False
        self.method = "FBA"  # Method name for solution creation
        
        if build:
            self.build()
            
        if attach:
            # TODO: Implement attach functionality if needed
            pass

    def _get_simulator(self):
        """Get the simulator from the model."""
        if hasattr(self.model, '_simulator'):
            return self.model._simulator
        elif hasattr(self.model, 'simulator'):
            return self.model.simulator
        else:
            # For native GERM models, we need to convert them to simulators
            from mewpy.simulation import get_simulator
            return get_simulator(self.model)

    def build(self):
        """
        Build the FBA problem using pure simulator approach.
        """
        # Get simulator from any model type
        simulator = self._get_simulator()
        
        # Create solver directly from simulator
        self._solver = solver_instance(simulator)
        
        # Set up the objective based on the model's objective
        self._linear_objective = {var.id: value for var, value in self.model.objective.items()}
        self._minimize = False
        
        # Mark as synchronized
        self._synchronized = True
        
        # Return self for chaining
        return self

    @property
    def synchronized(self):
        """Whether the solver is synchronized with the model."""
        return self._synchronized

    @property
    def solver(self):
        """Get the solver instance."""
        if self._solver is None:
            self.build()
        return self._solver

    def optimize(self, solver_kwargs: Dict = None, **kwargs) -> Solution:
        """
        Optimize the FBA problem using pure simulator approach.
        
        :param solver_kwargs: A dictionary of keyword arguments to be passed to the solver.
        :return: A Solution instance.
        """
        if not self.synchronized:
            self.build()
        
        if not solver_kwargs:
            solver_kwargs = {}
        
        # Make a copy to avoid modifying the original
        solver_kwargs_copy = solver_kwargs.copy()
        
        # Remove conflicting arguments that we set explicitly
        solver_kwargs_copy.pop('linear', None)
        solver_kwargs_copy.pop('minimize', None)
        
        # Pure simulator approach - clean and simple
        solution = self.solver.solve(
            linear=self._linear_objective,
            minimize=self._minimize,
            **solver_kwargs_copy
        )
        
        # Set the method attribute for compatibility
        solution._method = self.method
        solution._model = self.model
        
        return solution
