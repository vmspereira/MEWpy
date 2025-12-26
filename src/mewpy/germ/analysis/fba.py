"""
Internal base class for regulatory analysis methods.

This module provides a minimal base class for regulatory analysis methods
(RFBA, SRFBA, PROM, CoRegFlux). It is NOT intended for direct use by end users.

For pure FBA without regulatory networks, use the simulator directly:
    solution = model.simulator.optimize()

Or use COBRApy/reframed FBA implementations directly:
    solution = cobra_model.optimize()  # COBRApy
    solution = reframed_model.optimize()  # reframed
"""
from typing import Union, Dict

from mewpy.germ.models.regulatory_extension import RegulatoryExtension
from mewpy.solvers.solution import Solution
from mewpy.solvers.solver import Solver
from mewpy.solvers import solver_instance


class _RegulatoryAnalysisBase:
    """
    Internal base class for regulatory analysis methods.

    This class provides common functionality for regulatory analysis methods:
    - Solver management
    - Build pattern (create solver from simulator)
    - Basic optimization interface

    NOT intended for direct use. For pure FBA, use:
    - model.simulator.optimize() for RegulatoryExtension
    - cobra_model.optimize() for COBRApy
    - reframed_model.optimize() for reframed

    Subclasses: RFBA, SRFBA, PROM, CoRegFlux
    """

    def __init__(self,
                 model: RegulatoryExtension,
                 solver: Union[str, Solver, None] = None,
                 build: bool = False,
                 attach: bool = False):
        """
        Initialize regulatory analysis base.

        :param model: A RegulatoryExtension instance wrapping a simulator
        :param solver: A Solver instance or solver name. If None, a new solver is instantiated.
        :param build: Whether to build the problem upon instantiation. Default: False
        :param attach: Whether to attach the problem to the model upon instantiation. Default: False
        """
        self.model = model
        self.solver_name = solver
        self._solver = None
        self._linear_objective = None
        self._minimize = False
        self._synchronized = False
        self.method = "FBA"  # Subclasses should override

        if build:
            self.build()

        if attach:
            # TODO: Implement attach functionality if needed
            pass

    # Backwards compatibility helpers (work with both RegulatoryExtension and legacy models)
    def _has_regulatory_network(self) -> bool:
        """Check if model has a regulatory network (works with both model types)."""
        if hasattr(self.model, 'has_regulatory_network'):
            return self.model.has_regulatory_network()
        # Legacy model - check for interactions
        return hasattr(self.model, 'interactions') and len(self.model.interactions) > 0

    def _get_interactions(self):
        """Get interactions (works with both model types)."""
        if hasattr(self.model, 'yield_interactions'):
            return self.model.yield_interactions()
        # Legacy model
        return self.model.interactions.values() if hasattr(self.model, 'interactions') else []

    def _get_regulators(self):
        """Get regulators (works with both model types)."""
        if hasattr(self.model, 'yield_regulators'):
            # RegulatoryExtension yields (id, regulator) tuples
            # Legacy models yield single Regulator objects
            # We need to normalize to always return (id, regulator) tuples
            for item in self.model.yield_regulators():
                if isinstance(item, tuple) and len(item) == 2:
                    # Already a tuple (RegulatoryExtension)
                    yield item
                else:
                    # Single Regulator object (legacy model) - wrap in tuple with ID
                    yield (item.id, item)
        # Legacy model without yield_regulators (shouldn't happen but handle it)
        elif hasattr(self.model, 'regulators'):
            regulators = self.model.regulators
            # Check if it's a dict
            if isinstance(regulators, dict):
                for reg_id, regulator in regulators.items():
                    yield (reg_id, regulator)
        # No regulators available
        return

    def _get_gpr(self, rxn_id):
        """Get GPR expression (works with both model types)."""
        if hasattr(self.model, 'get_parsed_gpr'):
            return self.model.get_parsed_gpr(rxn_id)
        # Legacy model - get reaction and parse GPR
        from mewpy.germ.algebra import parse_expression, Expression, Symbol
        if hasattr(self.model, 'reactions') and rxn_id in self.model.reactions:
            rxn = self.model.reactions[rxn_id]
            if hasattr(rxn, 'gpr'):
                return rxn.gpr
        # No GPR available
        return Expression(Symbol('true'), {})

    def _get_reaction(self, rxn_id):
        """Get reaction data (works with both model types)."""
        if hasattr(self.model, 'get_reaction'):
            # RegulatoryExtension - returns dict
            return self.model.get_reaction(rxn_id)
        # Legacy model - reactions is a dict of Reaction objects
        if hasattr(self.model, 'reactions') and rxn_id in self.model.reactions:
            rxn = self.model.reactions[rxn_id]
            # Convert Reaction object to dict format
            return {
                'id': rxn.id,
                'lb': rxn.lower_bound if hasattr(rxn, 'lower_bound') else rxn.bounds[0] if hasattr(rxn, 'bounds') else -1000,
                'ub': rxn.upper_bound if hasattr(rxn, 'upper_bound') else rxn.bounds[1] if hasattr(rxn, 'bounds') else 1000,
                'gpr': str(rxn.gpr) if hasattr(rxn, 'gpr') else ''
            }
        return {'id': rxn_id, 'lb': -1000, 'ub': 1000, 'gpr': ''}

    def build(self):
        """
        Build the optimization problem.

        Creates the solver instance from the simulator and sets up the objective function.
        Subclasses should call super().build() and then add their specific constraints.
        """
        # Get simulator - support both RegulatoryExtension and legacy models
        if hasattr(self.model, 'simulator'):
            # RegulatoryExtension
            simulator = self.model.simulator
        else:
            # Legacy model or direct simulator
            # Try to get a simulator from it
            from mewpy.simulation import get_simulator
            try:
                simulator = get_simulator(self.model)
            except:
                # If that fails, assume it's already a simulator
                simulator = self.model

        # Create solver directly from simulator
        self._solver = solver_instance(simulator)

        # Set up objective
        if hasattr(self.model, 'objective'):
            objective = self.model.objective
            # Handle different objective formats
            if isinstance(objective, dict):
                # Already a dict - check if keys are objects or strings
                first_key = next(iter(objective.keys())) if objective else None
                if first_key and hasattr(first_key, 'id'):
                    # Keys are objects (legacy), convert to string keys
                    self._linear_objective = {var.id: value for var, value in objective.items()}
                else:
                    # Keys are already strings
                    self._linear_objective = dict(objective)
            else:
                # Some other format
                self._linear_objective = dict(objective)
        else:
            self._linear_objective = {}

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
        Optimize the problem.

        This basic implementation is used by some subclasses (e.g., SRFBA).
        Other subclasses (e.g., RFBA) override this completely.

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

        # Solve using simulator
        solution = self.solver.solve(
            linear=self._linear_objective,
            minimize=self._minimize,
            **solver_kwargs_copy
        )

        # Set the method attribute for compatibility
        solution._method = self.method
        solution._model = self.model

        return solution


# Alias for backwards compatibility during transition
# Users should not use this directly - use simulator.optimize() instead
FBA = _RegulatoryAnalysisBase
