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

from typing import Dict, Union

from mewpy.germ.models.regulatory_extension import RegulatoryExtension
from mewpy.solvers import solver_instance
from mewpy.solvers.solution import Solution
from mewpy.solvers.solver import Solver


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

    def __init__(
        self,
        model: RegulatoryExtension,
        solver: Union[str, Solver, None] = None,
        build: bool = False,
        attach: bool = False,
    ):
        """
        Initialize regulatory analysis base.

        :param model: A RegulatoryExtension instance wrapping a simulator
        :param solver: A Solver instance or solver name. If None, a new solver is instantiated.
        :param build: Whether to build the problem upon instantiation. Default: False
        :param attach: Whether to attach the problem to the model upon instantiation. Default: False
                       **Note:** This parameter is kept for backwards compatibility but is not used.
                       In the new architecture, analysis methods do not attach to models via observer pattern.
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

        # attach parameter is kept for backwards compatibility but is unused
        # In the new architecture with RegulatoryExtension, analysis methods do not
        # attach to models via observer pattern - they access data on-demand

    def __str__(self):
        """Simple string representation."""
        model_id = getattr(self.model, "id", str(self.model))
        return f"{self.method} for {model_id}"

    def __repr__(self):
        """
        Returns a formatted table representation of the analysis method.
        Displays method, model, objective, solver, and sync status.
        """
        # Get solver name
        if self._solver:
            solver_name = self._solver.__class__.__name__
        else:
            solver_name = "None"

        # Get model name/ID
        if hasattr(self.model, "id"):
            model_name = str(self.model.id)
        elif hasattr(self.model, "simulator") and hasattr(self.model.simulator, "model"):
            # RegulatoryExtension - get underlying model ID
            model_name = str(getattr(self.model.simulator.model, "id", "Unknown"))
        else:
            model_name = str(self.model)

        # Get model type info
        model_types = []
        if hasattr(self.model, "types"):
            model_types = list(self.model.types) if self.model.types else []
        elif hasattr(self.model, "is_metabolic") and self.model.is_metabolic():
            model_types.append("metabolic")
        if self._has_regulatory_network():
            if "regulatory" not in model_types:
                model_types.append("regulatory")
        model_type_str = ", ".join(model_types) if model_types else "metabolic"

        # Format objective
        if self._linear_objective:
            if len(self._linear_objective) == 1:
                key, val = next(iter(self._linear_objective.items()))
                objective_str = f"{key}: {val}"
            else:
                objective_str = f"{len(self._linear_objective)} objectives"
        else:
            objective_str = "None"

        # Get variables and constraints count if available
        vars_count = "N/A"
        constraints_count = "N/A"
        if self._solver:
            try:
                if hasattr(self._solver, "problem"):
                    problem = self._solver.problem
                    # Try SCIP methods first
                    if hasattr(problem, "getNVars"):
                        vars_count = problem.getNVars()
                    elif hasattr(problem, "getVars"):
                        vars_count = len(problem.getVars())
                    elif hasattr(problem, "variables"):
                        vars_count = len(problem.variables)

                    if hasattr(problem, "getNConss"):
                        constraints_count = problem.getNConss()
                    elif hasattr(problem, "getConss"):
                        constraints_count = len(problem.getConss())
                    elif hasattr(problem, "constraints"):
                        constraints_count = len(problem.constraints)
            except:
                pass

        # Build table
        lines = []
        lines.append("=" * 60)
        lines.append(f"{self.method}")
        lines.append("=" * 60)
        lines.append(f"{'Model:':<20} {model_name}")
        lines.append(f"{'Type:':<20} {model_type_str}")
        lines.append(f"{'Variables:':<20} {vars_count}")
        lines.append(f"{'Constraints:':<20} {constraints_count}")
        lines.append(f"{'Objective:':<20} {objective_str}")
        lines.append(f"{'Solver:':<20} {solver_name}")
        lines.append(f"{'Synchronized:':<20} {self.synchronized}")
        lines.append("=" * 60)

        return "\n".join(lines)

    def _repr_html_(self):
        """
        Returns an HTML table representation for Jupyter notebooks.
        """
        # Get solver name
        if self._solver:
            solver_name = self._solver.__class__.__name__
        else:
            solver_name = "None"

        # Get model name
        if hasattr(self.model, "id"):
            model_name = str(self.model.id)
        elif hasattr(self.model, "simulator") and hasattr(self.model.simulator, "model"):
            model_name = str(getattr(self.model.simulator.model, "id", "Unknown"))
        else:
            model_name = str(self.model)

        # Get model type
        model_types = []
        if hasattr(self.model, "types"):
            model_types = list(self.model.types) if self.model.types else []
        elif hasattr(self.model, "is_metabolic") and self.model.is_metabolic():
            model_types.append("metabolic")
        if self._has_regulatory_network():
            if "regulatory" not in model_types:
                model_types.append("regulatory")
        model_type_str = ", ".join(model_types) if model_types else "metabolic"

        # Format objective
        if self._linear_objective:
            if len(self._linear_objective) == 1:
                key, val = next(iter(self._linear_objective.items()))
                objective_str = f"{key}: {val}"
            else:
                objective_str = f"{len(self._linear_objective)} objectives"
        else:
            objective_str = "None"

        # Get variables and constraints count if available
        vars_count = "N/A"
        constraints_count = "N/A"
        if self._solver:
            try:
                if hasattr(self._solver, "problem"):
                    problem = self._solver.problem
                    # Try SCIP methods first
                    if hasattr(problem, "getNVars"):
                        vars_count = problem.getNVars()
                    elif hasattr(problem, "getVars"):
                        vars_count = len(problem.getVars())
                    elif hasattr(problem, "variables"):
                        vars_count = len(problem.variables)

                    if hasattr(problem, "getNConss"):
                        constraints_count = problem.getNConss()
                    elif hasattr(problem, "getConss"):
                        constraints_count = len(problem.getConss())
                    elif hasattr(problem, "constraints"):
                        constraints_count = len(problem.constraints)
            except:
                pass

        return f"""
        <table>
            <tr>
                <td><strong>Method</strong></td>
                <td>{self.method}</td>
            </tr>
            <tr>
                <td><strong>Model</strong></td>
                <td>{model_name}</td>
            </tr>
            <tr>
                <td><strong>Type</strong></td>
                <td>{model_type_str}</td>
            </tr>
            <tr>
                <td><strong>Variables</strong></td>
                <td>{vars_count}</td>
            </tr>
            <tr>
                <td><strong>Constraints</strong></td>
                <td>{constraints_count}</td>
            </tr>
            <tr>
                <td><strong>Objective</strong></td>
                <td>{objective_str}</td>
            </tr>
            <tr>
                <td><strong>Solver</strong></td>
                <td>{solver_name}</td>
            </tr>
            <tr>
                <td><strong>Synchronized</strong></td>
                <td>{self.synchronized}</td>
            </tr>
        </table>
        """

    # Backwards compatibility helpers (work with both RegulatoryExtension and legacy models)
    def _has_regulatory_network(self) -> bool:
        """Check if model has a regulatory network (works with both model types)."""
        if hasattr(self.model, "has_regulatory_network"):
            return self.model.has_regulatory_network()
        # Legacy model - check for interactions
        return hasattr(self.model, "interactions") and len(self.model.interactions) > 0

    def _get_interactions(self):
        """
        Get interactions (works with both model types).

        Yields just the Interaction objects, unpacking tuples from RegulatoryExtension.
        """
        if hasattr(self.model, "yield_interactions"):
            # RegulatoryExtension yields (id, interaction) tuples - unpack to get just interaction
            for item in self.model.yield_interactions():
                if isinstance(item, tuple) and len(item) == 2:
                    # New format: (id, interaction)
                    yield item[1]
                else:
                    # Legacy format: just interaction (shouldn't happen with RegulatoryExtension)
                    yield item
        # Legacy model - dict.values() yields just interaction objects
        elif hasattr(self.model, "interactions"):
            for interaction in self.model.interactions.values():
                yield interaction

    def _get_regulators(self):
        """Get regulators (works with both model types)."""
        if hasattr(self.model, "yield_regulators"):
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
        elif hasattr(self.model, "regulators"):
            regulators = self.model.regulators
            # Check if it's a dict
            if isinstance(regulators, dict):
                for reg_id, regulator in regulators.items():
                    yield (reg_id, regulator)
        # No regulators available
        return

    def _get_gpr(self, rxn_id):
        """Get GPR expression (works with both model types)."""
        if hasattr(self.model, "get_parsed_gpr"):
            return self.model.get_parsed_gpr(rxn_id)
        # Legacy model - get reaction and parse GPR
        from mewpy.germ.algebra import Expression, Symbol, parse_expression

        if hasattr(self.model, "reactions") and rxn_id in self.model.reactions:
            rxn = self.model.reactions[rxn_id]
            if hasattr(rxn, "gpr"):
                return rxn.gpr
        # No GPR available
        return Expression(Symbol("true"), {})

    def _get_reaction(self, rxn_id):
        """Get reaction data (works with both model types)."""
        if hasattr(self.model, "get_reaction"):
            # RegulatoryExtension - returns dict
            return self.model.get_reaction(rxn_id)
        # Legacy model - reactions is a dict of Reaction objects
        if hasattr(self.model, "reactions") and rxn_id in self.model.reactions:
            rxn = self.model.reactions[rxn_id]
            # Convert Reaction object to dict format
            return {
                "id": rxn.id,
                "lb": (
                    rxn.lower_bound
                    if hasattr(rxn, "lower_bound")
                    else rxn.bounds[0] if hasattr(rxn, "bounds") else -1000
                ),
                "ub": (
                    rxn.upper_bound
                    if hasattr(rxn, "upper_bound")
                    else rxn.bounds[1] if hasattr(rxn, "bounds") else 1000
                ),
                "gpr": str(rxn.gpr) if hasattr(rxn, "gpr") else "",
            }
        return {"id": rxn_id, "lb": -1000, "ub": 1000, "gpr": ""}

    def build(self):
        """
        Build the optimization problem.

        Creates the solver instance from the simulator and sets up the objective function.
        Subclasses should call super().build() and then add their specific constraints.
        """
        # Get simulator - support both RegulatoryExtension and legacy models
        if hasattr(self.model, "simulator"):
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
        if hasattr(self.model, "objective"):
            objective = self.model.objective
            # Handle different objective formats
            if isinstance(objective, dict):
                # Already a dict - check if keys are objects or strings
                first_key = next(iter(objective.keys())) if objective else None
                if first_key and hasattr(first_key, "id"):
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
        solver_kwargs_copy.pop("linear", None)
        solver_kwargs_copy.pop("minimize", None)

        # Solve using simulator
        solution = self.solver.solve(linear=self._linear_objective, minimize=self._minimize, **solver_kwargs_copy)

        # Set the method attribute for compatibility
        solution._method = self.method
        solution._model = self.model

        return solution


# Private alias - not for public use
# Users should use simulator.optimize() instead
_FBA = _RegulatoryAnalysisBase
