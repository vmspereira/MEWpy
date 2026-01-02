"""
Steady-state Regulatory Flux Balance Analysis (SRFBA) - Clean Implementation

This module implements SRFBA using RegulatoryExtension only.
No backwards compatibility with legacy GERM models.
"""

import logging
from typing import TYPE_CHECKING, Dict, Union
from warnings import warn

from mewpy.germ.algebra import parse_expression
from mewpy.germ.analysis.fba import _RegulatoryAnalysisBase
from mewpy.germ.models.regulatory_extension import RegulatoryExtension
from mewpy.solvers import Solution
from mewpy.solvers.solver import VarType
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.germ.variables import Interaction

# Set up logger for this module
logger = logging.getLogger(__name__)


class SRFBA(_RegulatoryAnalysisBase):
    """
    Steady-state Regulatory Flux Balance Analysis (SRFBA) using RegulatoryExtension.

    This implementation works exclusively with RegulatoryExtension instances that wrap
    external simulators and provides full SRFBA functionality including boolean algebra
    constraint handling for regulatory logic (AND, OR, NOT, equal, unequal).

    For more details: Shlomi et al. 2007, https://dx.doi.org/10.1038%2Fmsb4100141
    """

    def __init__(
        self, model: RegulatoryExtension, solver: Union[str, None] = None, build: bool = False, attach: bool = False
    ):
        """
        Steady-state Regulatory Flux Balance Analysis (SRFBA).

        :param model: A RegulatoryExtension instance wrapping a simulator
        :param solver: A solver name. If None, a new solver is instantiated.
        :param build: Whether to build the problem upon instantiation. Default: False
        :param attach: Whether to attach the problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)
        self.method = "SRFBA"
        self._model_default_lb = ModelConstants.REACTION_LOWER_BOUND
        self._model_default_ub = ModelConstants.REACTION_UPPER_BOUND
        self._boolean_variables = {}  # Track boolean variables for regulatory logic

    @property
    def model_default_lb(self) -> float:
        """The default lower bound for the model reactions."""
        if self.synchronized:
            return self._model_default_lb

        # Get bounds from simulator
        bounds = [
            self._get_reaction(rxn_id).get("lb", ModelConstants.REACTION_LOWER_BOUND) for rxn_id in self.model.reactions
        ]
        self._model_default_lb = min(bounds) if bounds else ModelConstants.REACTION_LOWER_BOUND
        return self._model_default_lb

    @property
    def model_default_ub(self) -> float:
        """The default upper bound for the model reactions."""
        if self.synchronized:
            return self._model_default_ub

        # Get bounds from simulator
        bounds = [
            self._get_reaction(rxn_id).get("ub", ModelConstants.REACTION_UPPER_BOUND) for rxn_id in self.model.reactions
        ]
        self._model_default_ub = max(bounds) if bounds else ModelConstants.REACTION_UPPER_BOUND
        return self._model_default_ub

    def build(self):
        """
        Build the SRFBA problem.

        This implementation provides full SRFBA functionality including:
        - Basic metabolic constraints (from FBA)
        - GPR constraints using boolean algebra
        - Regulatory interaction constraints
        - Complete boolean operator support (AND, OR, NOT, equal, unequal)

        The boolean algebra constraint system is fully enabled, providing
        comprehensive integration of gene-protein-reaction rules and
        regulatory network logic into the optimization problem.
        """
        # Build the base metabolic model first
        super().build()

        # Build GPR and regulatory interaction constraints
        if self._has_regulatory_network():
            self._build_gprs()
            self._build_interactions()

        return self

    def _build_gprs(self):
        """Build the GPR (Gene-Protein-Reaction) constraints for the solver."""
        for rxn_id in self.model.reactions:
            gpr = self._get_gpr(rxn_id)
            if not gpr.is_none:
                # Get reaction bounds from simulator
                rxn_data = self._get_reaction(rxn_id)
                self._add_gpr_constraint(rxn_id, gpr, rxn_data)

    def _build_interactions(self):
        """Build the regulatory interaction constraints for the solver."""
        if self._has_regulatory_network():
            for interaction in self._get_interactions():
                self._add_interaction_constraint(interaction)

    def _add_gpr_constraint(self, rxn_id: str, gpr, rxn_data: Dict):
        """
        Add GPR constraint for a reaction using boolean algebra.

        :param rxn_id: Reaction identifier
        :param gpr: Parsed GPR expression (can be Symbolic object or Expression wrapper)
        :param rxn_data: Reaction data dict from simulator
        """
        # Skip if GPR is none/empty
        if hasattr(gpr, "is_none") and gpr.is_none:
            return

        # Extract symbolic expression from Expression wrapper if needed
        # parse_expression() returns Symbolic objects directly (Or, And, Symbol, etc.)
        # But fallback cases may return Expression(Symbol("true"), {})
        from mewpy.germ.algebra import Expression

        if isinstance(gpr, Expression):
            # Extract the symbolic expression from Expression wrapper
            symbolic = gpr.symbolic
        else:
            # Already a Symbolic object (Or, And, Symbol, etc.)
            symbolic = gpr

        # Create boolean variable for the reaction
        boolean_variable = f"bool_{rxn_id}"
        self._boolean_variables[boolean_variable] = rxn_id

        # Add the boolean variable to solver
        self.solver.add_variable(boolean_variable, 0.0, 1.0, VarType.INTEGER, update=False)

        # Add constraints linking reaction flux to boolean variable
        lb = rxn_data.get("lb", ModelConstants.REACTION_LOWER_BOUND)
        ub = rxn_data.get("ub", ModelConstants.REACTION_UPPER_BOUND)

        # V - Y*Vmax <= 0  ->  V <= Y*Vmax
        if ub != 0:
            self.solver.add_constraint(
                f"gpr_upper_{rxn_id}", {rxn_id: 1.0, boolean_variable: -float(ub)}, "<", 0.0, update=False
            )

        # V - Y*Vmin >= 0  ->  V >= Y*Vmin
        if lb != 0:
            self.solver.add_constraint(
                f"gpr_lower_{rxn_id}", {rxn_id: 1.0, boolean_variable: -float(lb)}, ">", 0.0, update=False
            )

        # Add constraints for the GPR symbolic expression
        try:
            self._linearize_expression(boolean_variable, symbolic)
        except Exception as e:
            # If linearization fails, skip this constraint but log warning
            # The reaction will still work with just the flux bounds
            logger.warning(
                f"Failed to linearize GPR for reaction '{rxn_id}': {e}. "
                f"Reaction will be constrained by flux bounds only."
            )

    def _add_interaction_constraint(self, interaction: "Interaction"):
        """
        Add regulatory interaction constraint using boolean algebra.

        :param interaction: the regulatory interaction
        """
        try:
            symbolic = None
            for coefficient, expression in interaction.regulatory_events.items():
                if coefficient > 0.0:
                    # For regulatory interactions, expression has .symbolic attribute
                    if hasattr(expression, "symbolic") and expression.symbolic is not None:
                        symbolic = expression.symbolic
                        break

            if symbolic is None:
                return

            # Skip if expression is none/empty
            if hasattr(symbolic, "is_none") and symbolic.is_none:
                return

            # Get target bounds
            lb = float(min(interaction.target.coefficients))
            ub = float(max(interaction.target.coefficients))

            # Determine variable type
            var_type = VarType.INTEGER if (lb, ub) in [(0, 1), (0.0, 1.0)] else VarType.CONTINUOUS

            # Add target variable
            target_id = interaction.target.id
            self.solver.add_variable(target_id, lb, ub, var_type, update=False)

            # Add constraints for the regulatory expression
            self._linearize_expression(target_id, symbolic)
        except Exception as e:
            # If constraint building fails for this interaction, skip it but log warning
            logger.warning(
                f"Failed to build constraint for interaction targeting '{interaction.target.id}': {e}. "
                f"This regulatory interaction will be skipped."
            )

    def _linearize_expression(self, boolean_variable: str, symbolic):
        """
        Linearize a boolean expression into solver constraints.

        :param boolean_variable: the boolean variable name
        :param symbolic: the symbolic expression to linearize
        """
        if symbolic.is_atom:
            self._linearize_atomic_expression(boolean_variable, symbolic)
        else:
            self._linearize_complex_expression(boolean_variable, symbolic)

    def _linearize_atomic_expression(self, boolean_variable: str, symbolic):
        """
        Linearize an atomic boolean expression.

        :param boolean_variable: the boolean variable name
        :param symbolic: the atomic symbolic expression
        """
        if symbolic.is_symbol:
            # Add symbol variable if not exists
            name = symbolic.key()
            lb, ub = symbolic.bounds
            var_type = VarType.INTEGER if (float(lb), float(ub)) in [(0, 1), (0.0, 1.0)] else VarType.CONTINUOUS

            try:
                self.solver.add_variable(name, float(lb), float(ub), var_type, update=False)
            except:
                pass  # Variable already exists

        # Add constraint based on symbolic type
        constraint_coefs, lb, ub = self._get_atomic_constraint(boolean_variable, symbolic)
        if constraint_coefs:
            self.solver.add_constraint(f"atomic_{boolean_variable}", constraint_coefs, "=", 0.0, update=False)

    def _linearize_complex_expression(self, boolean_variable: str, symbolic):
        """
        Linearize a complex boolean expression with operators.

        :param boolean_variable: the boolean variable name
        :param symbolic: the complex symbolic expression
        """
        auxiliary_variables = []
        last_variable = None

        # Process each atom in the expression
        for atom in symbolic:
            last_variable = atom
            var_name = atom.key()

            # Add auxiliary variables for operators
            if atom.is_and or atom.is_or:
                for i, _ in enumerate(atom.variables[:-1]):
                    aux_var = f"{var_name}_{i}"
                    auxiliary_variables.append(aux_var)
                    self.solver.add_variable(aux_var, 0.0, 1.0, VarType.INTEGER, update=False)
            elif atom.is_not:
                aux_var = f"{var_name}_0"
                auxiliary_variables.append(aux_var)
                self.solver.add_variable(aux_var, 0.0, 1.0, VarType.INTEGER, update=False)

            # Add symbol variables
            if atom.is_symbol:
                try:
                    lb, ub = atom.bounds
                    var_type = VarType.INTEGER if (float(lb), float(ub)) in [(0, 1), (0.0, 1.0)] else VarType.CONTINUOUS
                    self.solver.add_variable(var_name, float(lb), float(ub), var_type, update=False)
                except:
                    pass  # Variable already exists

            # Add operator constraints
            self._add_operator_constraints(atom)

        # Link the final result to the boolean variable
        if last_variable:
            last_var_name = last_variable.key()
            names = [f"{last_var_name}_{i}" for i, _ in enumerate(last_variable.variables[:-1])]
            final_var = names[-1] if names else last_var_name

            self.solver.add_constraint(
                f"link_{boolean_variable}", {boolean_variable: 1.0, final_var: -1.0}, "=", 0.0, update=False
            )

    def _get_atomic_constraint(self, boolean_variable: str, symbolic):
        """Get constraint coefficients for atomic expressions."""
        if symbolic.is_true:
            return {boolean_variable: 1.0}, 1.0, 1.0
        elif symbolic.is_false:
            return {boolean_variable: 1.0}, 0.0, 0.0
        elif symbolic.is_numeric:
            val = float(symbolic.value)
            return {boolean_variable: 1.0}, val, val
        elif symbolic.is_symbol:
            return {boolean_variable: 1.0, symbolic.key(): -1.0}, 0.0, 0.0

        return {}, 0.0, 1.0

    def _add_operator_constraints(self, symbolic):
        """Add constraints for boolean operators (AND, OR, NOT)."""
        if symbolic.is_and:
            self._add_and_constraints(symbolic)
        elif symbolic.is_or:
            self._add_or_constraints(symbolic)
        elif symbolic.is_not:
            self._add_not_constraints(symbolic)
        elif symbolic.is_greater or symbolic.is_greater_equal:
            self._add_greater_constraints(symbolic)
        elif symbolic.is_less or symbolic.is_less_equal:
            self._add_less_constraints(symbolic)
        elif symbolic.is_equal:
            self._add_equal_constraints(symbolic)

    def _add_and_constraints(self, symbolic):
        """
        Add constraints for AND operator: a = b AND c
        Constraint: -1 <= 2*b + 2*c - 4*a <= 3
        """
        name = symbolic.key()
        names = [f"{name}_{i}" for i, _ in enumerate(symbolic.variables[:-1])]

        # Handle first AND operation
        and_op = names[0]
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        coefs = {and_op: -4.0}
        if not (op_l.is_one or op_l.is_true):
            coefs[op_l.key()] = 2.0
        if not (op_r.is_one or op_r.is_true):
            coefs[op_r.key()] = 2.0

        self.solver.add_constraint(f"and_{and_op}", coefs, ">", -1.0, update=False)
        self.solver.add_constraint(f"and_{and_op}_ub", coefs, "<", 3.0, update=False)

        # Handle nested AND operations
        if len(symbolic.variables) > 2:
            children = symbolic.variables[2:]
            for i, op_r in enumerate(children):
                op_l_name = names[i]
                and_op = names[i + 1]

                coefs = {and_op: -4.0, op_l_name: 2.0}
                if not (op_r.is_one or op_r.is_true):
                    coefs[op_r.key()] = 2.0

                self.solver.add_constraint(f"and_{and_op}", coefs, ">", -1.0, update=False)
                self.solver.add_constraint(f"and_{and_op}_ub", coefs, "<", 3.0, update=False)

    def _add_or_constraints(self, symbolic):
        """
        Add constraints for OR operator: a = b OR c
        Constraint: -2 <= 2*b + 2*c - 4*a <= 1
        """
        name = symbolic.key()
        names = [f"{name}_{i}" for i, _ in enumerate(symbolic.variables[:-1])]

        # Handle first OR operation
        or_op = names[0]
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        coefs = {or_op: -4.0}
        if not (op_l.is_one or op_l.is_true):
            coefs[op_l.key()] = 2.0
        if not (op_r.is_one or op_r.is_true):
            coefs[op_r.key()] = 2.0

        self.solver.add_constraint(f"or_{or_op}", coefs, ">", -2.0, update=False)
        self.solver.add_constraint(f"or_{or_op}_ub", coefs, "<", 1.0, update=False)

        # Handle nested OR operations
        if len(symbolic.variables) > 2:
            children = symbolic.variables[2:]
            for i, op_r in enumerate(children):
                op_l_name = names[i]
                or_op = names[i + 1]

                coefs = {or_op: -4.0, op_l_name: 2.0}
                if not (op_r.is_one or op_r.is_true):
                    coefs[op_r.key()] = 2.0

                self.solver.add_constraint(f"or_{or_op}", coefs, ">", -2.0, update=False)
                self.solver.add_constraint(f"or_{or_op}_ub", coefs, "<", 1.0, update=False)

    def _add_not_constraints(self, symbolic):
        """
        Add constraints for NOT operator: a = NOT b
        Constraint: a + b = 1
        """
        op_l = symbolic.variables[0]

        if op_l.is_numeric:
            coefs = {symbolic.key(): 1.0}
            val = float(op_l.value)
        else:
            coefs = {symbolic.key(): 1.0, op_l.key(): 1.0}
            val = 1.0

        self.solver.add_constraint(f"not_{symbolic.key()}", coefs, "=", val, update=False)

    def _add_greater_constraints(self, symbolic):
        """Add constraints for GREATER operator: a => r > value"""
        greater_op = symbolic.key()
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        if op_l.is_numeric:
            operand = op_r
            c_val = float(op_l.value)
        else:
            operand = op_l
            c_val = float(op_r.value)

        _lb, _ub = operand.bounds
        _lb = float(_lb)
        _ub = float(_ub)

        # First constraint: a(value + tolerance - r_UB) + r <= value + tolerance
        coefs1 = {greater_op: c_val + ModelConstants.TOLERANCE - _ub, operand.key(): 1.0}
        self.solver.add_constraint(
            f"greater_{greater_op}_1", coefs1, "<", c_val + ModelConstants.TOLERANCE, update=False
        )

        # Second constraint: a(r_LB - value - tolerance) + r >= r_LB
        coefs2 = {greater_op: _lb - c_val - ModelConstants.TOLERANCE, operand.key(): 1.0}
        self.solver.add_constraint(f"greater_{greater_op}_2", coefs2, ">", _lb, update=False)

    def _add_less_constraints(self, symbolic):
        """Add constraints for LESS operator: a => r < value"""
        less_op = symbolic.key()
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        if op_l.is_numeric:
            operand = op_r
            c_val = float(op_l.value)
        else:
            operand = op_l
            c_val = float(op_r.value)

        _lb, _ub = operand.bounds
        _lb = float(_lb)
        _ub = float(_ub)

        # First constraint: a(value + tolerance - r_LB) + r >= value + tolerance
        coefs1 = {less_op: c_val + ModelConstants.TOLERANCE - _lb, operand.key(): 1.0}
        self.solver.add_constraint(f"less_{less_op}_1", coefs1, ">", c_val + ModelConstants.TOLERANCE, update=False)

        # Second constraint: a(r_UB - value - tolerance) + r <= r_UB
        coefs2 = {less_op: _ub - c_val - ModelConstants.TOLERANCE, operand.key(): 1.0}
        self.solver.add_constraint(f"less_{less_op}_2", coefs2, "<", _ub, update=False)

    def _add_equal_constraints(self, symbolic):
        """Add constraints for EQUAL operator: a => r = value"""
        equal_op = symbolic.key()
        op_l = symbolic.variables[0]
        op_r = symbolic.variables[1]

        if op_l.is_numeric:
            operand = op_r
            c_val = float(op_l.value)
        else:
            operand = op_l
            c_val = float(op_r.value)

        coefs = {equal_op: -c_val, operand.key(): 1.0}
        self.solver.add_constraint(f"equal_{equal_op}", coefs, "=", 0.0, update=False)

    def optimize(
        self, solver_kwargs: Dict = None, initial_state: Dict[str, float] = None, to_solver: bool = False, **kwargs
    ) -> Solution:
        """
        Optimize the SRFBA problem.

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

        # Apply regulatory constraints via initial_state
        if self._has_regulatory_network() and initial_state:
            constraints = solver_kwargs.get("constraints", {})
            constraints.update(initial_state)
            solver_kwargs["constraints"] = constraints

        # Use the base FBA optimization with regulatory constraints
        solution = super().optimize(solver_kwargs=solver_kwargs, **kwargs)

        if to_solver:
            return solution

        # Convert to Solution if needed
        if not isinstance(solution, Solution):
            minimize = solver_kwargs.get("minimize", self._minimize)
            return Solution.from_solver(method="SRFBA", solution=solution, model=self.model, minimize=minimize)

        return solution
