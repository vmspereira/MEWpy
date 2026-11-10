"""
Steady-state Regulatory Flux Balance Analysis (SRFBA) - Clean Implementation

This module implements SRFBA using RegulatoryExtension only.
No backwards compatibility with legacy GERM models.
"""

import logging
from math import isfinite
from typing import TYPE_CHECKING, Dict, Union
from warnings import warn

from mewpy.germ.algebra import parse_expression
from mewpy.germ.analysis.analysis_utils import regulator_exchanges
from mewpy.germ.analysis.fba import _RegulatoryAnalysisBase
from mewpy.germ.models.regulatory_extension import RegulatoryExtension
from mewpy.solvers import Solution
from mewpy.solvers.solver import VarType
from mewpy.util.constants import ModelConstants

if TYPE_CHECKING:
    from mewpy.germ.variables import Interaction

# Set up logger for this module
logger = logging.getLogger(__name__)

# Margin of the indicator constraints (x > c, metabolite uptake): 10x the default feasibility tolerance of
# SCIP/CPLEX/Gurobi (1e-6), so that a solution within the solver tolerance cannot flip an indicator
INDICATOR_TOLERANCE = 1e-5


def _finite_bounds(lb: float, ub: float):
    """Bounds usable as big-M coefficients: infinite bounds are replaced by the model default bounds."""
    lb, ub = float(lb), float(ub)
    if not isfinite(lb):
        lb = float(ModelConstants.REACTION_LOWER_BOUND)
    if not isfinite(ub):
        ub = float(ModelConstants.REACTION_UPPER_BOUND)
    return lb, ub


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

        # Names already in the problem: GPRs and interactions share symbols and sub-expressions
        self._variable_names = set(self.solver.var_ids)
        self._constraint_names = set(self.solver.constr_ids)
        # Bounds of the variables known to the linearization (reaction fluxes first)
        self._bounds = {}
        for rxn_id in self.model.reactions:
            rxn_data = self._get_reaction(rxn_id)
            self._bounds[rxn_id] = _finite_bounds(rxn_data["lb"], rxn_data["ub"])
        # Target variable -> switch variable that relaxes its regulatory rule (see _linearize_expression)
        self._overrides = {}
        # Coefficients of the regulatory variables: RegulatoryExtension GPR symbols carry none, and ko()
        # fixes them to (0, 0)
        self._coefficients = self._regulatory_coefficients()

        # Build GPR and regulatory interaction constraints
        if self._has_regulatory_network():
            self._link_metabolite_regulators()
            self._build_gprs()
            self._build_interactions()

        # The rows above are added with update=False; push them into the solver problem
        self.solver.update()

        return self

    def _add_variable(self, var_id: str, lb: float, ub: float, vartype: VarType, update: bool = False):
        """Add a variable once. A repeated add would create a duplicate column (SCIP) or raise (CPLEX)."""
        if var_id in self._variable_names:
            return
        self._variable_names.add(var_id)
        self._bounds[var_id] = _finite_bounds(lb, ub)
        self.solver.add_variable(var_id, lb, ub, vartype, update=update)

    def _add_constraint(self, constr_id: str, lhs: Dict, sense: str, rhs: float, update: bool = False):
        """Add a constraint once. Equal names encode equal sub-expressions, so the repeat is redundant."""
        if constr_id in self._constraint_names:
            return
        self._constraint_names.add(constr_id)
        self.solver.add_constraint(constr_id, lhs, sense, rhs, update=update)

    def _regulatory_coefficients(self) -> Dict[str, tuple]:
        """Regulator and target id -> (min, max) of its coefficients."""
        coefficients = {}
        for variables in (getattr(self.model, "targets", None) or {}, dict(self._get_regulators())):
            for var_id, variable in variables.items():
                values = getattr(variable, "coefficients", None)
                if values:
                    coefficients[var_id] = (float(min(values)), float(max(values)))
        return coefficients

    def _link_metabolite_regulators(self):
        """
        A metabolite regulator m is 1 iff its exchange reaction v carries uptake (v < -tol):
        v + (tol + ub) m <= ub      (m = 1 -> v <= -tol)
        v - (lb + tol) m >= -tol    (m = 0 -> v >= -tol)
        A regulator fixed by ko() (coefficients (0, 0)) forbids the uptake.
        """
        tol = INDICATOR_TOLERANCE
        for reg_id, (exchange_id, _) in regulator_exchanges(self.model).items():
            if exchange_id not in self._bounds:
                continue
            lb, ub = self._bounds[exchange_id]
            m_lb, m_ub = self._coefficients.get(reg_id, (0.0, 1.0))
            fixed = m_lb == m_ub
            self._add_variable(reg_id, float(fixed and m_lb > 0), float(not fixed or m_ub > 0), VarType.INTEGER)
            self._add_constraint(f"uptake_{reg_id}_on", {exchange_id: 1.0, reg_id: tol + ub}, "<", ub)
            self._add_constraint(f"uptake_{reg_id}_off", {exchange_id: 1.0, reg_id: -(lb + tol)}, ">", -tol)

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
        self._add_variable(boolean_variable, 0.0, 1.0, VarType.INTEGER, update=False)

        # Add constraints linking reaction flux to boolean variable (big-M needs finite bounds)
        lb, ub = _finite_bounds(
            rxn_data.get("lb", ModelConstants.REACTION_LOWER_BOUND),
            rxn_data.get("ub", ModelConstants.REACTION_UPPER_BOUND),
        )

        # V - Y*Vmax <= 0  ->  V <= Y*Vmax
        if ub != 0:
            self._add_constraint(
                f"gpr_upper_{rxn_id}", {rxn_id: 1.0, boolean_variable: -float(ub)}, "<", 0.0, update=False
            )

        # V - Y*Vmin >= 0  ->  V >= Y*Vmin
        if lb != 0:
            self._add_constraint(
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
            self._add_variable(target_id, lb, ub, var_type, update=False)

            # Add constraints for the regulatory expression (relaxed when the target state is set explicitly)
            self._linearize_expression(target_id, symbolic, override=True)
        except Exception as e:
            # If constraint building fails for this interaction, skip it but log warning
            logger.warning(
                f"Failed to build constraint for interaction targeting '{interaction.target.id}': {e}. "
                f"This regulatory interaction will be skipped."
            )

    def _linearize_expression(self, boolean_variable: str, symbolic, override: bool = False):
        """
        Linearize a boolean expression into solver rows so that boolean_variable takes its value.

        :param boolean_variable: the variable that must take the value of the expression
        :param symbolic: the symbolic expression to linearize
        :param override: if True (regulatory targets), the link is relaxed when the variable is fixed explicitly
        through initial_state in optimize, so that a knockout wins over the regulatory rule
        """
        result = self._linearize_node(symbolic)

        # x = x (autoregulation): a tautology that constrains nothing
        if result == boolean_variable:
            return

        if not override:
            if isinstance(result, str):
                self._add_constraint(f"link_{boolean_variable}", {boolean_variable: 1.0, result: -1.0}, "=", 0.0)
            else:
                self._add_constraint(f"link_{boolean_variable}", {boolean_variable: 1.0}, "=", float(result))
            return

        # |target - result| <= M * switch; the switch is fixed to 0 unless the target state is set explicitly,
        # through initial_state in optimize or by fixing the target's coefficients (ko())
        lb, ub = self._bounds.get(boolean_variable, (0.0, 1.0))
        switch = f"override_{boolean_variable}"
        fixed = float(lb == ub)
        self._add_variable(switch, fixed, fixed, VarType.CONTINUOUS)
        self._overrides[boolean_variable] = switch

        if isinstance(result, str):
            r_lb, r_ub = self._bounds.get(result, (0.0, 1.0))
            big_m = 1.0 + max(abs(lb), abs(ub), abs(r_lb), abs(r_ub))
            self._add_constraint(
                f"link_{boolean_variable}_ub", {boolean_variable: 1.0, result: -1.0, switch: -big_m}, "<", 0.0
            )
            self._add_constraint(
                f"link_{boolean_variable}_lb", {boolean_variable: -1.0, result: 1.0, switch: -big_m}, "<", 0.0
            )
        else:
            value = float(result)
            big_m = 1.0 + max(abs(lb), abs(ub), abs(value))
            self._add_constraint(f"link_{boolean_variable}_ub", {boolean_variable: 1.0, switch: -big_m}, "<", value)
            self._add_constraint(f"link_{boolean_variable}_lb", {boolean_variable: -1.0, switch: -big_m}, "<", -value)

    def _add_binary(self, name: str) -> str:
        self._add_variable(name, 0.0, 1.0, VarType.INTEGER)
        return name

    def _linearize_node(self, node):
        """
        Add the rows of one node of an expression tree and return the name of the variable holding its value,
        or a float when the node is constant. Parents reference their children through these names, so a nested
        operator such as (a & b) | c is linked to the variable of (a & b) itself.

        :param node: a Symbolic node
        :return: variable name or constant value
        """
        if node.is_true or node.is_one:
            return 1.0

        if node.is_false or node.is_zero:
            return 0.0

        if node.is_symbol:
            name = node.key()
            # an existing variable (e.g. a reaction flux used as regulator) keeps its own bounds; a regulatory
            # variable takes its coefficients (fixed to (0, 0) by ko())
            lb, ub = (
                self._bounds.get(name)
                or (None if node.model_variable else self._coefficients.get(name))
                or tuple(float(bound) for bound in node.bounds)
            )
            var_type = VarType.INTEGER if (lb, ub) == (0.0, 1.0) else VarType.CONTINUOUS
            self._add_variable(name, lb, ub, var_type)
            return name

        if node.is_numeric:
            return float(node.value)

        if node.is_and or node.is_or:
            return self._linearize_and_or(node)

        if node.is_not:
            child = self._linearize_node(node.variables[0])
            if not isinstance(child, str):
                return 0.0 if child else 1.0
            name = self._add_binary(f"not_{node.key()}")
            # a = NOT b  ->  a + b = 1
            self._add_constraint(name, {name: 1.0, child: 1.0}, "=", 1.0)
            return name

        if node.is_relational:
            return self._linearize_relational(node)

        raise ValueError(f"unsupported expression node {node.to_string()}")

    def _linearize_and_or(self, node):
        """
        a = AND(x1..xn):  a <= xi for all i,  a >= sum(xi) - (n - 1)
        a = OR(x1..xn):   a >= xi for all i,  a <= sum(xi)
        """
        is_and = node.is_and
        children = [self._linearize_node(child) for child in node.variables]
        constants = [child for child in children if not isinstance(child, str)]
        names = list(dict.fromkeys(child for child in children if isinstance(child, str)))

        if is_and and any(not constant for constant in constants):
            return 0.0
        if not is_and and any(constant for constant in constants):
            return 1.0
        if not names:
            return 1.0 if is_and else 0.0
        if len(names) == 1:
            return names[0]

        name = self._add_binary(f"{'and' if is_and else 'or'}_{node.key()}")
        for i, child in enumerate(names):
            self._add_constraint(f"{name}_{i}", {name: 1.0, child: -1.0}, "<" if is_and else ">", 0.0)

        coefs = {name: 1.0, **{child: -1.0 for child in names}}
        if is_and:
            self._add_constraint(f"{name}_all", coefs, ">", -(len(names) - 1.0))
        else:
            self._add_constraint(f"{name}_any", coefs, "<", 0.0)
        return name

    def _linearize_relational(self, node):
        """
        Indicator a of a comparison between a variable r (bounds [lb, ub]) and a constant c:
        a = (r > c):  a(c + tol - ub) + r <= c + tol  and  a(lb - c - tol) + r >= lb
        a = (r < c):  a(c + tol - lb) + r >= c + tol  and  a(ub - c - tol) + r <= ub
        a = (r == c): only for binary r and c in {0, 1}
        """
        op_l, op_r = node.variables[0], node.variables[1]
        greater = node.is_greater or node.is_greater_equal or node.is_inequality
        less = node.is_less or node.is_less_equal

        if op_l.is_numeric and not op_r.is_numeric:
            # c > r is r < c
            operand, c_val = op_r, float(op_l.value)
            greater, less = less, greater
        else:
            operand, c_val = op_l, float(op_r.value)

        r = self._linearize_node(operand)
        if not isinstance(r, str):
            raise ValueError(f"comparison between constants {node.to_string()}")

        lb, ub = self._bounds.get(r, (0.0, 1.0))
        tol = INDICATOR_TOLERANCE

        if node.is_equal:
            if (lb, ub) != (0.0, 1.0) or c_val not in (0.0, 1.0):
                raise ValueError(f"equality on a non-binary operand {node.to_string()}")
            if c_val == 1.0:
                return r
            name = self._add_binary(f"rel_{node.key()}")
            self._add_constraint(name, {name: 1.0, r: 1.0}, "=", 1.0)
            return name

        if not (greater or less):
            raise ValueError(f"unsupported comparison {node.to_string()}")

        name = self._add_binary(f"rel_{node.key()}")
        if greater:
            self._add_constraint(f"{name}_1", {name: c_val + tol - ub, r: 1.0}, "<", c_val + tol)
            self._add_constraint(f"{name}_2", {name: lb - c_val - tol, r: 1.0}, ">", lb)
        else:
            self._add_constraint(f"{name}_1", {name: c_val + tol - lb, r: 1.0}, ">", c_val + tol)
            self._add_constraint(f"{name}_2", {name: ub - c_val - tol, r: 1.0}, "<", ub)
        return name

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

        # Ids may be given in either SBML form or as original TRN ids
        initial_state = self._ids.keys(initial_state)

        # Apply regulatory constraints via initial_state (new dicts: the caller's must not be mutated).
        # An explicit state of a regulatory target (e.g. a gene knockout) overrides its regulatory rule.
        if self._has_regulatory_network() and initial_state:
            switches = {self._overrides[key]: (1.0, 1.0) for key in initial_state if key in self._overrides}
            constraints = {**(solver_kwargs.get("constraints") or {}), **initial_state, **switches}
            solver_kwargs = {**solver_kwargs, "constraints": constraints}

        # Use the base FBA optimization with regulatory constraints
        solution = super().optimize(solver_kwargs=solver_kwargs, **kwargs)

        if to_solver:
            return solution

        # Convert to Solution if needed
        if not isinstance(solution, Solution):
            minimize = solver_kwargs.get("minimize", self._minimize)
            return Solution.from_solver(method="SRFBA", solution=solution, model=self.model, minimize=minimize)

        return solution
