# Copyright (C) 2025 Vitor Pereira
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
PySCIPOpt solver interface

Author: Vitor Pereira
##############################################################################
"""
from math import inf
from warnings import warn

from pyscipopt import Model as SCIPModel

from .solution import Solution, Status
from .solver import Parameter, Solver, VarType, default_parameters


class PySCIPOptSolver(Solver):
    """Implements the solver interface using PySCIPOpt (SCIP)."""

    def __init__(self, model=None):
        Solver.__init__(self)
        self.problem = SCIPModel()

        # Map MEWpy status to SCIP status
        self.status_mapping = {
            "optimal": Status.OPTIMAL,
            "unbounded": Status.UNBOUNDED,
            "infeasible": Status.INFEASIBLE,
            "inforunbd": Status.INF_OR_UNB,
        }

        # Map MEWpy variable types to SCIP variable types
        self.vartype_mapping = {VarType.BINARY: "B", VarType.INTEGER: "I", VarType.CONTINUOUS: "C"}

        # SCIP variables and constraints objects
        self._vars = {}
        self._constrs = {}

        # Cache constraint data for reconstruction (needed for SCIP's change_coefficients limitation)
        self._constr_data = {}  # {constr_id: (lhs, sense, rhs)}

        # Caching for efficient updates
        self._cached_lin_obj = {}
        self._cached_sense = None
        self._cached_lower_bounds = {}
        self._cached_upper_bounds = {}
        self._cached_vars = []
        self._cached_constrs = []

        self.set_parameters(default_parameters)
        self.set_logging(False)

        # Additional SCIP parameters for numerical stability and better performance
        # These help with repeated constraint modifications common in deletion analyses
        try:
            # Numerical stability parameters
            self.problem.setParam("numerics/feastol", 1e-6)  # Feasibility tolerance
            self.problem.setParam("numerics/dualfeastol", 1e-7)  # Dual feasibility tolerance
            self.problem.setParam("numerics/epsilon", 1e-9)  # General epsilon for comparisons

            # LP solver parameters for stability
            self.problem.setParam("lp/threads", 1)  # Single-threaded LP for consistency

            # For models with many constraints/variables, increase limits
            self.problem.setParam("limits/memory", 8192)  # Memory limit in MB (8GB)
        except:
            # Older SCIP versions may not support all parameters
            pass

        if model:
            self.build_problem(model)

    def add_variable(self, var_id, lb=-inf, ub=inf, vartype=VarType.CONTINUOUS, update=True):
        """Add a variable to the current problem.

        Arguments:
            var_id (str): variable identifier
            lb (float): lower bound
            ub (float): upper bound
            vartype (VarType): variable type (default: CONTINUOUS)
            update (bool): update problem immediately
        """

        if update:
            self.add_variables([var_id], [lb], [ub], [vartype])
        else:
            self._cached_vars.append((var_id, lb, ub, vartype))

    def add_variables(self, var_ids, lbs, ubs, vartypes):
        """Add multiple variables to the current problem.

        Arguments:
            var_ids (list): variable identifiers
            lbs (list): lower bounds
            ubs (list): upper bounds
            vartypes (list): variable types
        """

        for var_id, lb, ub, vartype in zip(var_ids, lbs, ubs, vartypes):
            # Handle infinities
            lb = None if lb == -inf else lb
            ub = None if ub == inf else ub

            vtype = self.vartype_mapping[vartype]
            var = self.problem.addVar(name=var_id, lb=lb, ub=ub, vtype=vtype)

            self._vars[var_id] = var
            self.var_ids.append(var_id)
            self._cached_lower_bounds[var_id] = lb if lb is not None else -inf
            self._cached_upper_bounds[var_id] = ub if ub is not None else inf
            self._cached_lin_obj[var_id] = 0.0

    def set_variable_bounds(self, var_id, lb, ub):
        """Modify a variable bounds

        Args:
            var_id (str): variable identifier
            lb (float): lower bound
            ub (float): upper bound

        Note:
            SCIP has a strict state machine. After solving, the problem is in a "transformed" state
            where modifications are not allowed. We must call freeTransform() to return to the
            "problem building" state before making changes. This has some performance overhead
            compared to CPLEX/Gurobi which allow modifications in any state.
        """
        if var_id in self._vars:
            # Free the transformed problem to allow modifications
            # SCIP limitation: Can't modify bounds after solving without this
            try:
                self.problem.freeTransform()
            except:
                pass  # Might not be transformed yet

            var = self._vars[var_id]
            if lb is not None:
                lb_val = None if lb == -inf else lb
                self.problem.chgVarLb(var, lb_val if lb_val is not None else -self.problem.infinity())
                self._cached_lower_bounds[var_id] = lb
            if ub is not None:
                ub_val = None if ub == inf else ub
                self.problem.chgVarUb(var, ub_val if ub_val is not None else self.problem.infinity())
                self._cached_upper_bounds[var_id] = ub

    def add_constraint(self, constr_id, lhs, sense="=", rhs=0, update=True):
        """Add a constraint to the current problem.

        Arguments:
            constr_id (str): constraint identifier
            lhs (dict): variables and respective coefficients
            sense (str): constraint sense (any of: '<', '=', '>'; default '=')
            rhs (float): right-hand side of equation (default: 0)
            update (bool): update problem immediately
        """

        if update:
            self.add_constraints([constr_id], [lhs], [sense], [rhs])
        else:
            self._cached_constrs.append((constr_id, lhs, sense, rhs))

    def add_constraints(self, constr_ids, lhs, senses, rhs):
        """Add a list of constraints to the current problem.

        Arguments:
            constr_ids (list): constraint identifiers
            lhs (list): variables and respective coefficients
            senses (list): constraint senses
            rhs (list): right-hand side of equations
        """

        for constr_id, lh, sense, rh in zip(constr_ids, lhs, senses, rhs):
            # Cache constraint data for potential reconstruction
            self._constr_data[constr_id] = (lh.copy(), sense, rh)

            # Build the linear expression
            expr = sum(coeff * self._vars[var_id] for var_id, coeff in lh.items() if var_id in self._vars)

            # Add constraint based on sense
            if sense == "=":
                constr = self.problem.addCons(expr == rh, name=constr_id)
            elif sense == "<":
                constr = self.problem.addCons(expr <= rh, name=constr_id)
            elif sense == ">":
                constr = self.problem.addCons(expr >= rh, name=constr_id)
            else:
                raise ValueError(f"Invalid constraint sense: {sense}")

            self._constrs[constr_id] = constr
            self.constr_ids.append(constr_id)

    def remove_variable(self, var_id):
        """Remove a variable from the current problem.

        Arguments:
            var_id (str): variable identifier
        """
        self.remove_variables([var_id])

    def remove_variables(self, var_ids):
        """Remove variables from the current problem.

        Arguments:
            var_ids (list): variable identifiers
        """

        for var_id in var_ids:
            if var_id in self._vars:
                var = self._vars[var_id]
                self.problem.delVar(var)
                del self._vars[var_id]
                self.var_ids.remove(var_id)
                del self._cached_lower_bounds[var_id]
                del self._cached_upper_bounds[var_id]
                if var_id in self._cached_lin_obj:
                    del self._cached_lin_obj[var_id]

    def remove_constraint(self, constr_id):
        """Remove a constraint from the current problem.

        Arguments:
            constr_id (str): constraint identifier
        """
        self.remove_constraints([constr_id])

    def remove_constraints(self, constr_ids):
        """Remove constraints from the current problem.

        Arguments:
            constr_ids (list): constraint identifiers
        """

        for constr_id in constr_ids:
            if constr_id in self._constrs:
                constr = self._constrs[constr_id]
                self.problem.delCons(constr)
                del self._constrs[constr_id]
                self.constr_ids.remove(constr_id)
                if constr_id in self._constr_data:
                    del self._constr_data[constr_id]

    def update(self):
        """Update internal structure. Used for efficient lazy updating."""

        if self._cached_vars:
            var_ids = [x[0] for x in self._cached_vars]
            lbs = [x[1] for x in self._cached_vars]
            ubs = [x[2] for x in self._cached_vars]
            vartypes = [x[3] for x in self._cached_vars]
            self.add_variables(var_ids, lbs, ubs, vartypes)
            self._cached_vars = []

        if self._cached_constrs:
            constr_ids = [x[0] for x in self._cached_constrs]
            lhs = [x[1] for x in self._cached_constrs]
            senses = [x[2] for x in self._cached_constrs]
            rhs = [x[3] for x in self._cached_constrs]
            self.add_constraints(constr_ids, lhs, senses, rhs)
            self._cached_constrs = []

    def set_objective(self, linear=None, quadratic=None, minimize=True):
        """Set a predefined objective for this problem.

        Args:
            linear (str or dict): linear coefficients (or a single variable to optimize)
            quadratic (dict): quadratic coefficients (optional)
            minimize (bool): solve a minimization problem (default: True)

        Notes:
            Setting the objective is optional. It can also be passed directly when calling **solve**.
        """

        if quadratic:
            warn("PySCIPOpt solver does not fully support quadratic objectives in this interface.")

        if linear:
            if isinstance(linear, str):
                linear = {linear: 1.0}

            # Free the transformed problem to allow modifications
            try:
                self.problem.freeTransform()
            except:
                pass  # Might not be transformed yet

            # Build objective expression
            obj_expr = sum(coeff * self._vars[var_id] for var_id, coeff in linear.items() if var_id in self._vars)

            # Set objective
            sense = "minimize" if minimize else "maximize"
            self.problem.setObjective(obj_expr, sense)

            self._cached_lin_obj.update(linear)
            self._cached_sense = minimize

            # Check for undeclared variables
            for var_id in linear:
                if var_id not in self._vars:
                    warn(f"Objective variable not previously declared: {var_id}")

    def solve(
        self,
        linear=None,
        quadratic=None,
        minimize=None,
        model=None,
        constraints=None,
        get_values=True,
        shadow_prices=False,
        reduced_costs=False,
        pool_size=0,
        pool_gap=None,
    ):
        """Solve the optimization problem.

        Arguments:
            linear (str or dict): linear objective (optional)
            quadratic (dict): quadratic objective (optional)
            minimize (bool): solve a minimization problem (default: True)
            model: model (optional, leave blank to reuse previous model structure)
            constraints (dict): additional constraints (optional)
            get_values (bool or list): set to false for speedup (default: True)
            shadow_prices (bool): return shadow prices if available (default: False)
            reduced_costs (bool): return reduced costs if available (default: False)
            pool_size (int): calculate solution pool (SCIP supports this)
            pool_gap (float): maximum relative gap for solutions in pool (optional)

        Returns:
            Solution: solution
        """

        if model:
            self.build_problem(model)

        if constraints:
            temp_constrs = self._apply_temporary_constraints(constraints)

        if minimize is not None or linear is not None:
            self.set_objective(linear, quadratic, minimize if minimize is not None else True)

        # Solve the problem
        self.problem.optimize()

        # Get status
        status_str = self.problem.getStatus()
        status = self.status_mapping.get(status_str, Status.UNKNOWN)
        message = status_str

        if status == Status.OPTIMAL:
            fobj = self.problem.getObjVal()
            values, s_prices, r_costs = None, None, None

            if get_values:
                try:
                    if isinstance(get_values, list):
                        values = {
                            var_id: self.problem.getVal(self._vars[var_id])
                            for var_id in get_values
                            if var_id in self._vars
                        }
                    else:
                        values = {var_id: self.problem.getVal(var) for var_id, var in self._vars.items()}
                except Exception:
                    values = {var_id: self.problem.getVal(var) for var_id, var in self._vars.items()}

            if shadow_prices:
                # SCIP provides dual values for linear constraints
                s_prices = {}
                for constr_id, constr in self._constrs.items():
                    try:
                        s_prices[constr_id] = self.problem.getDualsolLinear(constr)
                    except:
                        s_prices[constr_id] = 0.0

            if reduced_costs:
                # SCIP provides reduced costs for variables
                r_costs = {}
                for var_id, var in self._vars.items():
                    try:
                        r_costs[var_id] = self.problem.getVarRedcost(var)
                    except:
                        r_costs[var_id] = 0.0

            solution = Solution(status, message, fobj, values, s_prices, r_costs)
        else:
            solution = Solution(status, message)

        if constraints:
            self._remove_temporary_constraints(temp_constrs)

        return solution

    def _apply_temporary_constraints(self, constraints):
        """Apply temporary constraints and return them for later removal."""
        temp_constrs = []

        for var_id, bounds in constraints.items():
            if var_id in self._vars:
                lb, ub = bounds if isinstance(bounds, tuple) else (bounds, bounds)

                # Store original bounds
                orig_lb = self._cached_lower_bounds[var_id]
                orig_ub = self._cached_upper_bounds[var_id]

                # Apply new bounds
                self.set_variable_bounds(var_id, lb, ub)
                temp_constrs.append((var_id, orig_lb, orig_ub))
            else:
                warn(f"Constrained variable not previously declared: {var_id}")

        return temp_constrs

    def _remove_temporary_constraints(self, temp_constrs):
        """Restore original bounds after temporary constraints."""
        for var_id, orig_lb, orig_ub in temp_constrs:
            self.set_variable_bounds(var_id, orig_lb, orig_ub)

    def set_parameter(self, parameter, value):
        """Set a parameter value for this optimization problem

        Arguments:
            parameter (Parameter): parameter type
            value (float): parameter value
        """

        parameter_mapping = {
            Parameter.TIME_LIMIT: ("limits/time", value),
            Parameter.FEASIBILITY_TOL: ("numerics/feastol", value),
            Parameter.OPTIMALITY_TOL: ("numerics/dualfeastol", value),
            Parameter.MIP_REL_GAP: ("limits/gap", value),
        }

        if parameter in parameter_mapping:
            param_name, param_value = parameter_mapping[parameter]
            self.problem.setParam(param_name, param_value)
        else:
            warn(f"Parameter {parameter} not yet supported for PySCIPOpt.")

    def set_logging(self, enabled=False):
        """Enable or disable log output:

        Arguments:
            enabled (bool): turn logging on (default: False)
        """

        if not enabled:
            self.problem.hideOutput()
        else:
            self.problem.hideOutput(False)

    def write_to_file(self, filename):
        """Write problem to file:

        Arguments:
            filename (str): file path
        """

        self.problem.writeProblem(filename)

    def change_coefficients(self, coefficients):
        """Changes variables coefficients in constraints

        :param coefficients: A list of tuples (constraint name, variable name, new value)
        :type coefficients: list

        Note: SCIP doesn't support modifying constraints after solving,
        so we free the transform, delete and recreate constraints with new coefficients.
        """
        # Free the transformed problem to allow modifications
        try:
            self.problem.freeTransform()
        except:
            pass  # Might not be transformed yet

        # Group changes by constraint
        changes_by_constr = {}
        for constr_id, var_id, new_value in coefficients:
            if constr_id not in changes_by_constr:
                changes_by_constr[constr_id] = {}
            changes_by_constr[constr_id][var_id] = new_value

        # For each constraint that needs modification
        for constr_id, var_changes in changes_by_constr.items():
            if constr_id not in self._constrs or constr_id not in self._constr_data:
                continue

            # Get the cached constraint data
            lhs, sense, rhs = self._constr_data[constr_id]

            # Update the coefficients in the LHS
            new_lhs = lhs.copy()
            for var_id, new_value in var_changes.items():
                new_lhs[var_id] = new_value

            # Delete the old constraint
            old_constr = self._constrs[constr_id]
            self.problem.delCons(old_constr)

            # Update cache
            self._constr_data[constr_id] = (new_lhs, sense, rhs)

            # Build new expression
            expr = sum(coeff * self._vars[var_id] for var_id, coeff in new_lhs.items() if var_id in self._vars)

            # Recreate constraint
            if sense == "=":
                new_constr = self.problem.addCons(expr == rhs, name=constr_id)
            elif sense == "<":
                new_constr = self.problem.addCons(expr <= rhs, name=constr_id)
            elif sense == ">":
                new_constr = self.problem.addCons(expr >= rhs, name=constr_id)
            else:
                raise ValueError(f"Invalid constraint sense: {sense}")

            # Update constraint reference
            self._constrs[constr_id] = new_constr
