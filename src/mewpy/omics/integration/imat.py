# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

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
iMAT algorithm

Author: Vitor Pereira
Contributors: Paulo Carvalhais
##############################################################################
"""
from copy import deepcopy
from math import inf

from mewpy.simulation import get_simulator
from mewpy.simulation.simulation import Simulator
from mewpy.solvers import solver_instance
from mewpy.solvers.solution import to_simulation_result
from mewpy.solvers.solver import VarType

from .. import ExpressionSet, Preprocessing


def iMAT(model, expr, constraints=None, cutoff=(25, 75), condition=0, epsilon=1, build_model=False):
    """
    iMAT (Integrative Metabolic Analysis Tool) algorithm.

    Integrates gene expression data using MILP to maximize consistency between
    fluxes and expression levels. Can generate tissue-specific models by removing
    inactive reactions.

    :param model: a REFRAMED or COBRApy model or a MEWpy Simulator
    :param expr: ExpressionSet or tuple of (low_coeffs, high_coeffs) dicts
    :param constraints: additional constraints (optional)
    :param cutoff: tuple of (low_percentile, high_percentile) for classification.
                   Default (25, 75) means reactions below 25th percentile are
                   "lowly expressed" and above 75th are "highly expressed"
    :param condition: condition index to use from ExpressionSet
    :param epsilon: threshold for considering a reaction "active" (flux > epsilon)
    :param build_model: if True, returns a tissue-specific model with inactive reactions removed.
                        if False, returns only flux predictions (default: False)
    :return: Solution (or tuple of (solution, model) if build_model=True)
    """
    # Validate cutoff parameter
    if not isinstance(cutoff, tuple) or len(cutoff) != 2:
        raise ValueError(f"cutoff must be a tuple of (low, high) percentiles, got: {cutoff}")

    low_cutoff, high_cutoff = cutoff

    if not (0 <= low_cutoff < high_cutoff <= 100):
        raise ValueError(f"cutoff must be (low, high) with 0 <= low < high <= 100, got: ({low_cutoff}, {high_cutoff})")

    # Use deepcopy if building a tissue-specific model (will modify structure)
    if not build_model:
        sim = get_simulator(model)
    else:
        sim = get_simulator(deepcopy(model))

    if isinstance(expr, ExpressionSet):
        pp = Preprocessing(sim, expr)
        coeffs, _ = pp.percentile(condition, cutoff=cutoff)
        low_coeffs, high_coeffs = coeffs
    else:
        low_coeffs, high_coeffs = expr

    solver = solver_instance(sim)

    if not constraints:
        constraints = {}

    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + "_p", r_id + "_n"
            solver.add_variable(pos, 0, inf, update=False)
            solver.add_variable(neg, 0, inf, update=False)
    solver.update()

    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + "_p", r_id + "_n"
            solver.add_constraint("c" + pos, {r_id: -1, pos: 1}, ">", 0, update=False)
            solver.add_constraint("c" + neg, {r_id: 1, neg: 1}, ">", 0, update=False)
    solver.update()

    objective = list()

    # For highly expressed reactions, add binary variables to reward activity
    # We want to maximize the number of highly expressed reactions that carry significant flux
    for r_id, val in high_coeffs.items():
        lb, ub = sim.get_reaction_bounds(r_id)

        # Binary variable for positive direction activity (flux away from lower bound)
        # y_pos = 1 is rewarded when flux > lb + epsilon
        pos_cons = lb - epsilon
        pos = "y_" + r_id + "_p"
        objective.append(pos)
        solver.add_variable(pos, 0, 1, vartype=VarType.BINARY, update=True)
        # Constraint: r_id + (lb - epsilon) * y_pos > lb
        # When y_pos = 1: r_id + lb - epsilon > lb => r_id > epsilon (for lb=0)
        # When y_pos = 0: r_id > lb (always satisfied)
        solver.add_constraint("c" + pos, {r_id: 1, pos: pos_cons}, ">", lb, update=False)

        # Binary variable for negative direction activity (flux toward upper bound)
        # y_neg = 1 is rewarded when flux < ub - epsilon
        neg_cons = ub + epsilon
        neg = "y_" + r_id + "_n"
        objective.append(neg)
        solver.add_variable(neg, 0, 1, vartype=VarType.BINARY, update=True)
        # Constraint: r_id + (ub + epsilon) * y_neg < ub
        # When y_neg = 1: r_id + ub + epsilon < ub => r_id < -epsilon (for ub=0)
        # When y_neg = 0: r_id < ub (always satisfied)
        solver.add_constraint("c" + neg, {r_id: 1, neg: neg_cons}, "<", ub, update=False)

    solver.update()

    # For lowly expressed reactions, add binary variables to reward inactivity
    # We want to maximize the number of lowly expressed reactions with near-zero flux
    for r_id, val in low_coeffs.items():
        lb, ub = sim.get_reaction_bounds(r_id)
        x_var = "x_" + r_id
        objective.append(x_var)
        solver.add_variable(x_var, 0, 1, vartype=VarType.BINARY, update=True)

        # Constraints to reward x_var = 1 when flux is near zero
        # Constraint 1: r_id + lb * x_var > lb
        # When x_var = 1: r_id + lb > lb => r_id > 0 (if lb < 0, forces positive flux)
        # When x_var = 0: r_id > lb (always satisfied)
        solver.add_constraint("c" + x_var + "_pos", {r_id: 1, x_var: lb}, ">", lb, update=False)

        # Constraint 2: r_id + ub * x_var < ub
        # When x_var = 1: r_id + ub < ub => r_id < 0 (if ub > 0, forces negative flux)
        # When x_var = 0: r_id < ub (always satisfied)
        # Together: when x_var = 1, forces lb < r_id < 0 or 0 < r_id < ub (near zero)
        solver.add_constraint("c" + x_var + "_neg", {r_id: 1, x_var: ub}, "<", ub, update=False)

    solver.update()

    objective_dict = {x: 1 for x in objective}

    solution = solver.solve(objective_dict, minimize=False, constraints=constraints)

    # Build tissue-specific model if requested
    if build_model:
        # Remove reactions with near-zero flux (inactive reactions)
        # A reaction is considered active if |flux| >= epsilon
        rx_to_delete = []
        for r_id in sim.reactions:
            flux = solution.values.get(r_id, 0)
            if abs(flux) < epsilon:
                rx_to_delete.append(r_id)

        sim.remove_reactions(rx_to_delete)

    res = to_simulation_result(model, None, constraints, sim, solution)

    if build_model:
        return res, sim
    else:
        return res
