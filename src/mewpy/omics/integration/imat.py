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
    iMAT (Integrative Metabolic Analysis Tool) algorithm [1]_.

    Integrates gene expression data using MILP with binary variables to maximize
    consistency between fluxes and expression levels. Uses the big-M method for
    indicator constraints.

    The algorithm maximizes the number of:
    - Highly expressed reactions with |flux| >= epsilon (active)
    - Lowly expressed reactions with |flux| < epsilon (inactive)

    :param model: a REFRAMED or COBRApy model or a MEWpy Simulator
    :param expr: ExpressionSet or tuple of (low_coeffs, high_coeffs) dicts
    :param constraints: additional constraints (optional)
    :param cutoff: tuple of (low_percentile, high_percentile) for classification.
                   Default (25, 75) means reactions below 25th percentile are
                   "lowly expressed" and above 75th are "highly expressed"
    :param condition: condition index to use from ExpressionSet
    :param epsilon: threshold for considering a reaction "active" (default: 1).
                   A reaction is active if |flux| >= epsilon
    :param build_model: if True, returns a tissue-specific model with inactive reactions removed.
                        if False, returns only flux predictions (default: False)
    :return: Solution (or tuple of (solution, model) if build_model=True)

    Notes:
        - Uses big-M method with M = max(|lb|, |ub|) + 100 for each reaction
        - Handles reversible and irreversible reactions differently
        - Binary variables: y_* for highly expressed (activity), x_* for lowly expressed (inactivity)
        - MILP can be computationally expensive for large models

    Mathematical Formulation:
        Maximize: Σ y_r (highly expressed) + Σ x_r (lowly expressed)

        Subject to:
        - Standard FBA constraints
        - For highly expressed reactions:
          * Reversible: y_fwd=1 forces flux >= ε, y_rev=1 forces flux <= -ε
          * Irreversible: y=1 forces |flux| >= ε
        - For lowly expressed reactions:
          * x=1 forces -ε < flux < ε (near zero)

    References
    ----------
    .. [1] Shlomi, T., Cabili, M. N., Herrgård, M. J., Palsson, B. Ø., & Ruppin, E. (2008).
           Network-based prediction of human tissue-specific metabolism.
           Nature Biotechnology, 26(9), 1003-1010.
           doi:10.1038/nbt.1487
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

    # ========================================================================
    # CORRECTED IMAT FORMULATION USING BIG-M METHOD
    # ========================================================================
    #
    # For highly expressed reactions: y = 1 indicates |flux| >= epsilon (active)
    # For lowly expressed reactions: x = 1 indicates |flux| < epsilon (inactive)
    #
    # Big-M method: Use M > max possible flux to create indicator constraints
    # ========================================================================

    # For highly expressed reactions, add binary variables to reward activity
    # Goal: Maximize number of highly expressed reactions with |flux| >= epsilon
    for r_id, val in high_coeffs.items():
        lb, ub = sim.get_reaction_bounds(r_id)

        # Compute big-M: larger than maximum possible flux
        M = max(abs(lb), abs(ub)) + 100

        # Reversible reaction: can carry flux in either direction
        if lb < 0 and ub > 0:
            # y_forward = 1 indicates forward flux >= epsilon
            y_fwd = "y_" + r_id + "_fwd"
            objective.append(y_fwd)
            solver.add_variable(y_fwd, 0, 1, vartype=VarType.BINARY, update=True)
            # Constraint: flux >= epsilon - M*(1 - y_fwd)
            # When y_fwd = 1: flux >= epsilon (forces forward activity)
            # When y_fwd = 0: flux >= epsilon - M (always satisfied)
            # Using ">" instead of ">=" (equivalent for MILP)
            solver.add_constraint("c" + y_fwd, {r_id: 1, y_fwd: M}, ">", epsilon + M - 0.001, update=False)

            # y_reverse = 1 indicates reverse flux <= -epsilon
            y_rev = "y_" + r_id + "_rev"
            objective.append(y_rev)
            solver.add_variable(y_rev, 0, 1, vartype=VarType.BINARY, update=True)
            # Constraint: flux <= -epsilon + M*(1 - y_rev)
            # When y_rev = 1: flux <= -epsilon (forces reverse activity)
            # When y_rev = 0: flux <= -epsilon + M (always satisfied)
            # Using "<" instead of "<=" (equivalent for MILP)
            solver.add_constraint("c" + y_rev, {r_id: 1, y_rev: -M}, "<", -epsilon - M + 0.001, update=False)

        # Irreversible forward reaction (lb >= 0)
        elif lb >= 0:
            y = "y_" + r_id
            objective.append(y)
            solver.add_variable(y, 0, 1, vartype=VarType.BINARY, update=True)
            # Constraint: flux >= epsilon - M*(1 - y)
            # When y = 1: flux >= epsilon (forces activity)
            # When y = 0: flux >= epsilon - M (always satisfied)
            # Using ">" instead of ">=" (equivalent for MILP)
            solver.add_constraint("c" + y, {r_id: 1, y: M}, ">", epsilon + M - 0.001, update=False)

        # Irreversible reverse reaction (ub <= 0)
        else:  # ub <= 0
            y = "y_" + r_id
            objective.append(y)
            solver.add_variable(y, 0, 1, vartype=VarType.BINARY, update=True)
            # Constraint: flux <= -epsilon + M*(1 - y)
            # When y = 1: flux <= -epsilon (forces activity)
            # When y = 0: flux <= -epsilon + M (always satisfied)
            # Using "<" instead of "<=" (equivalent for MILP)
            solver.add_constraint("c" + y, {r_id: 1, y: -M}, "<", -epsilon - M + 0.001, update=False)

    solver.update()

    # For lowly expressed reactions, add binary variables to reward inactivity
    # Goal: Maximize number of lowly expressed reactions with |flux| < epsilon
    for r_id, val in low_coeffs.items():
        lb, ub = sim.get_reaction_bounds(r_id)

        # Compute big-M: larger than maximum possible flux
        M = max(abs(lb), abs(ub)) + 100

        x_var = "x_" + r_id
        objective.append(x_var)
        solver.add_variable(x_var, 0, 1, vartype=VarType.BINARY, update=True)

        # x = 1 should enforce: -epsilon < flux < epsilon (inactive)
        # Using big-M method:

        # Constraint 1: flux <= epsilon - epsilon*(1 - x) + M*(1 - x)
        # Simplifies to: flux <= M - (M - epsilon)*(1 - x)
        # When x = 1: flux <= epsilon (upper bound for inactivity)
        # When x = 0: flux <= M (always satisfied)
        # Rearranged: flux - M*x <= epsilon - M + M*x => flux <= epsilon + M*(1-x)
        # In solver form: flux + M*x <= epsilon + M
        # Using "<" instead of "<=" (equivalent for MILP)
        solver.add_constraint("c" + x_var + "_upper", {r_id: 1, x_var: -M}, "<", epsilon - M + 0.001, update=False)

        # Constraint 2: flux >= -epsilon + epsilon*(1 - x) - M*(1 - x)
        # Simplifies to: flux >= -M + (M - epsilon)*(1 - x)
        # When x = 1: flux >= -epsilon (lower bound for inactivity)
        # When x = 0: flux >= -M (always satisfied)
        # Rearranged: flux + M*x >= -epsilon + M
        # Using ">" instead of ">=" (equivalent for MILP)
        solver.add_constraint("c" + x_var + "_lower", {r_id: 1, x_var: M}, ">", -epsilon + M - 0.001, update=False)

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
