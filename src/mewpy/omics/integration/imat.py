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
from mewpy.util.constants import ModelConstants

from .. import ExpressionSet, Preprocessing


def _big_m(lb, ub, epsilon):
    """Big-M value larger than any feasible |flux| of a reaction.

    Bounds beyond ModelConstants.REACTION_UPPER_BOUND, including infinite ones and placeholders
    such as +-999999, are capped at that value: a larger M lets the solver's integrality tolerance
    leak flux through the indicator constraints, and REFRAMED reads such placeholders as +-inf
    while COBRApy keeps them, so capping also gives both backends the same MILP.
    """
    cap = ModelConstants.REACTION_UPPER_BOUND
    bound = max(min(abs(b), cap) for b in (lb, ub))
    return bound + abs(epsilon) + 100


def iMAT(
    model,
    expr,
    constraints=None,
    cutoff=(25, 75),
    condition=0,
    epsilon=1,
    build_model=False,
    biomass=None,
    growth_frac=None,
):
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
    :param build_model: if True, also returns a tissue-specific model that keeps the objective
                        reactions and every reaction carrying flux in the iMAT solution, so the
                        solution remains feasible in it. If False, returns only flux predictions
                        (default: False)
    :param biomass: biomass reaction used by growth_frac (default: the model's biomass reaction)
    :param growth_frac: minimum growth as a fraction of the wild-type growth. Defaults to 0 (no
                        requirement, as in the original iMAT) when build_model is False, and to 0.9
                        (as in GIMME) when build_model is True, so that the built model can grow.
    :return: Solution (or tuple of (solution, model) if build_model=True). The solution's
             objective_value is the MILP objective: the number of highly expressed reactions
             found active plus lowly expressed reactions found inactive.

    Notes:
        - Uses big-M method with M = max(|lb|, |ub|) + epsilon + 100 for each reaction,
          where bounds beyond ModelConstants.REACTION_UPPER_BOUND (e.g. infinite) are capped at it
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
          * x=1 forces -ε <= flux <= ε (near zero)

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
        # percentile() returns, for each cutoff, the reactions BELOW it: that is the lowly expressed
        # set for the low cutoff, but the highly expressed set lies ABOVE the high cutoff
        (low_coeffs, _), (_, high_threshold) = pp.percentile(condition, cutoff=cutoff)
        rxn_exp = pp.reactions_expression(condition)
        high_coeffs = {r_id: val - high_threshold for r_id, val in rxn_exp.items() if val > high_threshold}
    else:
        low_coeffs, high_coeffs = expr

    constraints = dict(constraints) if constraints else {}

    if growth_frac is None:
        growth_frac = 0.9 if build_model else 0
    if growth_frac > 0:
        if biomass is None:
            biomass = sim.biomass_reaction
        wt_solution = sim.simulate(constraints=constraints)
        constraints[biomass] = (growth_frac * wt_solution.fluxes[biomass], inf)

    solver = solver_instance(sim)

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
        M = _big_m(lb, ub, epsilon)

        # Reversible reaction: can carry flux in either direction
        if lb < 0 and ub > 0:
            # y_forward = 1 indicates forward flux >= epsilon
            y_fwd = "y_" + r_id + "_fwd"
            objective.append(y_fwd)
            solver.add_variable(y_fwd, 0, 1, vartype=VarType.BINARY, update=True)
            # flux >= epsilon - M*(1 - y_fwd)  <=>  flux - M*y_fwd >= epsilon - M
            # y_fwd = 1: flux >= epsilon; y_fwd = 0: flux >= epsilon - M (always satisfied)
            solver.add_constraint("c" + y_fwd, {r_id: 1, y_fwd: -M}, ">", epsilon - M, update=False)

            # y_reverse = 1 indicates reverse flux <= -epsilon
            y_rev = "y_" + r_id + "_rev"
            objective.append(y_rev)
            solver.add_variable(y_rev, 0, 1, vartype=VarType.BINARY, update=True)
            # flux <= -epsilon + M*(1 - y_rev)  <=>  flux + M*y_rev <= M - epsilon
            # y_rev = 1: flux <= -epsilon; y_rev = 0: flux <= M - epsilon (always satisfied)
            solver.add_constraint("c" + y_rev, {r_id: 1, y_rev: M}, "<", M - epsilon, update=False)

        # Irreversible forward reaction (lb >= 0)
        elif lb >= 0:
            y = "y_" + r_id
            objective.append(y)
            solver.add_variable(y, 0, 1, vartype=VarType.BINARY, update=True)
            # flux >= epsilon - M*(1 - y)  <=>  flux - M*y >= epsilon - M
            solver.add_constraint("c" + y, {r_id: 1, y: -M}, ">", epsilon - M, update=False)

        # Irreversible reverse reaction (ub <= 0)
        else:  # ub <= 0
            y = "y_" + r_id
            objective.append(y)
            solver.add_variable(y, 0, 1, vartype=VarType.BINARY, update=True)
            # flux <= -epsilon + M*(1 - y)  <=>  flux + M*y <= M - epsilon
            solver.add_constraint("c" + y, {r_id: 1, y: M}, "<", M - epsilon, update=False)

    solver.update()

    # For lowly expressed reactions, add binary variables to reward inactivity
    # Goal: Maximize number of lowly expressed reactions with |flux| <= epsilon
    for r_id, val in low_coeffs.items():
        lb, ub = sim.get_reaction_bounds(r_id)
        M = _big_m(lb, ub, epsilon)

        x_var = "x_" + r_id
        objective.append(x_var)
        solver.add_variable(x_var, 0, 1, vartype=VarType.BINARY, update=True)

        # x = 1 enforces -epsilon <= flux <= epsilon (inactive); x = 0 leaves the flux free
        # flux <= epsilon + M*(1 - x)  <=>  flux + M*x <= epsilon + M
        solver.add_constraint("c" + x_var + "_upper", {r_id: 1, x_var: M}, "<", epsilon + M, update=False)
        # flux >= -epsilon - M*(1 - x)  <=>  flux - M*x >= -epsilon - M
        solver.add_constraint("c" + x_var + "_lower", {r_id: 1, x_var: -M}, ">", -epsilon - M, update=False)

    solver.update()

    objective_dict = {x: 1 for x in objective}

    solution = solver.solve(objective_dict, minimize=False, constraints=constraints)

    # Build tissue-specific model if requested
    if build_model:
        # Remove only reactions without flux: removing reactions that carry flux (e.g. |flux| < epsilon)
        # makes the solution infeasible in the built model. The objective reactions are always kept:
        # removing them leaves COBRApy with a stale objective that breaks every later simulation.
        keep = set(sim.objective)
        rx_to_delete = [
            r_id
            for r_id in sim.reactions
            if r_id not in keep and abs(solution.values.get(r_id, 0)) <= ModelConstants.TOLERANCE
        ]
        sim.remove_reactions(rx_to_delete)

    res = to_simulation_result(model, solution.fobj, constraints, sim, solution)

    if build_model:
        return res, sim
    else:
        return res
