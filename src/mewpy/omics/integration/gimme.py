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
GIMME algorithm

Author: Vitor Pereira
Contributors: Paulo Carvalhais
##############################################################################
"""
from copy import deepcopy
from math import inf

from mewpy.cobra.util import convert_to_irreversible
from mewpy.simulation import get_simulator
from mewpy.solvers import solver_instance
from mewpy.solvers.solution import to_simulation_result

from .. import ExpressionSet, Preprocessing


def GIMME(
    model,
    expr,
    biomass=None,
    condition=0,
    cutoff=25,
    growth_frac=0.9,
    constraints=None,
    parsimonious=False,
    build_model=False,
    **kwargs,
):
    """ Run a GIMME simulation [1]_.

    GIMME minimizes usage of lowly expressed reactions while maintaining growth,
    enabling context-specific flux predictions or tissue-specific model generation.

    Arguments:
        model: a REFRAMED or COBRApy model or a MEWpy Simulator.
        expr (ExpressionSet): transcriptomics data or preprocessed coefficients.
        biomass: the biomass reaction identifier (default: uses model's biomass reaction)
        condition: the condition to use in the simulation\
            (default:0, the first condition is used if more than one.)
        cutoff (int): percentile cutoff for low expression (default: 25).
                     Reactions below this percentile are considered lowly expressed.
        growth_frac (float): minimum growth requirement as fraction of wild-type (default: 0.9)
        constraints (dict): additional constraints (optional)
        parsimonious (bool): compute a parsimonious solution (default: False).
                            If True, performs secondary minimization of total flux.
        build_model (bool): if True, returns a tissue-specific model with inactive reactions removed.
                           if False, returns only flux predictions (default: False)

    Returns:
        Solution: solution (or tuple of (solution, model) if build_model=True)

    Notes:
        The algorithm handles reversible reactions differently depending on build_model:
        - build_model=False: Uses solver variables (_p, _n) preserving original model
        - build_model=True: Physically splits reactions for structural model modification

    References
    ----------
    .. [1] Becker, S. and Palsson, B. O. (2008).
           Context-specific metabolic networks are consistent with experiments.
           PLoS Computational Biology, 4(5), e1000082.
           doi:10.1371/journal.pcbi.1000082
    """
    if not build_model:
        sim = get_simulator(model)
    else:
        sim = get_simulator(deepcopy(model))

    if isinstance(expr, ExpressionSet):
        pp = Preprocessing(sim, expr)
        coeffs, threshold = pp.percentile(condition, cutoff=cutoff)
    else:
        coeffs = expr
        threshold = cutoff

    solver = solver_instance(sim)

    if biomass is None:
        biomass = sim.biomass_reaction

    wt_solution = sim.simulate(constraints=constraints)

    if not constraints:
        constraints = {}
    # add growth constraint
    constraints[biomass] = (growth_frac * wt_solution.fluxes[biomass], inf)

    # Make model irreversible to handle expression coefficients properly
    # Two strategies are used depending on whether we're building a tissue-specific model:
    #
    # Strategy 1 (build_model=False): Use solver variables for irreversibility
    #   - Adds _p and _n variables to the solver without modifying the model
    #   - Preserves original model structure
    #   - Solution values need to be reconstructed (net = forward - reverse)
    #   - Used when we only want flux predictions, not a modified model
    #
    # Strategy 2 (build_model=True): Modify model structure
    #   - Physically splits reversible reactions into forward/reverse reactions
    #   - Creates a new model structure with only irreversible reactions
    #   - Reactions can be deleted without affecting constraint definitions
    #   - Used when building a tissue-specific model for further analysis
    if not build_model:
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

    else:
        convert_to_irreversible(sim, inline=True)

    # define the objective
    objective = dict()
    for r_id, val in coeffs.items():
        lb, _ = sim.get_reaction_bounds(r_id)
        if not build_model and lb < 0:
            pos, neg = r_id + "_p", r_id + "_n"
            objective[pos] = val
            objective[neg] = val
        else:
            objective[r_id] = val

    solution = solver.solve(objective, minimize=True, constraints=constraints)

    if not build_model and parsimonious:
        pre_solution = solution

        solver.add_constraint("obj", objective, "=", pre_solution.fobj)
        objective = dict()

        for r_id in sim.reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + "_p", r_id + "_n"
                objective[pos] = 1
                objective[neg] = 1
            else:
                objective[r_id] = 1

        solution = solver.solve(objective, minimize=True, constraints=constraints)
        solver.remove_constraint("obj")
        solution.pre_solution = pre_solution

    if build_model:
        # Build tissue-specific model by removing inactive reactions
        # Activity classification:
        #   0 = Inactive (lowly expressed AND no flux in solution) -> REMOVE
        #   1 = Highly expressed (above threshold) -> KEEP
        #   2 = Active despite low expression (required for biomass) -> KEEP

        # Get original reaction expression for comparison
        if isinstance(expr, ExpressionSet):
            rxn_exp = pp.reactions_expression(condition)
        else:
            # If expr is already coefficients, we need the original expression
            # This is a limitation - coeffs don't contain original expression values
            rxn_exp = {}

        activity = dict()
        for rx_id in sim.reactions:
            activity[rx_id] = 0  # Default: inactive
            # Check if reaction is highly expressed (above threshold)
            if rx_id in rxn_exp and rxn_exp[rx_id] > threshold:
                activity[rx_id] = 1  # Highly expressed
            elif solution.values[rx_id] > 0:
                activity[rx_id] = 2  # Active despite low/unknown expression (needed for growth)

        # Remove reactions with activity = 0 (inactive and not required)
        rx_to_delete = [rx_id for rx_id, v in activity.items() if v == 0]
        sim.remove_reactions(rx_to_delete)
    else:
        # Reconstruct net flux for reversible reactions before deleting split variables
        for r_id in sim.reactions:
            lb, _ = sim.get_reaction_bounds(r_id)
            if lb < 0:
                pos, neg = r_id + "_p", r_id + "_n"
                # Calculate net flux: forward - reverse
                net_flux = solution.values.get(pos, 0) - solution.values.get(neg, 0)
                solution.values[r_id] = net_flux
                # Remove split variables
                if pos in solution.values:
                    del solution.values[pos]
                if neg in solution.values:
                    del solution.values[neg]

    res = to_simulation_result(model, solution.fobj, constraints, sim, solution)
    if hasattr(solution, "pre_solution"):
        res.pre_solution = solution.pre_solution

    if build_model:
        return res, sim
    else:
        return res
