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
E-Flux algorithm

Author: Vitor Pereira
Contributors: Paulo Carvalhais
##############################################################################
"""
from copy import deepcopy

from mewpy.simulation import get_simulator

from .. import ExpressionSet, Preprocessing


def eFlux(
    model,
    expr,
    condition=0,
    scale_rxn=None,
    scale_value=1,
    constraints=None,
    parsimonious=False,
    max_exp=None,
    build_model=False,
    flux_threshold=1e-6,
    **kwargs,
):
    """ Run an E-Flux simulation (Colijn et al, 2009).

    E-Flux scales reaction bounds based on expression levels, enabling
    context-specific flux predictions or tissue-specific model generation.

    :param model: a REFRAMED or COBRApy model or a MEWpy Simulator.
    :param expr (ExpressionSet): transcriptomics data.
    :param condition: the condition to use in the simulation\
                (default:0, the first condition is used if more than one.)
    :param scale_rxn (str): reaction to scale flux vector (optional)
    :param scale_value (float): scaling factor (mandatory if scale_rxn\
            is specified)
    :param constraints (dict): additional constraints (optional)
    :param parsimonious (bool): compute a parsimonious solution (default: False)
    :param max_exp (float): maximum expression value for normalization.\
            If None, uses max from expression data (optional)
    :param build_model (bool): if True, returns a tissue-specific model with lowly\
            expressed reactions removed. if False, returns only flux predictions (default: False)
    :param flux_threshold (float): threshold for removing reactions when build_model=True.\
            Reactions with scaled bounds below this threshold are removed (default: 1e-6)

    :return: Solution (or tuple of (solution, model) if build_model=True)
    """

    # Use deepcopy if building a tissue-specific model (will modify structure)
    if not build_model:
        sim = get_simulator(model)
    else:
        sim = get_simulator(deepcopy(model))

    if isinstance(expr, ExpressionSet):
        pp = Preprocessing(sim, expr)
        rxn_exp = pp.reactions_expression(condition)
    else:
        rxn_exp = expr

    if max_exp is None:
        max_exp = max(rxn_exp.values())

    # Protection against division by zero (all expression values are zero)
    if max_exp == 0:
        # Treat all-zero expression as uniform expression (no scaling)
        max_exp = 1.0

    bounds = {}

    for r_id in sim.reactions:
        val = rxn_exp[r_id] / max_exp if r_id in rxn_exp else 1
        lb, ub = sim.get_reaction_bounds(r_id)
        lb2 = -val if lb < 0 else 0
        ub2 = val if ub > 0 else 0
        bounds[r_id] = (lb2, ub2)

    # User constraints override expression-based bounds
    # These are NOT scaled by expression (applied as absolute values)
    if constraints:
        for r_id, x in constraints.items():
            lb, ub = x if isinstance(x, tuple) else (x, x)
            bounds[r_id] = (lb, ub)

    if parsimonious:
        sol = sim.simulate(constraints=bounds, method="pFBA")
    else:
        sol = sim.simulate(constraints=bounds)

    if scale_rxn is not None:

        if sol.fluxes[scale_rxn] != 0:
            k = abs(scale_value / sol.fluxes[scale_rxn])
        else:
            k = 0

        for r_id, val in sol.fluxes.items():
            sol.fluxes[r_id] = val * k

    # Build tissue-specific model if requested
    if build_model:
        # Remove reactions with very low expression-scaled bounds
        # These are reactions with expression so low that their flux capacity is negligible
        rx_to_delete = []
        for r_id in sim.reactions:
            if r_id in bounds:
                lb, ub = bounds[r_id]
                # Consider reaction inactive if both bounds are near zero
                if max(abs(lb), abs(ub)) < flux_threshold:
                    rx_to_delete.append(r_id)

        sim.remove_reactions(rx_to_delete)

    if build_model:
        return sol, sim
    else:
        return sol
