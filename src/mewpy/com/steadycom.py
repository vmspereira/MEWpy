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
SteadyCom adapted from REFRAMED to be used with both COBRApy and REFRAMED

Authors: Vitor Pereira
##############################################################################
"""
from math import inf, isinf
from warnings import warn

from mewpy.solvers import solver_instance
from mewpy.solvers.solution import Status, print_balance, print_values
from mewpy.util.constants import ModelConstants
from mewpy.util.utilities import molecular_weight


def calculate_bigM(community, min_value=1000, max_value=1e6, safety_factor=10):
    """
    Calculate an appropriate BigM value for SteadyCom based on model characteristics.

    BigM is used in SteadyCom constraints to handle reactions with infinite bounds.
    The value must be:
    - Large enough to not artificially constrain fluxes
    - Small enough to avoid numerical instability in LP solvers
    - Appropriate for the specific models in the community

    Algorithm:
    1. Find maximum finite flux bound across all organisms
    2. Apply safety factor (default 10x)
    3. Clamp between min_value and max_value

    Args:
        community (CommunityModel): The community model
        min_value (float): Minimum BigM value (default 1000)
        max_value (float): Maximum BigM value to avoid numerical issues (default 1e6)
        safety_factor (float): Multiplier for max bound (default 10)

    Returns:
        float: Calculated BigM value

    Example:
        If max reaction bound is 100, with safety_factor=10:
        BigM = min(1000 * 10, 1e6) = min(10000, 1e6) = 10000
    """
    max_bound = 0.0

    for org_id, organism in community.organisms.items():
        for r_id in organism.reactions:
            reaction = organism.get_reaction(r_id)

            # Check lower bound
            if not isinf(reaction.lb) and abs(reaction.lb) > max_bound:
                max_bound = abs(reaction.lb)

            # Check upper bound
            if not isinf(reaction.ub) and abs(reaction.ub) > max_bound:
                max_bound = abs(reaction.ub)

    # Apply safety factor
    calculated_bigM = max_bound * safety_factor

    # Clamp to reasonable range
    bigM = max(min_value, min(calculated_bigM, max_value))

    return bigM


def SteadyCom(community, constraints=None, solver=None):
    """Implementation of SteadyCom (Chan et al 2017). Adapted from REFRAMED
    Args:
        community (CommunityModel): community model
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
    Returns:
        CommunitySolution: solution object
    """
    # set the proper community building configuration
    community.add_compartments = False
    community.merged_biomasses = False
    community.balance_exchanges = False

    if solver is None:
        solver = build_problem(community)

    objective = community.biomass
    sol = binary_search(solver, objective, minimize=False, constraints=constraints)

    solution = CommunitySolution(community, sol.values)
    solution.solver = solver
    return solution


def SteadyComVA(community, obj_frac=1.0, constraints=None, solver=None):
    """Abundance Variability Analysis using SteadyCom (Chan et al 2017). Adapated from REFRAMED
    Args:
        community (CommunityModel): community model
        obj_frac (float): minimum fraction of the maximum growth rate (default 1.0)
        constraints (dict): environmental or additional constraints (optional)
        solver (Solver): solver instance instantiated with the model, for speed (optional)
    Returns:
        dict: species abundance variability
    """

    community.merged_biomasses = False

    if solver is None:
        solver = build_problem(community)

    objective = {community.biomass: 1}

    sol = binary_search(solver, objective, constraints=constraints)
    growth = obj_frac * sol.values[community.biomass]
    solver.update_growth(growth)

    variability = {org_id: [None, None] for org_id in community.organisms}

    for org_id in community.organisms:
        sol2 = solver.solve({f"x_{org_id}": 1}, minimize=True, get_values=False, constraints=constraints)
        variability[org_id][0] = sol2.fobj

    for org_id in community.organisms:
        sol2 = solver.solve({f"x_{org_id}": 1}, minimize=False, get_values=False, constraints=constraints)
        variability[org_id][1] = sol2.fobj

    return variability


def build_problem(community, growth=1, bigM=None):
    """
    Build the SteadyCom optimization problem.

    Constructs the LP/MILP formulation for SteadyCom as described in Chan et al. 2017.
    The formulation uses Big-M constraints to enforce abundance-scaled flux bounds:
    lb_ij * X_i <= v_ij <= ub_ij * X_i

    Args:
        community (CommunityModel): The community model to optimize
        growth (float): Initial growth rate value for binary search. Defaults to 1.
        bigM (float, optional): Big-M value for reactions with infinite bounds.
            If None (default), automatically calculates based on model characteristics
            using calculate_bigM(). Manual values should be chosen carefully:
            - Too small: artificially constrains fluxes, may cause infeasibility
            - Too large: numerical instability in LP solver
            Recommended: Use automatic calculation (bigM=None)

    Returns:
        Solver: Configured solver instance with SteadyCom problem formulation
            - Variables: x_{org_id} (abundances), reaction fluxes
            - Constraints: abundance sum = 1, mass balance, growth coupling, flux bounds
            - Method: solver.update_growth(value) to update growth parameter

    Note:
        The BigM value is critical for correct results. Different BigM values can
        yield different abundance predictions. The automatic calculation (bigM=None)
        analyzes the model to choose an appropriate value.

    Reference:
        Chan, S. H. J., et al. (2017). SteadyCom: Predicting microbial abundances
        while ensuring community stability. PLoS Computational Biology, 13(5), e1005539.
    """
    # Calculate BigM automatically if not provided
    if bigM is None:
        bigM = calculate_bigM(community)

    # Validate BigM is reasonable
    if bigM < 100:
        warn(
            f"BigM value ({bigM}) is very small and may artificially constrain fluxes. "
            "This could lead to incorrect abundance predictions or infeasibility. "
            "Consider using a larger value or automatic calculation (bigM=None).",
            UserWarning,
            stacklevel=2,
        )
    elif bigM > 1e7:
        warn(
            f"BigM value ({bigM}) is very large and may cause numerical instability. "
            "This could lead to inaccurate solutions. "
            "Consider using a smaller value or automatic calculation (bigM=None).",
            UserWarning,
            stacklevel=2,
        )

    solver = solver_instance()
    community.add_compartments = False
    sim = community.get_community_model()

    # create biomass variables
    for org_id in community.organisms.keys():
        solver.add_variable(f"x_{org_id}", 0, 1, update=False)

    # create all community reactions
    for r_id in sim.reactions:
        reaction = sim.get_reaction(r_id)
        if r_id in sim.get_exchange_reactions():
            solver.add_variable(r_id, reaction.lb, reaction.ub, update=False)
        else:
            # For internal reactions, use tighter bounds based on original reaction bounds
            # Since fluxes are scaled by abundance (v = X * flux) and X <= 1,
            # we can use the original bounds directly (multiplied by max abundance = 1)
            # This provides better numerical conditioning than using (-inf, inf)
            lb = reaction.lb if not isinf(reaction.lb) else (-bigM if reaction.lb < 0 else 0)
            ub = reaction.ub if not isinf(reaction.ub) else (bigM if reaction.ub > 0 else 0)
            solver.add_variable(r_id, lb, ub, update=False)

    solver.update()

    # sum biomass = 1
    solver.add_constraint("abundance", {f"x_{org_id}": 1 for org_id in community.organisms.keys()}, rhs=1, update=False)

    # S.v = 0
    table = sim.metabolite_reaction_lookup()
    for m_id in sim.metabolites:
        solver.add_constraint(m_id, table[m_id], update=False)
    solver.update()
    # organism-specific constraints
    for org_id, organism in community.organisms.items():

        for r_id in organism.reactions:

            if (org_id, r_id) not in community.reaction_map:
                continue

            new_id = community.reaction_map[(org_id, r_id)]
            reaction = organism.get_reaction(r_id)

            # growth = mu * X
            if r_id in organism.objective.keys():
                solver.add_constraint(f"g_{org_id}", {f"x_{org_id}": growth, new_id: -1}, update=False)
            # lb * X < R < ub * X
            else:
                lb = -bigM if isinf(reaction.lb) else reaction.lb
                ub = bigM if isinf(reaction.ub) else reaction.ub

                if lb != 0:
                    solver.add_constraint(f"lb_{new_id}", {f"x_{org_id}": lb, new_id: -1}, "<", 0, update=False)

                if ub != 0:
                    solver.add_constraint(f"ub_{new_id}", {f"x_{org_id}": ub, new_id: -1}, ">", 0, update=False)

    solver.update()

    def update_growth(value):
        coefficients = [(f"g_{x}", f"x_{x}", value) for x in community.organisms]
        solver.change_coefficients(coefficients)

    solver.update_growth = update_growth
    return solver


def binary_search(
    solver,
    objective,
    obj_frac=1,
    minimize=False,
    max_iters=30,
    abs_tol=1e-6,
    rel_tol=1e-4,
    constraints=None,
    raise_on_fail=False,
):
    """
    Binary search to find maximum community growth rate.

    Args:
        solver: Solver instance with update_growth method
        objective: Objective dictionary
        obj_frac: Fraction of optimal growth to use (default 1.0)
        minimize: Whether to minimize objective (default False)
        max_iters: Maximum iterations (default 30)
        abs_tol: Absolute tolerance for convergence (default 1e-6)
        rel_tol: Relative tolerance for convergence (default 1e-4, i.e., 0.01%)
        constraints: Additional constraints
        raise_on_fail: Raise exception on non-convergence (default False for backward compatibility)

    Returns:
        Solution object

    Raises:
        RuntimeError: If raise_on_fail=True and binary search does not converge
        ValueError: If community has no feasible growth (all attempts infeasible)
    """
    previous_value = 0
    value = 1
    fold = 2
    feasible = False
    last_feasible = 0
    converged = False

    for i in range(max_iters):
        diff = value - previous_value

        # Check convergence with both absolute and relative tolerance
        # Use absolute tolerance for small growth rates, relative for large
        if last_feasible > 0:
            rel_diff = abs(diff) / last_feasible
            if abs(diff) < abs_tol or rel_diff < rel_tol:
                converged = True
                break
        elif abs(diff) < abs_tol:
            converged = True
            break

        if feasible:
            last_feasible = value
            previous_value = value
            value = fold * diff + value
        else:
            if i > 0:
                fold = 0.5
            value = fold * diff + previous_value

        solver.update_growth(value)
        sol = solver.solve(objective, get_values=False, minimize=minimize, constraints=constraints)

        feasible = sol.status == Status.OPTIMAL

    # Check if we found any feasible solution
    if last_feasible == 0:
        raise ValueError(
            "Community has no viable growth rate (all attempts infeasible). "
            "Check that organisms can grow and have compatible metabolic capabilities."
        )

    # Final solve at optimal growth rate
    if feasible:
        solver.update_growth(obj_frac * value)
    else:
        solver.update_growth(obj_frac * last_feasible)
    sol = solver.solve(objective, minimize=minimize, constraints=constraints)

    # Handle non-convergence
    if not converged:
        msg = (
            f"Binary search did not converge in {max_iters} iterations. "
            f"Last feasible growth: {last_feasible:.6f}, "
            f"difference: {abs(diff):.2e}, "
            f"relative difference: {abs(diff)/last_feasible if last_feasible > 0 else 'N/A'}. "
            f"Consider increasing max_iters or adjusting tolerance."
        )
        if raise_on_fail:
            raise RuntimeError(msg)
        else:
            warn(msg)

    return sol


class CommunitySolution(object):

    def __init__(self, community, values):
        self.community = community
        self.values = values
        self.growth = None
        self.abundance = None
        self.exchange = None
        self.internal = None
        self.normalized = None
        self._exchange_map = None
        self.parse_values()

    def __str__(self):
        growth = f"Community growth: {self.growth}\n"
        abundance = "\n".join(f"{org_id}\t{val}" for org_id, val in self.abundance.items())
        return growth + abundance

    def __repr__(self):
        return str(self)

    @property
    def exchange_map(self):
        if self._exchange_map is None:
            self._exchange_map = self.compute_exchanges()

        return self._exchange_map

    def parse_values(self):
        model = self.community.merged_model
        reaction_map = self.community.reaction_map

        self.growth = self.values[self.community.biomass]

        self.abundance = {}
        for org_id, organism in self.community.organisms.items():
            growth_i = self.community.organisms_biomass[org_id]
            self.abundance[org_id] = self.values[growth_i] / self.growth

        self.exchange = {r_id: self.values[r_id] for r_id in model.get_exchange_reactions()}

        self.internal = {}
        self.normalized = {}

        for org_id, organism in self.community.organisms.items():

            abundance = self.abundance[org_id]

            fluxes = {
                r_id: self.values[reaction_map[(org_id, r_id)]]
                for r_id in organism.reactions
                if (org_id, r_id) in reaction_map
            }

            rates = {
                r_id: fluxes[r_id] / abundance if abundance > 0 else 0 for r_id in organism.reactions if r_id in fluxes
            }

            self.internal[org_id] = fluxes
            self.normalized[org_id] = rates

    # calculate overall exchanges (organism x metabolite) -> rate

    def compute_exchanges(self):
        sim = self.community.merged_model
        reaction_map = self.community.reaction_map
        exchanges = {}

        for m_id in sim.get_external_metabolites():

            for org_id, organism in self.community.organisms.items():
                rate = 0
                if m_id not in organism.metabolites:
                    continue
                for r_id in organism.get_metabolite_reactions(m_id):
                    if (org_id, r_id) not in reaction_map:
                        continue

                    flux = self.values[reaction_map[(org_id, r_id)]]

                    if flux != 0:
                        coeff = organism.get_reaction(r_id).stoichiometry[m_id]
                        rate += coeff * flux

                if rate != 0:
                    exchanges[(org_id, m_id)] = rate

        return exchanges

    def cross_feeding(self, as_df=True, abstol=1e-6):
        exchanges = self.compute_exchanges()
        cross_all = []

        for m_id in self.community.merged_model.get_external_metabolites():
            r_out = {x: r for (x, m), r in exchanges.items() if m == m_id and r > abstol}
            r_in = {x: -r for (x, m), r in exchanges.items() if m == m_id and -r > abstol}

            total_in = sum(r_in.values())
            total_out = sum(r_out.values())
            total = max(total_in, total_out)

            if total_in > total_out:
                r_out[None] = total_in - total_out
            if total_out > total_in:
                r_in[None] = total_out - total_in

            cross = [(o1, o2, m_id, r1 * r2 / total) for o1, r1 in r_out.items() for o2, r2 in r_in.items()]
            cross_all.extend(cross)

        if as_df:
            from pandas import DataFrame

            cross_all = DataFrame(cross_all, columns=["donor", "receiver", "compound", "rate"])

        return cross_all

    def mass_flow(self, element=None, as_df=False, abstol=1e-6):

        def get_mass(x):
            formula = x.formula
            mw = molecular_weight(formula, element=element)
            return 0.001 * mw

        entities = list(self.community.organisms) + [None]
        flow = {(o1, o2): 0 for o1 in entities for o2 in entities}

        for o1, o2, m_id, rate in self.cross_feeding(as_df=False):
            flow[(o1, o2)] += get_mass(m_id) * rate

        flow = {key: val for key, val in flow.items() if val > abstol}

        if as_df:
            from pandas import DataFrame

            flow = [(o1, o2, val) for (o1, o2), val in flow.items()]
            flow = DataFrame(flow, columns=["donor", "receiver", "flow"])

        return flow

    def print_external_fluxes(self, pattern=None, sort=False, abstol=1e-9):
        print_values(self.exchange, pattern=pattern, sort=sort, abstol=abstol)

    def print_internal_fluxes(self, org_id, normalized=False, pattern=None, sort=False, abstol=1e-9):
        if normalized:
            print_values(self.normalized[org_id], pattern=pattern, sort=sort, abstol=abstol)
        else:
            print_values(self.internal[org_id], pattern=pattern, sort=sort, abstol=abstol)

    def print_external_balance(self, m_id, sort=False, percentage=False, abstol=1e-9):

        print_balance(self.values, m_id, self.community.merged_model, sort=sort, percentage=percentage, abstol=abstol)

    def print_exchanges(self, m_id=None, abstol=1e-9):

        model = self.community.merged_model
        exchange_rxns = set(model.get_exchange_reactions())

        if m_id:
            mets = [m_id]
        else:
            mets = model.get_external_metabolites()

        for m_id in mets:

            entries = []

            ext_rxns = set(model.get_metabolite_reactions(m_id)) & exchange_rxns

            for r_id in ext_rxns:

                flux = self.values[r_id]
                coeff = model.reactions[r_id].stoichiometry[m_id]
                rate = coeff * flux

                if rate > abstol:
                    entries.append(("=> *   ", "in", rate))
                elif rate < -abstol:
                    entries.append(("   * =>", "out", rate))

            for org_id in self.community.organisms:

                if (org_id, m_id) not in self.exchange_map:
                    continue

                rate = self.exchange_map[(org_id, m_id)]
                if rate > abstol:
                    entries.append(("O --> *", org_id, rate))
                elif rate < -abstol:
                    entries.append(("* --> O", org_id, rate))

            if entries:
                print(m_id)
                entries.sort(key=lambda x: x[2])

                for sense, org_id, rate in entries:
                    print(f"[ {sense} ] {org_id:<12} {rate:< 10.6g}")
