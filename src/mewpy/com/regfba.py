"""
##############################################################################
Regularized Flux Balance Analysis for communities
Author: Vitor Pereira
##############################################################################
"""

from warnings import warn

from mewpy.simulation import Simulator, SStatus, get_simulator
from mewpy.solvers import solver_instance
from mewpy.solvers.solution import Status, to_simulation_result

from . import CommunityModel


def regComFBA(cmodel, objective=None, maximize=True, constraints=None, obj_frac=0.99):
    """Run a Regularized Flux Balance Analysis simulation:

    Arguments:
        model (CommunityModel): a constraint-based model
        objective (dict: objective coefficients (optional)
        minimize (bool): minimize objective function (False by default)
        constraints (dict): environmental or additional constraints (optional)

    Returns:
        Solution: solution

    Raises:
        RuntimeError: if the LP solver cannot optimize the quadratic objective
    """
    if isinstance(cmodel, CommunityModel):
        sim = cmodel.get_community_model()
    elif isinstance(cmodel, Simulator):
        sim = cmodel
    else:
        sim = get_simulator(cmodel)

    if not hasattr(sim, "community"):
        raise Exception("The model does not seem to be a community model")

    if not objective:
        objective = sim.objective
        if len(objective) == 0:
            warn("Model objective undefined.")

    solver = solver_instance(sim)

    constraints = sim.resolve_keys(constraints, "reaction") if constraints else {}

    if not objective:
        objective = sim.get_objective()

    pre_solution = sim.simulate(objective, maximize=maximize, constraints=constraints)
    if pre_solution.status != SStatus.OPTIMAL:
        return pre_solution

    solver.add_constraint("obj", objective, ">", obj_frac * pre_solution.objective_value)

    solver.update()
    org_bio = list(sim.community.organisms_biomass.values())
    qobjective = {(rid, rid): 1 for rid in org_bio}

    no_qp = (
        f"regComFBA minimizes a quadratic objective, which the {type(solver).__name__} solver interface "
        "cannot optimize; select a solver with quadratic objectives (SCIP, CPLEX or Gurobi) with "
        "mewpy.solvers.set_default_solver"
    )
    try:
        solution = solver.solve(quadratic=qobjective, minimize=True, constraints=constraints)
    except (NotImplementedError, ValueError) as error:
        raise RuntimeError(no_qp) from error
    # an interface that drops the quadratic term returns objective 0 although the growth rates are not 0
    if solution.status == Status.OPTIMAL and not solution.fobj:
        if any(abs(solution.values.get(r_id, 0)) > 1e-9 for r_id in org_bio):
            raise RuntimeError(no_qp)

    result = to_simulation_result(sim, solution.fobj, constraints, sim, solution, regComFBA)

    return result
