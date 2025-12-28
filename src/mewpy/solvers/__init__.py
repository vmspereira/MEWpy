"""
##############################################################################
Linear programming and ODE solvers

Author: Vitor Pereira
##############################################################################
"""

from .ode import KineticConfigurations, ODEMethod, ODEStatus, SolverConfigurations
from .sglobal import __MEWPY_ode_solvers__, __MEWPY_solvers__
from .solution import Solution, Status

# #################################################
# Linear Programming Solvers
# #################################################

default_solver = None


def get_default_solver():

    global default_solver

    if default_solver:
        return default_solver

    solver_order = ["cplex", "gurobi", "scip", "optlang"]

    for solver in solver_order:
        if solver in list(__MEWPY_solvers__.keys()):
            default_solver = solver
            break

    if not default_solver:
        raise RuntimeError("No solver available.")

    return default_solver


def set_default_solver(solvername):
    """Sets default solver.

    Arguments:
        solvername : (str) solver name (currently available: 'gurobi', 'cplex', 'scip', 'optlang')
    """

    global default_solver

    if solvername.lower() in list(__MEWPY_solvers__.keys()):
        default_solver = solvername.lower()
    else:
        raise RuntimeError(f"Solver {solvername} not available.")


def solvers():
    return list(__MEWPY_solvers__.keys())


def solver_instance(model=None):
    """Returns a new instance of the currently selected solver.

    Arguments:
        model : COBRApy/REFRAMED model or a Simulator (optional) -- immediatly instantiate problem with given model

    Returns:
        Solver
    """

    solver = get_default_solver()

    if solver:
        return __MEWPY_solvers__[solver](model)


def is_scip_solver():
    """Check if the current default solver is SCIP.

    SCIP has different performance characteristics than commercial solvers:
    - Requires freeTransform() before modifying problems after solving
    - May benefit from fresh solver instances in repeated optimization scenarios

    Returns:
        bool: True if SCIP is the default solver
    """
    return get_default_solver() == "scip"


def solver_prefers_fresh_instance():
    """Check if the current solver benefits from fresh instances in repeated optimizations.

    Some solvers (like SCIP) have state machine constraints that make repeated
    modifications less efficient. For these solvers, creating fresh instances
    per optimization can be faster and more stable.

    Returns:
        bool: True if solver benefits from fresh instances (currently only SCIP)
    """
    # Currently only SCIP benefits from fresh instances due to freeTransform() overhead
    # CPLEX and Gurobi handle repeated modifications efficiently
    return is_scip_solver()


# #################################################
# ODE solvers
# #################################################


try:
    from .scikits_solver import ScikitsODESolver

    __MEWPY_ode_solvers__["scikits"] = ScikitsODESolver
except ImportError:
    pass


try:
    from .scipy_solver import ScipySolver

    __MEWPY_ode_solvers__["scipy"] = ScipySolver
except ImportError:
    pass


try:
    from .odespy_solver import ODESpySolver

    __MEWPY_ode_solvers__["odespy"] = ODESpySolver
except ImportError:
    pass


default_ode_solver = None


def get_default_ode_solver():
    global default_ode_solver

    if default_ode_solver:
        return default_ode_solver

    ode_solver_order = ["scikits", "scipy", "odespy"]

    for solver in ode_solver_order:
        if solver in list(__MEWPY_ode_solvers__.keys()):
            default_ode_solver = solver
            break

    if not default_ode_solver:
        raise RuntimeError("No solver ODE available.")

    return default_ode_solver


def set_default_ode_solver(solvername):
    """Sets default solver.

    Arguments:
        solvername : (str) solver name (currently available: 'gurobi', 'cplex')
    """

    global default_ode_solver

    if solvername.lower() in list(__MEWPY_ode_solvers__.keys()):
        default_ode_solver = solvername.lower()
    else:
        raise RuntimeError(f"ODE solver {solvername} not available.")


def ode_solvers():
    return list(__MEWPY_ode_solvers__.keys())


def ode_solver_instance(func, method: ODEMethod):
    """Returns a new instance of the currently selected solver.

    Arguments:
        func : a function
        method: a method

    Returns:
        Solver
    """

    solver = get_default_ode_solver()
    if solver:
        return __MEWPY_ode_solvers__[solver](func, method)
