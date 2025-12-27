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
Author: Vítor Pereira
##############################################################################
"""

# isort: off
# Import order matters to avoid circular imports
from .simulator import get_container, get_simulator
from .simulation import SimulationMethod, SimulationResult, Simulator, SStatus
from .environment import Environment
from .sglobal import __MEWPY_sim_solvers__

# isort: on

default_solver = None


def get_default_solver():
    """
    Get the currently configured default solver.

    Returns:
        str: Name of the default solver (e.g., 'cplex', 'gurobi', 'glpk')
    """
    global default_solver

    if default_solver:
        return default_solver

    solver_order = ["cplex", "gurobi", "scip", "glpk"]

    for solver in solver_order:
        if solver in __MEWPY_sim_solvers__:
            default_solver = solver
            break

    if not default_solver:
        raise RuntimeError("No solver available.")

    return default_solver


def set_default_solver(solvername):
    """Sets default solver.
    Arguments:
        solvername : (str) solver name (currently available: 'gurobi', 'cplex', 'scip', 'glpk')
    """

    global default_solver

    if solvername.lower() in __MEWPY_sim_solvers__:
        default_solver = solvername.lower()
        # try to also set the local solvers interfaces
        # implementation to the selected solver
        try:
            import mewpy.solvers as msolvers

            msolvers.set_default_solver(solvername)
        except (ImportError, AttributeError) as e:
            import warnings

            warnings.warn(f"Failed to set default solver '{solvername}': {e}")
    else:
        raise RuntimeError(f"Solver {solvername} not available.")
