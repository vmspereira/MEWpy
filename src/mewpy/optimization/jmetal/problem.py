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
JMetal Problems

Authors: Vitor Pereira
##############################################################################
"""
import logging
import random
from typing import List, Tuple

from jmetal.core.problem import Problem
from jmetal.core.solution import Solution

from ...util.process import Evaluable
from ..ea import SolutionInterface, dominance_test

logger = logging.getLogger(__name__)

# define EA representation for OU
IntTupple = Tuple[int]


class KOSolution(Solution[int], SolutionInterface):
    """Class representing a KO solution"""

    def __init__(
        self,
        lower_bound: int,
        upper_bound: int,
        number_of_variables: int,
        number_of_objectives: int,
        number_of_constraints: int = 0,
    ):
        super(KOSolution, self).__init__(number_of_variables, number_of_objectives, number_of_constraints)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        # Initialize variables list for jmetalpy 1.9+ compatibility
        self._variables: List[int] = [[] for _ in range(number_of_variables)]

    @property
    def variables(self) -> List[int]:
        """Get the decision variables (jmetalpy 1.9+ compatibility)"""
        return self._variables

    @variables.setter
    def variables(self, values: List[int]):
        """Set the decision variables (jmetalpy 1.9+ compatibility)"""
        self._variables = values

    def __eq__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return self.variables.sort() == solution.variables.sort()
        return False

    # JMetal consideres all problems as minimization
    # Based on pareto dominance

    def __gt__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=False) == 1
        return False

    def __lt__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=False) == -1
        return False

    def __ge__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=False) != -1
        return False

    def __le__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=False) != 1
        return False

    def __copy__(self):
        new_solution = KOSolution(
            self.lower_bound, self.upper_bound, self.number_of_variables, self.number_of_objectives
        )
        new_solution.objectives = self.objectives[:]
        new_solution.variables = self.variables[:]
        new_solution.constraints = self.constraints[:]
        new_solution.attributes = self.attributes.copy()

        return new_solution

    def get_representation(self):
        """
        Returns a set representation of the candidate
        """
        return set(self.variables)

    def get_fitness(self):
        """
        Returns the candidate fitness list
        """
        return self.objectives

    def __str__(self):
        return " ".join((self.variables))


class OUSolution(Solution[IntTupple], SolutionInterface):
    """
    Class representing a Over/Under expression solution.
    """

    def __init__(
        self, lower_bound: List[int], upper_bound: List[int], number_of_variables: int, number_of_objectives: int
    ):
        super(OUSolution, self).__init__(number_of_variables, number_of_objectives)
        self.upper_bound = upper_bound
        self.lower_bound = lower_bound
        # Initialize variables list for jmetalpy 1.9+ compatibility
        self._variables: List[IntTupple] = [[] for _ in range(number_of_variables)]

    @property
    def variables(self) -> List[IntTupple]:
        """Get the decision variables (jmetalpy 1.9+ compatibility)"""
        return self._variables

    @variables.setter
    def variables(self, values: List[IntTupple]):
        """Set the decision variables (jmetalpy 1.9+ compatibility)"""
        self._variables = values

    def __eq__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return self.variables.sort() == solution.variables.sort()
        return False

    # JMetal consideres all problems as minimization

    def __gt__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=False) == 1
        return False

    def __lt__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=False) == -1
        return False

    def __ge__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=False) != -1
        return False

    def __le__(self, solution) -> bool:
        if isinstance(solution, self.__class__):
            return dominance_test(self, solution, maximize=False) != 1
        return False

    def __copy__(self):
        new_solution = OUSolution(
            self.lower_bound, self.upper_bound, self.number_of_variables, self.number_of_objectives
        )
        new_solution.objectives = self.objectives[:]
        new_solution.variables = self.variables[:]
        new_solution.constraints = self.constraints[:]
        new_solution.attributes = self.attributes.copy()

        return new_solution

    def get_fitness(self):
        """
        Returns the candidate fitness list
        """
        return self.objectives

    def __str__(self):
        return " ".join((self.variables))


class JMetalKOProblem(Problem[KOSolution], Evaluable):

    def __init__(self, problem, initial_population):
        """JMetal OU problem. Encapsulates a MEWpy problem
        so that it can be used in jMetal.
        """
        self.problem = problem
        self._number_of_objectives = len(self.problem.fevaluation)
        # Handle different bounder types
        try:
            if hasattr(self.problem.bounder, "upper_bound") and hasattr(self.problem.bounder, "lower_bound"):
                if isinstance(self.problem.bounder.upper_bound, int) and isinstance(
                    self.problem.bounder.lower_bound, int
                ):
                    self._number_of_variables = self.problem.bounder.upper_bound - self.problem.bounder.lower_bound + 1
                else:
                    self._number_of_variables = 100  # Default fallback
            else:
                self._number_of_variables = 100  # Default fallback
        except (AttributeError, TypeError):
            # Bounder doesn't have expected attributes or values aren't compatible
            self._number_of_variables = 100  # Default fallback
        self._number_of_constraints = 0
        self.obj_directions = []
        self.obj_labels = []
        for f in self.problem.fevaluation:
            self.obj_labels.append(str(f))
            if f.maximize:
                self.obj_directions.append(self.MAXIMIZE)
            else:
                self.obj_directions.append(self.MINIMIZE)
        self.initial_population = initial_population
        self.__next_ini_sol = 0

    @property
    def name(self) -> str:
        return self.problem.get_name()

    def number_of_objectives(self) -> int:
        return self._number_of_objectives

    def number_of_variables(self) -> int:
        return self._number_of_variables

    def number_of_constraints(self) -> int:
        return self._number_of_constraints

    def create_solution(self) -> KOSolution:
        solution = None
        flag = False
        while self.__next_ini_sol < len(self.initial_population) and not flag:
            s = self.initial_population[self.__next_ini_sol]
            try:
                solution = self.problem.encode(s)
                flag = True
                self.__next_ini_sol += 1
            except ValueError as e:
                logger.warning("Skipping seed: %s - %s", s, e)
                self.__next_ini_sol += 1
        if not solution:
            solution = self.problem.generator(random)
        new_solution = KOSolution(
            self.problem.bounder.lower_bound,
            self.problem.bounder.upper_bound,
            len(solution),
            self._number_of_objectives,
        )
        new_solution.variables = list(solution)[:]
        return new_solution

    def reset_initial_population_counter(self):
        """Resets the pointer to the next initial population element.
        This strategy is used to overcome the unavailable seeding API in jMetal.
        """
        import random

        random.shuffle(self.initial_population)
        self.__next_ini_sol = 0

    def get_constraints(self, solution):
        return self.problem.decode(set(solution.variables))

    def evaluate(self, solution: KOSolution) -> KOSolution:
        candidate = set(solution.variables)
        p = self.problem.evaluate_solution(candidate)
        for i in range(len(p)):
            # JMetalPy only deals with minimization problems
            if self.obj_directions[i] == self.MAXIMIZE:
                solution.objectives[i] = -1 * p[i]
            else:
                solution.objectives[i] = p[i]
        return solution

    def evaluator(self, candidates, *args):
        res = []
        for candidate in candidates:
            res.append(self.evaluate(candidate))
        return res

    def get_name(self) -> str:
        return self.problem.get_name()

    def build_operators(self):
        from .operators import build_ko_operators

        return build_ko_operators(self.problem)


class JMetalOUProblem(Problem[OUSolution], Evaluable):

    def __init__(self, problem, initial_population=[]):
        """JMetal OU problem. Encapsulates a MEWpy problem
        so that it can be used in jMetal.
        """
        self.problem = problem
        self._number_of_objectives = len(self.problem.fevaluation)
        # Handle different bounder types for OU problems
        try:
            if hasattr(self.problem.bounder, "lower_bound") and isinstance(
                self.problem.bounder.lower_bound, (list, tuple)
            ):
                self._number_of_variables = len(self.problem.bounder.lower_bound)
            else:
                self._number_of_variables = 100  # Default fallback
        except (AttributeError, TypeError):
            # Bounder doesn't have expected attributes or values aren't compatible
            self._number_of_variables = 100  # Default fallback
        self._number_of_constraints = 0
        self.obj_directions = []
        self.obj_labels = []
        for f in self.problem.fevaluation:
            self.obj_labels.append(str(f))
            if f.maximize:
                self.obj_directions.append(self.MAXIMIZE)
            else:
                self.obj_directions.append(self.MINIMIZE)
        self.initial_population = initial_population
        self.__next_ini_sol = 0

    @property
    def name(self) -> str:
        return self.problem.get_name()

    def number_of_objectives(self) -> int:
        return self._number_of_objectives

    def number_of_variables(self) -> int:
        return self._number_of_variables

    def number_of_constraints(self) -> int:
        return self._number_of_constraints

    def create_solution(self) -> OUSolution:
        solution = None
        flag = False
        while self.__next_ini_sol < len(self.initial_population) and not flag:
            s = self.initial_population[self.__next_ini_sol]
            try:
                solution = self.problem.encode(s)
                flag = True
                self.__next_ini_sol += 1
            except ValueError as e:
                logger.warning("Skipping seed: %s - %s", s, e)
                self.__next_ini_sol += 1
        if not solution:
            solution = self.problem.generator(random)
        new_solution = OUSolution(
            self.problem.bounder.lower_bound,
            self.problem.bounder.upper_bound,
            len(solution),
            self._number_of_objectives,
        )
        new_solution.variables = list(solution)[:]
        return new_solution

    def reset_initial_population_counter(self):
        import random

        random.shuffle(self.initial_population)
        self.__next_ini_sol = 0

    def get_constraints(self, solution):
        return self.problem.decode(set(solution.variables))

    def evaluate(self, solution: KOSolution) -> KOSolution:
        candidate = set(solution.variables)
        p = self.problem.evaluate_solution(candidate)
        for i in range(len(p)):
            # JMetalPy only deals with minimization problems
            if self.obj_directions[i] == self.MAXIMIZE:
                solution.objectives[i] = -1 * p[i]
            else:
                solution.objectives[i] = p[i]
        return solution

    def evaluator(self, candidates, *args):
        res = []
        for candidate in candidates:
            res.append(self.evaluate(candidate))
        return res

    def get_name(self) -> str:
        return self.problem.get_name()

    def build_operators(self):
        from .operators import build_ou_operators

        return build_ou_operators(self.problem)
