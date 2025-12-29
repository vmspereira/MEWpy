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
Abstract evaluators

Author: Vitor Pereira
##############################################################################
"""
import math
from abc import ABCMeta, abstractmethod


class EvaluationFunction:
    __metaclass__ = ABCMeta

    def __init__(self, maximize: bool = True, worst_fitness: float = 0.0):
        """This abstract class should be extended by all evaluation functions.

        :param maximize: Wether to maximize (True) or minimize (False), defaults to True
        :type maximize: bool, optional
        :param worst_fitness: The worst fitness value, defaults to 0.0
        :type worst_fitness: float, optional
        """
        self.worst_fitness = worst_fitness
        self.maximize = maximize

    @abstractmethod
    def get_fitness(self, simul_results, candidate, **kwargs):
        """Evaluates a candidate

        :param simul_results: (dic) A dictionary of phenotype SimulationResult objects
        :param candidate:  Candidate being evaluated
        :returns: A fitness value.

        """
        raise NotImplementedError

    @abstractmethod
    def method_str(self):
        raise NotImplementedError

    def short_str(self) -> str:
        return self.method_str()

    def __str__(self):
        return self.method_str()

    def __repr__(self):
        """Rich representation showing evaluation function details."""
        lines = []
        lines.append("=" * 60)
        lines.append(f"Evaluation Function: {self.__class__.__name__}")
        lines.append("=" * 60)

        # Method description
        try:
            method = self.method_str()
            if method and len(method) > 0:
                # Wrap long descriptions
                if len(method) > 50:
                    lines.append(f"{'Description:':<20} {method[:50]}...")
                    lines.append(f"{'':<20} (Use .method_str() for full)")
                else:
                    lines.append(f"{'Description:':<20} {method}")
        except:
            pass

        # Optimization direction
        try:
            direction = "Maximize" if self.maximize else "Minimize"
            lines.append(f"{'Direction:':<20} {direction}")
        except:
            pass

        # Worst fitness
        try:
            if self.worst_fitness is not None:
                lines.append(f"{'Worst fitness:':<20} {self.worst_fitness}")
        except:
            pass

        # Required simulations
        try:
            req_sims = self.required_simulations()
            if req_sims and len(req_sims) > 0:
                sims_str = ", ".join(str(s) for s in req_sims)
                lines.append(f"{'Required methods:':<20} {sims_str}")
        except:
            pass

        # Type-specific attributes
        try:
            # For phenotype evaluation functions with targets
            if hasattr(self, "biomass") and self.biomass:
                lines.append(f"{'Biomass:':<20} {self.biomass}")
            if hasattr(self, "product") and self.product:
                lines.append(f"{'Product:':<20} {self.product}")
            if hasattr(self, "target") and self.target:
                lines.append(f"{'Target:':<20} {self.target}")
            if hasattr(self, "substrate") and self.substrate:
                lines.append(f"{'Substrate:':<20} {self.substrate}")
        except:
            pass

        lines.append("=" * 60)
        return "\n".join(lines)

    def _repr_html_(self):
        """Pandas-like HTML representation for Jupyter notebooks."""
        from mewpy.util.html_repr import render_html_table

        rows = []

        # Method description
        try:
            method = self.method_str()
            if method and len(method) > 0:
                # Wrap long descriptions
                if len(method) > 50:
                    rows.append(("Description", method[:50] + "..."))
                    rows.append(("", "(Use .method_str() for full)"))
                else:
                    rows.append(("Description", method))
        except:
            pass

        # Optimization direction
        try:
            direction = "Maximize" if self.maximize else "Minimize"
            rows.append(("Direction", direction))
        except:
            pass

        # Worst fitness
        try:
            if self.worst_fitness is not None:
                rows.append(("Worst fitness", str(self.worst_fitness)))
        except:
            pass

        # Required simulations
        try:
            req_sims = self.required_simulations()
            if req_sims and len(req_sims) > 0:
                sims_str = ", ".join(str(s) for s in req_sims)
                rows.append(("Required methods", sims_str))
        except:
            pass

        # Type-specific attributes
        try:
            # For phenotype evaluation functions with targets
            if hasattr(self, "biomass") and self.biomass:
                rows.append(("Biomass", self.biomass))
            if hasattr(self, "product") and self.product:
                rows.append(("Product", self.product))
            if hasattr(self, "target") and self.target:
                rows.append(("Target", self.target))
            if hasattr(self, "substrate") and self.substrate:
                rows.append(("Substrate", self.substrate))
        except:
            pass

        return render_html_table(f"Evaluation Function: {self.__class__.__name__}", rows)

    @abstractmethod
    def required_simulations(self):
        return None

    @property
    def no_solution(self):
        """
        Value to be returned for worst case evaluation
        """
        if self.worst_fitness is not None:
            res = self.worst_fitness
        elif self.maximize:
            res = -math.inf
        else:
            res = math.inf
        return res

    def __call__(self, simulationResult, candidate, **kwargs):
        return self.get_fitness(simulationResult, candidate, **kwargs)


class PhenotypeEvaluationFunction(EvaluationFunction):

    def __init__(self, maximize=True, worst_fitness=0.0):
        super(PhenotypeEvaluationFunction, self).__init__(maximize=maximize, worst_fitness=worst_fitness)


class KineticEvaluationFunction(EvaluationFunction):

    def __init__(self, maximize=True, worst_fitness=0.0):
        super(KineticEvaluationFunction, self).__init__(maximize=maximize, worst_fitness=worst_fitness)
