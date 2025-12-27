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
Community evaluators

Author: Vitor Pereira
##############################################################################
"""
from mewpy.cobra.com.analysis import *

from .evaluator import PhenotypeEvaluationFunction


class Interaction(PhenotypeEvaluationFunction):

    def __init__(self, community, environment, maximize=True, worst_fitness=0):
        super().__init__(maximize, worst_fitness)
        self.community = community
        self.environment = environment

    def get_fitness(self, simul_results, candidate, **kwargs):
        constraints = kwargs.get("constraints", None)
        if constraints is None:
            return self.worst_fitness

        # Identify modifications by organism
        mutations = {org_id: dict() for org_id in self.community.organisms}
        for k, v in constraints.items():
            org_id, original_id = self.community.reverse_map[k]
            mutations[org_id][original_id] = v

        # Apply the modifications
        community = self.community.copy()
        for org_id, org_mutations in mutations.items():
            sim_org = community.organism[org_id]
            for k, v in org_mutations.items():
                lb, ub = v if isinstance(v, tuple) else v, v
                sim_org.set_reaction_bounds(k, lb, ub)

        sc_scores = sc_score(community, environment=self.environment)  # noqa: F841
        mu_scores = mu_score(community, environment=self.environment)  # noqa: F841
        mp_scores = mp_score(community, environment=self.environment)  # noqa: F841
        mro_scores = mro_score(community, environment=self.environment)  # noqa: F841

        # TODO: Implement scoring function to combine all the score metrics into a single value.
        # Currently returns 0 as placeholder. The combined score should aggregate:
        # - sc_scores: Species contribution scores
        # - mu_scores: Maximum growth scores
        # - mp_scores: Metabolic productivity scores
        # - mro_scores: Metabolic resource overlap scores
        # Consider using weighted sum, min/max aggregation, or other combination strategy.
        score = 0

        return score
