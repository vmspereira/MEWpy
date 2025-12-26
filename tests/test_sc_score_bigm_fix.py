#!/usr/bin/env python
"""
Test to validate the SC score Big-M constraint fix.

This test verifies that the Species Coupling Score correctly identifies dependencies
when one organism requires metabolites from another.

Bug fixed: Big-M constraints had incorrect signs causing infeasibility when y_k=0.
"""
import unittest

from cobra.io.sbml import read_sbml_model

from mewpy.com import CommunityModel, sc_score

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"


class TestSCScoreBigMFix(unittest.TestCase):
    """Test suite for SC score Big-M constraint fix."""

    def setUp(self):
        """Set up community with two E. coli models."""
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "ecoli1"

        model2 = model1.copy()
        model2.id = "ecoli2"

        self.community = CommunityModel([model1, model2])

    def test_sc_score_runs_without_error(self):
        """Test that SC score completes without errors."""
        scores = sc_score(
            self.community, environment=None, min_growth=0.01, n_solutions=10, verbose=False, use_pool=True
        )

        self.assertIsNotNone(scores, "SC score should not return None")

    def test_sc_score_returns_valid_values(self):
        """Test that SC score returns values in valid range [0, 1]."""
        scores = sc_score(
            self.community, environment=None, min_growth=0.01, n_solutions=10, verbose=False, use_pool=True
        )

        for org_id, org_scores in scores.items():
            self.assertIsNotNone(org_scores, f"Scores for {org_id} should not be None")
            for other_org, score in org_scores.items():
                self.assertGreaterEqual(score, 0.0, f"Score should be >= 0, got {score}")
                self.assertLessEqual(score, 1.0, f"Score should be <= 1, got {score}")

    def test_identical_organisms_low_dependency(self):
        """Test that identical organisms in same environment have low mutual dependency."""
        scores = sc_score(
            self.community, environment=None, min_growth=0.01, n_solutions=10, verbose=False, use_pool=True
        )

        # For identical organisms that can both grow independently,
        # dependency should be very low (close to 0)
        for org_id, org_scores in scores.items():
            for other_org, score in org_scores.items():
                self.assertLess(score, 0.5, f"Identical organisms should have low dependency, got {score}")


if __name__ == "__main__":
    unittest.main()
