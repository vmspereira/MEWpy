#!/usr/bin/env python
"""
Test to validate the exchange balancing deprecation and demonstrate mass balance issues.

This test verifies that balance_exchange=True triggers a deprecation warning
and documents why this feature violates conservation of mass.
"""
import unittest
import warnings

from cobra.io.sbml import read_sbml_model

from mewpy.com import CommunityModel

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"


class TestExchangeBalanceDeprecation(unittest.TestCase):
    """Test suite for balance_exchange deprecation."""

    def test_balance_exchange_default_is_false(self):
        """Test that balance_exchange defaults to False."""
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "ecoli1"

        model2 = model1.copy()
        model2.id = "ecoli2"

        # Create community without specifying balance_exchange
        community = CommunityModel([model1, model2])

        # Should default to False
        self.assertFalse(community.balance_exchanges, "balance_exchange should default to False")

    def test_balance_exchange_explicit_false_no_warning(self):
        """Test that balance_exchange=False does not trigger warning."""
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "ecoli1"

        model2 = model1.copy()
        model2.id = "ecoli2"

        # Explicitly set to False - should not warn
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            community = CommunityModel([model1, model2], balance_exchange=False)
            self.assertIsNotNone(community)  # Use the variable

            # Check no DeprecationWarning was raised
            deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            self.assertEqual(
                len(deprecation_warnings), 0, "No deprecation warning should be raised when balance_exchange=False"
            )

    def test_balance_exchange_true_triggers_warning(self):
        """Test that balance_exchange=True triggers DeprecationWarning."""
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "ecoli1"

        model2 = model1.copy()
        model2.id = "ecoli2"

        # Set to True - should warn
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            community = CommunityModel([model1, model2], balance_exchange=True)
            self.assertIsNotNone(community)  # Use the variable

            # Check that DeprecationWarning was raised
            deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            self.assertGreater(
                len(deprecation_warnings), 0, "DeprecationWarning should be raised when balance_exchange=True"
            )

            # Check warning message content
            warning_msg = str(deprecation_warnings[0].message)
            self.assertIn("deprecated", warning_msg.lower())
            self.assertIn("mass", warning_msg.lower())

    def test_balance_exchange_setter_triggers_warning(self):
        """Test that setting balance_exchanges property to True triggers warning."""
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "ecoli1"

        model2 = model1.copy()
        model2.id = "ecoli2"

        # Create with False (no warning)
        community = CommunityModel([model1, model2], balance_exchange=False)

        # Now set to True - should warn
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            community.balance_exchanges = True

            # Check that DeprecationWarning was raised
            deprecation_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            self.assertGreater(
                len(deprecation_warnings), 0, "DeprecationWarning should be raised when setting balance_exchanges=True"
            )

    def test_community_fba_works_with_default_false(self):
        """Test that community FBA works correctly with default balance_exchange=False."""
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "ecoli1"

        model2 = model1.copy()
        model2.id = "ecoli2"

        # Create community with default settings
        community = CommunityModel([model1, model2])
        sim = community.get_community_model()

        # Run FBA
        result = sim.simulate()

        # Should work and produce positive growth
        self.assertIsNotNone(result)
        self.assertGreater(result.objective_value, 0, "Community should grow with balance_exchange=False")

    def test_mass_balance_documentation(self):
        """
        Document the mass balance violation that occurs with balance_exchange=True.

        This test doesn't actually run balance_exchange=True, but documents the issue.
        """
        # DOCUMENTATION OF THE BUG:
        # ========================
        #
        # When balance_exchange=True and add_compartments=True, the code creates
        # transport reactions between organism-specific and shared compartments:
        #
        # Example: Transport glucose from organism A to shared environment
        # Original stoichiometry:
        #   glc_A + (-1) → glc_shared + (1)
        #   (1 molecule consumed produces 1 molecule - MASS BALANCED)
        #
        # With balance_exchange=True and abundance_A=0.3:
        #   glc_A + (-1) → glc_shared + (0.3)
        #   (1 molecule consumed produces only 0.3 molecules - MASS VIOLATED!)
        #
        # Where did the other 0.7 molecules go? They vanished!
        #
        # This violates the fundamental law of conservation of mass and makes
        # the model thermodynamically inconsistent.
        #
        # The correct approach:
        # - Keep stoichiometry 1:1 (mass balanced)
        # - Scale fluxes through bounds or constraints, not stoichiometry
        # - Abundance scaling is already handled by the merged biomass equation

        self.assertTrue(True, "Documentation test - see comments for explanation")


if __name__ == "__main__":
    unittest.main()
