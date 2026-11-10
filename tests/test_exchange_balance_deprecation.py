#!/usr/bin/env python
"""
Test to validate the exchange balancing deprecation and demonstrate mass balance issues.

This test verifies that balance_exchange=True triggers a deprecation warning
and documents why this feature violates conservation of mass.
"""

import unittest
import warnings

from cobra.io.sbml import read_sbml_model

# the community oracle: identical members grow as one alone in the open community medium
from test_g_com import check_community_fba

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
            self.assertEqual(
                len(deprecation_warnings), 1, "DeprecationWarning should be raised when setting balance_exchanges=True"
            )
            self.assertIn("deprecated", str(deprecation_warnings[0].message).lower())
            self.assertIn("mass", str(deprecation_warnings[0].message).lower())
        self.assertTrue(community.balance_exchanges)

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

        # the growth of one member alone; shared exchanges balance the members' exchanges
        check_community_fba(self, community, result)

    def test_mass_balance_documentation(self):
        """
        Document the mass balance violation that occurs with balance_exchange=True.

        A member's exchange moves a metabolite from its own compartment to the shared environment. With the
        default (False) one molecule out is one molecule in (coefficients -1 and 1). With balance_exchange=True
        and abundance 0.3, 1 molecule consumed produces only 0.3 molecules: mass is not conserved.
        Abundance scaling is already handled by the merged biomass equation.
        """
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "ecoli1"

        model2 = model1.copy()
        model2.id = "ecoli2"

        abundances = {"ecoli1": 0.3, "ecoli2": 0.7}
        balanced = CommunityModel([model1, model2], abundances=list(abundances.values())).get_community_model()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            community = CommunityModel([model1, model2], abundances=list(abundances.values()), balance_exchange=True)
            unbalanced = community.get_community_model()

        for org_id, abundance in abundances.items():
            exchange = f"R_EX_glc__D_e_{org_id}"
            self.assertEqual(
                dict(balanced.get_reaction(exchange)["stoichiometry"]), {f"M_glc__D_e_{org_id}": -1, "M_glc__D_e": 1}
            )
            self.assertEqual(
                dict(unbalanced.get_reaction(exchange)["stoichiometry"]),
                {f"M_glc__D_e_{org_id}": -1, "M_glc__D_e": abundance},
            )


if __name__ == "__main__":
    unittest.main()
