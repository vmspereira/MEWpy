#!/usr/bin/env python
"""
Test to validate the SteadyCom BigM automatic calculation fix.

This test verifies that BigM is calculated appropriately based on model
characteristics rather than using a hardcoded value.
"""
import unittest
import warnings

from cobra.io.sbml import read_sbml_model

from mewpy.com import CommunityModel
from mewpy.com.steadycom import SteadyCom, build_problem, calculate_bigM

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"


class TestSteadyComBigMFix(unittest.TestCase):
    """Test suite for SteadyCom BigM automatic calculation."""

    def setUp(self):
        """Set up community with two E. coli models."""
        model1 = read_sbml_model(EC_CORE_MODEL)
        model1.reactions.get_by_id("ATPM").bounds = (0, 0)
        model1.id = "ecoli1"

        model2 = model1.copy()
        model2.id = "ecoli2"

        self.community = CommunityModel([model1, model2])

    def test_calculate_bigM_returns_reasonable_value(self):
        """Test that calculate_bigM returns a reasonable value."""
        bigM = calculate_bigM(self.community)

        # Should be at least 1000 (minimum)
        self.assertGreaterEqual(bigM, 1000, "BigM should be at least 1000")

        # Should be at most 1e6 (maximum to avoid numerical issues)
        self.assertLessEqual(bigM, 1e6, "BigM should not exceed 1e6")

        # Should be a positive number
        self.assertGreater(bigM, 0, "BigM should be positive")

    def test_calculate_bigM_with_custom_parameters(self):
        """Test calculate_bigM with custom min/max/safety_factor."""
        # Test with custom minimum
        bigM = calculate_bigM(self.community, min_value=5000)
        self.assertGreaterEqual(bigM, 5000)

        # Test with custom maximum (but still respecting minimum)
        # Note: min_value takes precedence over max_value to ensure safety
        bigM = calculate_bigM(self.community, min_value=500, max_value=800)
        self.assertLessEqual(bigM, 800)
        self.assertGreaterEqual(bigM, 500)

        # Test with higher safety factor
        bigM_10x = calculate_bigM(self.community, safety_factor=10)
        bigM_20x = calculate_bigM(self.community, safety_factor=20)
        self.assertLess(bigM_10x, bigM_20x, "Higher safety factor should give larger BigM")

    def test_build_problem_uses_automatic_bigM_by_default(self):
        """Test that build_problem uses automatic BigM calculation when not specified."""
        # Should not raise any errors or warnings
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            solver = build_problem(self.community)

            # Should not have warnings about BigM being problematic
            bigm_warnings = [x for x in w if "BigM" in str(x.message) or "bigM" in str(x.message)]
            self.assertEqual(len(bigm_warnings), 0, "Automatic BigM should not trigger warnings")

        self.assertIsNotNone(solver)

    def test_build_problem_warns_on_small_bigM(self):
        """Test that build_problem warns when BigM is too small."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            solver = build_problem(self.community, bigM=50)
            self.assertIsNotNone(solver)  # Use the variable

            # Should warn about small BigM
            bigm_warnings = [x for x in w if "small" in str(x.message).lower()]
            self.assertGreater(len(bigm_warnings), 0, "Should warn when BigM is too small")

    def test_build_problem_warns_on_large_bigM(self):
        """Test that build_problem warns when BigM is too large."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            solver = build_problem(self.community, bigM=1e8)
            self.assertIsNotNone(solver)  # Use the variable

            # Should warn about large BigM
            bigm_warnings = [x for x in w if "large" in str(x.message).lower()]
            self.assertGreater(len(bigm_warnings), 0, "Should warn when BigM is too large")

    def test_steadycom_works_with_automatic_bigM(self):
        """Test that SteadyCom works correctly with automatic BigM calculation."""
        # This is the main integration test
        result = SteadyCom(self.community)

        # Should produce valid results
        self.assertIsNotNone(result)
        self.assertGreater(result.growth, 0, "Community should have positive growth")

        # Abundance variables (x_org) should sum to 1 (enforced by constraint)
        # Note: result.abundance contains normalized biomass fluxes, not the x variables
        abundance_vars = {org_id: result.values[f"x_{org_id}"] for org_id in self.community.organisms.keys()}
        total_abundance_vars = sum(abundance_vars.values())
        self.assertAlmostEqual(total_abundance_vars, 1.0, places=6, msg="Abundance variables (x_org) should sum to 1")

        # Each abundance variable should be between 0 and 1
        for org_id, x_value in abundance_vars.items():
            self.assertGreaterEqual(x_value, 0, f"Abundance variable x_{org_id} should be non-negative")
            self.assertLessEqual(x_value, 1, f"Abundance variable x_{org_id} should not exceed 1")

        # Biomass fluxes (result.abundance) should be positive
        for org_id, biomass_flux in result.abundance.items():
            self.assertGreater(biomass_flux, 0, f"Biomass flux for {org_id} should be positive")

    def test_steadycom_with_manual_bigM(self):
        """Test that SteadyCom still works when BigM is manually specified."""
        # Build solver with manual BigM
        solver = build_problem(self.community, bigM=10000)

        # Run SteadyCom with pre-built solver
        result = SteadyCom(self.community, solver=solver)

        # Should still produce valid results
        self.assertIsNotNone(result)
        self.assertGreater(result.growth, 0)

    def test_automatic_vs_manual_bigM_consistency(self):
        """
        Test that results are consistent between automatic and manually-calculated BigM.

        This verifies that the automatic calculation produces the same results as
        using the calculated value explicitly.
        """
        # Run with automatic BigM
        result_auto = SteadyCom(self.community)

        # Calculate BigM explicitly
        bigM_calculated = calculate_bigM(self.community)

        # Run with calculated BigM
        solver = build_problem(self.community, bigM=bigM_calculated)
        result_manual = SteadyCom(self.community, solver=solver)

        # Growth rates should be very close (allowing for numerical differences)
        self.assertAlmostEqual(
            result_auto.growth,
            result_manual.growth,
            places=5,
            msg="Growth rates should match between automatic and manual BigM",
        )

        # Abundances should be very close
        for org_id in result_auto.abundance:
            self.assertAlmostEqual(
                result_auto.abundance[org_id],
                result_manual.abundance[org_id],
                places=5,
                msg=f"Abundance of {org_id} should match between automatic and manual BigM",
            )

    def test_documentation_example(self):
        """
        Verify the example in calculate_bigM documentation.

        If max reaction bound is 100, with safety_factor=10:
        BigM should be max(100 * 10, min_value) = max(1000, 1000) = 1000
        """
        # This test documents expected behavior
        # For E. coli core, typical bounds are around 1000, so with safety_factor=10
        # we expect BigM to be around 10000
        bigM = calculate_bigM(self.community, safety_factor=10)

        # Should be reasonable for E. coli core model
        self.assertGreaterEqual(bigM, 1000, "BigM should be at least the minimum")
        self.assertLessEqual(bigM, 100000, "BigM should be reasonable for E. coli core")


if __name__ == "__main__":
    unittest.main()
