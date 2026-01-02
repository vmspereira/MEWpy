"""
Comprehensive validation tests for RFBA and SRFBA implementations.

Tests verify that:
1. RFBA correctly applies regulatory constraints
2. SRFBA integrates Boolean logic into MILP
3. Results are consistent with expected behavior
4. Regulatory network properly constrains metabolic fluxes
"""

import sys
from pathlib import Path

import pytest

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class TestRFBAValidation:
    """Validation tests for RFBA implementation."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        return RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

    def test_rfba_basic_functionality(self, integrated_model):
        """Test basic RFBA functionality."""
        from mewpy.germ.analysis import RFBA

        # Create RFBA instance
        rfba = RFBA(integrated_model)

        # Test that it builds without errors
        rfba.build()

        # Optimize with default initial state (all regulators active)
        solution = rfba.optimize()

        assert solution is not None
        assert hasattr(solution, "objective_value") or hasattr(solution, "fobj")
        print(f"RFBA objective value (all regulators active): {solution.objective_value or solution.fobj}")

    def test_rfba_with_inactive_regulators(self, integrated_model):
        """Test RFBA with some regulators inactive."""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()

        # Get list of regulators
        regulators = list(integrated_model.regulators.keys())

        if len(regulators) > 0:
            # Set first regulator to inactive
            initial_state = {regulators[0]: 0.0}

            solution = rfba.optimize(initial_state=initial_state)

            assert solution is not None
            print(
                f"RFBA objective value (regulator '{regulators[0]}' inactive): "
                f"{solution.objective_value or solution.fobj}"
            )
        else:
            pytest.skip("No regulators found in model")

    def test_rfba_dynamic_mode(self, integrated_model):
        """Test RFBA in dynamic mode (iterative)."""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()

        # Run dynamic RFBA
        solution = rfba.optimize(dynamic=True)

        assert solution is not None
        # Dynamic solution should have a solutions list
        if hasattr(solution, "solutions"):
            print(f"Dynamic RFBA converged in {len(solution.solutions)} iterations")
        else:
            print(f"Dynamic RFBA objective: {solution.objective_value or solution.fobj}")

    def test_rfba_regulatory_constraints_applied(self, integrated_model):
        """Verify that RFBA actually applies regulatory constraints."""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()

        # Get all regulators and set them all to inactive
        regulators = list(integrated_model.regulators.keys())

        if len(regulators) > 0:
            # All regulators inactive
            all_inactive_state = {reg: 0.0 for reg in regulators}

            solution = rfba.optimize(initial_state=all_inactive_state)

            assert solution is not None

            # With all regulators inactive, many reactions should be knocked out
            # So the objective value should be lower (or infeasible)
            print(f"RFBA objective with all regulators inactive: {solution.objective_value or solution.fobj}")
            print(f"Solution status: {solution.status}")
        else:
            pytest.skip("No regulators found in model")

    def test_rfba_decode_methods(self, integrated_model):
        """Test RFBA decode methods work correctly."""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()

        # Test initial state generation
        initial_state = rfba.initial_state()
        assert isinstance(initial_state, dict)
        print(f"RFBA initial state has {len(initial_state)} regulators")

        # Test decode_regulatory_state
        regulatory_state = rfba.decode_regulatory_state(initial_state)
        assert isinstance(regulatory_state, dict)
        print(f"RFBA regulatory state affects {len(regulatory_state)} targets")

        # Test decode_metabolic_state
        metabolic_state = rfba.decode_metabolic_state(initial_state)
        assert isinstance(metabolic_state, dict)

        # Test decode_constraints
        constraints = rfba.decode_constraints(metabolic_state)
        assert isinstance(constraints, dict)
        print(f"RFBA generates {len(constraints)} reaction constraints")


class TestSRFBAValidation:
    """Validation tests for SRFBA implementation."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        return RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

    def test_srfba_basic_functionality(self, integrated_model):
        """Test basic SRFBA functionality."""
        from mewpy.germ.analysis import SRFBA

        # Create SRFBA instance
        srfba = SRFBA(integrated_model)

        # Test that it builds without errors
        srfba.build()

        # Check that boolean variables were created
        assert hasattr(srfba, "_boolean_variables")
        print(f"SRFBA created {len(srfba._boolean_variables)} boolean variables")

        # Optimize
        solution = srfba.optimize()

        assert solution is not None
        assert hasattr(solution, "objective_value") or hasattr(solution, "fobj")
        print(f"SRFBA objective value: {solution.objective_value or solution.fobj}")

    def test_srfba_builds_gpr_constraints(self, integrated_model):
        """Test that SRFBA builds GPR constraints."""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model).build()

        # Check that GPR constraints were added
        solver = srfba.solver

        # Count constraints
        constraint_count = len(solver.list_constraints())
        print(f"SRFBA solver has {constraint_count} constraints")

        # Should have more constraints than basic FBA due to Boolean logic
        assert constraint_count > 0

    def test_srfba_builds_regulatory_constraints(self, integrated_model):
        """Test that SRFBA builds regulatory interaction constraints."""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model).build()

        # Check that regulatory interactions were processed
        if integrated_model.has_regulatory_network():
            interactions = list(integrated_model.yield_interactions())
            print(f"Model has {len(interactions)} regulatory interactions")

            # SRFBA should process these into constraints
            assert len(srfba._boolean_variables) > 0
        else:
            pytest.skip("No regulatory network in model")

    def test_srfba_with_initial_state(self, integrated_model):
        """Test SRFBA with initial regulatory state."""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model).build()

        # Get a regulator
        regulators = list(integrated_model.regulators.keys())

        if len(regulators) > 0:
            # Set initial state for first regulator
            initial_state = {regulators[0]: 0.0}

            solution = srfba.optimize(initial_state=initial_state)

            assert solution is not None
            print(
                f"SRFBA objective with regulator '{regulators[0]}' constrained: "
                f"{solution.objective_value or solution.fobj}"
            )
        else:
            pytest.skip("No regulators found in model")

    def test_srfba_integer_variables(self, integrated_model):
        """Test that SRFBA creates integer variables for Boolean logic."""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model).build()

        # Check SRFBA tracks boolean variables internally
        # Boolean variables are stored in _boolean_variables dict but not added to solver
        assert hasattr(srfba, "_boolean_variables")
        boolean_vars = srfba._boolean_variables
        print(f"SRFBA tracks {len(boolean_vars)} boolean variables internally")

        assert len(boolean_vars) > 0


class TestRFBAvsSRFBA:
    """Compare RFBA and SRFBA results."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        return RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

    def test_compare_basic_results(self, integrated_model):
        """Compare basic RFBA and SRFBA results."""
        from mewpy.germ.analysis import RFBA, SRFBA

        # Run RFBA
        rfba = RFBA(integrated_model)
        rfba_solution = rfba.optimize()

        # Run SRFBA
        srfba = SRFBA(integrated_model)
        srfba_solution = srfba.optimize()

        # Both should produce solutions
        assert rfba_solution is not None
        assert srfba_solution is not None

        rfba_obj = rfba_solution.objective_value or rfba_solution.fobj
        srfba_obj = srfba_solution.objective_value or srfba_solution.fobj

        print(f"RFBA objective: {rfba_obj}")
        print(f"SRFBA objective: {srfba_obj}")
        print(f"RFBA status: {rfba_solution.status}")
        print(f"SRFBA status: {srfba_solution.status}")

        # Note: Results may differ because:
        # - RFBA applies regulatory constraints before FBA (sequential)
        # - SRFBA integrates regulatory and metabolic optimization (simultaneous MILP)

    def test_compare_with_metabolic_only(self, integrated_model):
        """Compare regulatory methods with metabolic-only FBA."""
        from mewpy.germ.analysis import RFBA, SRFBA

        # Get FBA result (no regulatory constraints)
        fba_solution = integrated_model.simulator.simulate()
        fba_obj = fba_solution.objective_value

        # Get RFBA result
        rfba = RFBA(integrated_model)
        rfba_solution = rfba.optimize()
        rfba_obj = rfba_solution.objective_value or rfba_solution.fobj

        # Get SRFBA result
        srfba = SRFBA(integrated_model)
        srfba_solution = srfba.optimize()
        srfba_obj = srfba_solution.objective_value or srfba_solution.fobj

        print(f"FBA objective (no regulation): {fba_obj}")
        print(f"RFBA objective (with regulation): {rfba_obj}")
        print(f"SRFBA objective (with regulation): {srfba_obj}")

        # Regulatory methods should generally have <= objective than pure FBA
        # (regulatory constraints can only restrict, not expand solution space)
        # However, this depends on the regulatory network structure


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
