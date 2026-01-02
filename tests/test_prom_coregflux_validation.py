"""
Comprehensive tests for PROM and CoRegFlux implementations.

Tests verify that:
1. PROM correctly applies probabilistic regulatory constraints
2. CoRegFlux integrates gene expression predictions
3. Results are consistent with expected behavior
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class TestPROMValidation:
    """Validation tests for PROM implementation."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        return RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

    def test_prom_basic_functionality(self, integrated_model):
        """Test basic PROM functionality."""
        from mewpy.germ.analysis import PROM

        # Create PROM instance
        prom = PROM(integrated_model)

        # Test that it builds without errors
        prom.build()

        assert prom.method == "PROM"
        assert prom.synchronized
        print("PROM builds successfully")

    def test_prom_with_probabilities(self, integrated_model):
        """Test PROM with interaction probabilities."""
        from mewpy.germ.analysis import PROM

        prom = PROM(integrated_model).build()

        # Create some sample probabilities (target, regulator): probability
        # Probability of 1.0 means no effect, < 1.0 means reduced flux
        probabilities = {}

        # Get some regulators
        regulators = list(integrated_model.regulators.keys())[:3]
        targets = list(integrated_model.targets.keys())[:5]

        for target in targets:
            for regulator in regulators:
                probabilities[(target, regulator)] = 0.5

        print(f"Testing PROM with {len(probabilities)} interaction probabilities")

        # Run PROM with first regulator knockout
        result = prom.optimize(initial_state=probabilities, regulators=[regulators[0]])

        assert result is not None
        print(f"PROM result type: {type(result).__name__}")
        print(f"Number of solutions: {len(result.solutions) if hasattr(result, 'solutions') else 1}")

    def test_prom_single_regulator_ko(self, integrated_model):
        """Test PROM with single regulator knockout."""
        from mewpy.germ.analysis import PROM

        prom = PROM(integrated_model).build()

        # Get first regulator
        regulators = list(integrated_model.regulators.keys())
        if len(regulators) > 0:
            regulator = regulators[0]

            # Run without probabilities (default to 1.0)
            result = prom.optimize(regulators=[regulator])

            assert result is not None
            assert hasattr(result, "solutions")

            # Get the solution for this regulator
            sol = result.solutions.get(f"ko_{regulator}")
            if sol:
                print(f"Regulator {regulator} knockout:")
                print(f"  Status: {sol.status}")
                print(f"  Objective: {sol.objective_value}")
        else:
            pytest.skip("No regulators found in model")

    def test_prom_multiple_regulator_ko(self, integrated_model):
        """Test PROM with multiple regulator knockouts."""
        from mewpy.germ.analysis import PROM

        prom = PROM(integrated_model).build()

        # Get first 3 regulators
        regulators = list(integrated_model.regulators.keys())[:3]

        if len(regulators) > 0:
            result = prom.optimize(regulators=regulators)

            assert result is not None
            assert hasattr(result, "solutions")
            print(f"Tested {len(regulators)} regulator knockouts")
            print(f"Got {len(result.solutions)} solutions")
        else:
            pytest.skip("Not enough regulators in model")

    def test_prom_probability_calculation(self, integrated_model):
        """Test PROM probability calculation function."""
        pytest.importorskip("scipy", reason="scipy is required for PROM probability calculation")
        from mewpy.germ.analysis import target_regulator_interaction_probability

        # Create mock expression data
        genes = list(integrated_model.targets.keys())[:10]

        # Create random expression matrix (genes x samples)
        n_samples = 20
        expression = pd.DataFrame(np.random.randn(len(genes), n_samples), index=genes)

        # Create binary expression (thresholded)
        binary_expression = (expression > 0).astype(int)

        # Calculate probabilities
        probs, missed = target_regulator_interaction_probability(integrated_model, expression, binary_expression)

        assert isinstance(probs, dict)
        assert isinstance(missed, dict)
        print(f"Calculated {len(probs)} interaction probabilities")
        print(f"Missed {sum(missed.values())} interactions")

        # Check that probabilities are between 0 and 1
        for (target, reg), prob in probs.items():
            assert 0 <= prob <= 1, f"Probability {prob} not in [0, 1] for {target}-{reg}"


class TestCoRegFluxValidation:
    """Validation tests for CoRegFlux implementation."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        return RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

    def test_coregflux_basic_functionality(self, integrated_model):
        """Test basic CoRegFlux functionality."""
        from mewpy.germ.analysis import CoRegFlux

        # Create CoRegFlux instance
        coregflux = CoRegFlux(integrated_model)

        # Test that it builds without errors
        coregflux.build()

        assert coregflux.synchronized
        print("CoRegFlux builds successfully")

    def test_coregflux_with_gene_state(self, integrated_model):
        """Test CoRegFlux with gene state."""
        from mewpy.germ.analysis import CoRegFlux

        coregflux = CoRegFlux(integrated_model).build()

        # Create gene state (all genes active)
        # model.genes is a list of gene IDs
        genes = integrated_model.genes[:10] if hasattr(integrated_model, "genes") and integrated_model.genes else []
        if not genes:
            # Try targets instead
            genes = list(integrated_model.targets.keys())[:10]

        initial_state = {gene: 1.0 for gene in genes}

        print(f"Testing CoRegFlux with {len(initial_state)} genes")

        # Run CoRegFlux
        result = coregflux.optimize(initial_state=initial_state)

        assert result is not None
        assert hasattr(result, "status")
        assert hasattr(result, "objective_value")
        print("CoRegFlux result:")
        print(f"  Status: {result.status}")
        print(f"  Objective: {result.objective_value}")

    def test_coregflux_dynamic_simulation(self, integrated_model):
        """Test CoRegFlux dynamic simulation with multiple time steps."""
        from mewpy.germ.analysis import CoRegFlux

        coregflux = CoRegFlux(integrated_model).build()

        # Create multiple gene states for time steps
        genes = list(integrated_model.targets.keys())[:10]

        # Create 3 time steps with varying gene expression
        initial_states = [
            {gene: 1.0 for gene in genes},  # Time 0: all active
            {gene: 0.8 for gene in genes},  # Time 1: reduced
            {gene: 0.6 for gene in genes},  # Time 2: further reduced
        ]

        time_steps = [0.1, 0.2, 0.3]

        print(f"Testing dynamic CoRegFlux with {len(initial_states)} time steps")

        # Run dynamic simulation
        result = coregflux.optimize(initial_state=initial_states, time_steps=time_steps)

        assert result is not None
        print(f"Dynamic result type: {type(result).__name__}")

        if hasattr(result, "solutions"):
            print(f"Number of time points: {len(result.solutions)}")

    def test_coregflux_gene_expression_prediction(self, integrated_model):
        """Test CoRegFlux gene expression prediction function."""
        pytest.importorskip("sklearn", reason="sklearn is required for CoRegFlux gene expression prediction")
        from mewpy.germ.analysis import predict_gene_expression

        # Create mock data
        targets = list(integrated_model.targets.keys())[:10]
        regulators = list(integrated_model.regulators.keys())[:5]

        # Influence matrix (regulators x samples)
        n_samples = 20
        influence = pd.DataFrame(np.random.randn(len(regulators), n_samples), index=regulators)

        # Expression matrix (targets x samples)
        expression = pd.DataFrame(np.random.randn(len(targets), n_samples), index=targets)

        # Experiments (regulators x test_conditions)
        n_experiments = 5
        experiments = pd.DataFrame(np.random.randn(len(regulators), n_experiments), index=regulators)

        print("Testing gene expression prediction")
        print(f"  Regulators: {len(regulators)}")
        print(f"  Targets: {len(targets)}")
        print(f"  Training samples: {n_samples}")
        print(f"  Test experiments: {n_experiments}")

        # Predict gene expression
        predictions = predict_gene_expression(integrated_model, influence, expression, experiments)

        assert isinstance(predictions, pd.DataFrame)
        print(f"Predicted expression for {predictions.shape[0]} genes in {predictions.shape[1]} experiments")

    def test_coregflux_with_metabolites(self, integrated_model):
        """Test CoRegFlux with metabolite concentrations."""
        from mewpy.germ.analysis import CoRegFlux

        coregflux = CoRegFlux(integrated_model).build()

        # Create gene state
        genes = list(integrated_model.targets.keys())[:10]
        initial_state = {gene: 1.0 for gene in genes}

        # Create metabolite concentrations using external metabolites that have exchange reactions
        # External metabolites typically end with '_e'
        external_mets = [m for m in integrated_model.metabolites if m.endswith("_e")][:5]
        metabolites = {met_id: 1.0 for met_id in external_mets}

        print("Testing CoRegFlux with metabolites")

        # Run CoRegFlux with metabolites
        result = coregflux.optimize(initial_state=initial_state, metabolites=metabolites)

        assert result is not None
        print(f"  Status: {result.status}")
        print(f"  Objective: {result.objective_value}")


class TestPROMvsCoRegFlux:
    """Compare PROM and CoRegFlux results."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        return RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

    def test_compare_with_fba(self, integrated_model):
        """Compare PROM and CoRegFlux with pure FBA."""
        from mewpy.germ.analysis import PROM, CoRegFlux

        # Get FBA baseline
        fba_result = integrated_model.simulator.simulate()
        fba_obj = fba_result.objective_value

        # Test PROM (with single regulator)
        prom = PROM(integrated_model).build()
        regulators = list(integrated_model.regulators.keys())[:1]
        if regulators:
            prom_result = prom.optimize(regulators=regulators)
            prom_sol = list(prom_result.solutions.values())[0]
            prom_obj = prom_sol.objective_value

            print(f"FBA objective: {fba_obj}")
            print(f"PROM objective (1 regulator KO): {prom_obj}")

        # Test CoRegFlux
        coregflux = CoRegFlux(integrated_model).build()
        genes = list(integrated_model.targets.keys())[:10]
        initial_state = {gene: 1.0 for gene in genes}
        coregflux_result = coregflux.optimize(initial_state=initial_state)
        coregflux_obj = coregflux_result.objective_value

        print(f"CoRegFlux objective: {coregflux_obj}")

        # All should return valid objectives
        assert fba_obj > 0
        assert coregflux_obj >= 0


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
