"""
Comprehensive tests for RegulatoryExtension class.

Tests factory methods, regulatory network management, and integration with analysis methods.
"""

import sys
from pathlib import Path

import pytest

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class TestRegulatoryExtensionFactoryMethods:
    """Test factory methods for creating RegulatoryExtension instances."""

    def test_from_sbml_metabolic_only(self):
        """Test from_sbml() with metabolic model only."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"

        model = RegulatoryExtension.from_sbml(str(model_path), flavor="reframed")

        assert model is not None
        assert len(model.reactions) > 0
        assert len(model.genes) > 0
        assert not model.has_regulatory_network()

    def test_from_sbml_with_regulatory(self):
        """Test from_sbml() with metabolic + regulatory network."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        model = RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

        assert model is not None
        assert len(model.reactions) > 0
        assert model.has_regulatory_network()
        assert len(model.regulators) > 0
        assert len(model.interactions) > 0

    def test_from_model_cobrapy(self):
        """Test from_model() with COBRApy model."""
        import cobra

        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        cobra_model = cobra.io.read_sbml_model(str(model_path))
        model = RegulatoryExtension.from_model(cobra_model, str(reg_path), regulatory_format="csv", sep=",")

        assert model is not None
        assert len(model.reactions) > 0
        assert model.has_regulatory_network()


class TestRegulatoryExtensionAPI:
    """Test RegulatoryExtension API methods."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        return RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

    def test_yield_interactions_returns_tuples(self, integrated_model):
        """Test that yield_interactions() returns (id, interaction) tuples."""
        interactions = list(integrated_model.yield_interactions())

        assert len(interactions) > 0

        # Check first interaction is a tuple
        first = interactions[0]
        assert isinstance(first, tuple)
        assert len(first) == 2

        int_id, interaction = first
        assert isinstance(int_id, str)
        assert hasattr(interaction, "target")

    def test_yield_regulators_returns_tuples(self, integrated_model):
        """Test that yield_regulators() returns (id, regulator) tuples."""
        regulators = list(integrated_model.yield_regulators())

        assert len(regulators) > 0

        first = regulators[0]
        assert isinstance(first, tuple)
        assert len(first) == 2

        reg_id, regulator = first
        assert isinstance(reg_id, str)

    def test_yield_targets_returns_tuples(self, integrated_model):
        """Test that yield_targets() returns (id, target) tuples."""
        targets = list(integrated_model.yield_targets())

        assert len(targets) > 0

        first = targets[0]
        assert isinstance(first, tuple)
        assert len(first) == 2

    def test_get_parsed_gpr_caching(self, integrated_model):
        """Test that get_parsed_gpr() caches results."""
        # Get a reaction with GPR
        rxn_id = list(integrated_model.reactions)[0]

        # First call
        gpr1 = integrated_model.get_parsed_gpr(rxn_id)

        # Second call should return cached version (same object)
        gpr2 = integrated_model.get_parsed_gpr(rxn_id)

        assert gpr1 is gpr2  # Same object reference = cached


class TestRegulatoryExtensionWithAnalysis:
    """Test RegulatoryExtension with analysis methods."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        from mewpy.germ.models import RegulatoryExtension

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        return RegulatoryExtension.from_sbml(
            str(model_path), str(reg_path), regulatory_format="csv", sep=",", flavor="reframed"
        )

    def test_fba_analysis(self, integrated_model):
        """Test FBA works with RegulatoryExtension via simulator."""
        # Use simulator directly which works with RegulatoryExtension
        result = integrated_model.simulator.simulate()

        assert result.objective_value is not None
        assert result.objective_value > 0

    def test_rfba_analysis(self, integrated_model):
        """Test RFBA works with RegulatoryExtension."""
        from mewpy.germ.analysis import RFBA
        from mewpy.solvers.solution import Status

        rfba = RFBA(integrated_model)
        solution = rfba.optimize()

        assert solution is not None
        # RFBA may return infeasible with default initial state due to regulatory constraints
        # This is correct behavior - we just verify it runs without errors
        assert solution.status in [Status.OPTIMAL, Status.INFEASIBLE]

    def test_srfba_analysis(self, integrated_model):
        """Test SRFBA works with RegulatoryExtension."""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model)
        solution = srfba.optimize()

        assert solution is not None
        assert solution.objective_value is not None or solution.fobj is not None


class TestBackwardsCompatibility:
    """Test backwards compatibility with legacy models."""

    def test_analysis_with_legacy_model(self):
        """Test that analysis methods work with legacy read_model()."""
        from mewpy.io import Engines, Reader, read_model
        from mewpy.simulation import get_simulator

        model_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
        reg_path = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"

        metabolic_reader = Reader(Engines.MetabolicSBML, str(model_path))
        regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, str(reg_path), sep=",")

        legacy_model = read_model(metabolic_reader, regulatory_reader, warnings=False)

        # Should work with simulator
        simulator = get_simulator(legacy_model)
        result = simulator.simulate()

        assert result.objective_value is not None
        assert result.objective_value > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
