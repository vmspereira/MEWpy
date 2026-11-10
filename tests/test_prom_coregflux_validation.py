"""
Comprehensive tests for PROM and CoRegFlux implementations.

Tests verify that:
1. PROM correctly applies probabilistic regulatory constraints
2. CoRegFlux integrates gene expression predictions
3. Results are consistent with expected behavior

The expected growth rates are computed with COBRApy on the same SBML file, applying the bounds the
methods are defined to set (PROM: probability * maximum flux; CoRegFlux: softplus of the continuous GPR).
"""

import ast
import math
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

MODEL_PATH = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
REG_PATH = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"
# e_coli_core_trn.csv: header row, columns bNum, Gene (alias), Rule
TRN_CSV = dict(regulatory_format="csv", sep=",", id_col=0, rule_col=2, aliases_cols=[1], header=0)
# Orth et al. (2010): e_coli_core grows at 0.8739 /h on 10 mmol/gDW/h glucose, aerobic
AEROBIC_GROWTH = 0.8739
# CoRegFlux defaults: time step and the TOLERANCE below which a bound or a state counts as zero
TIME_STEP, TOLERANCE = 0.1, 1e-10


def _integrated_model():
    from mewpy.germ.models import RegulatoryExtension

    return RegulatoryExtension.from_sbml(str(MODEL_PATH), str(REG_PATH), flavor="reframed", **TRN_CSV)


def _cobra_id(identifier):
    """The COBRApy id of a REFRAMED reaction, gene or metabolite id"""
    return identifier[2:] if identifier[:2] in ("R_", "G_", "M_") else identifier


@lru_cache(maxsize=None)
def _cobra_model():
    import cobra

    return cobra.io.read_sbml_model(str(MODEL_PATH))


def _wild_type():
    return _cobra_model().slim_optimize()


def _prom_growth(probabilities, regulator, model):
    """PROM knockout of a regulator: reactions that its metabolic targets switch off are limited to
    probability * maximum flux (FVA at 99% growth), in the direction of the wild-type flux"""
    from cobra.flux_analysis import flux_variability_analysis

    cobra_model = _cobra_model()
    targets = [t.id for t in model.get_regulator(regulator).yield_targets() if t.id in model.genes]
    knockouts = {_cobra_id(t) for t in targets}
    limits = {}
    for target in targets:
        p = probabilities.get((target, regulator))
        if p is None or p >= 1:
            continue
        for rxn in cobra_model.genes.get_by_id(_cobra_id(target)).reactions:
            if not rxn.gpr.eval(knockouts=knockouts):
                limits[rxn.id] = min(p, limits.get(rxn.id, 1))
    if not limits:
        return _wild_type()
    reference = cobra_model.optimize().fluxes
    fva = flux_variability_analysis(cobra_model, reaction_list=list(limits), fraction_of_optimum=0.99)
    with cobra_model:
        for r_id, p in limits.items():
            rxn = cobra_model.reactions.get_by_id(r_id)
            if reference[r_id] > 0:
                rxn.upper_bound = max(min(rxn.upper_bound, p * fva.loc[r_id, "maximum"]), TOLERANCE)
            elif reference[r_id] < 0:
                rxn.lower_bound = min(max(rxn.lower_bound, p * fva.loc[r_id, "minimum"]), -TOLERANCE)
        return cobra_model.slim_optimize()


def _gpr_level(node, state):
    """Continuous GPR: min over AND, max over OR"""
    if isinstance(node, ast.Name):
        return state[node.id]
    values = [_gpr_level(v, state) for v in node.values]
    return min(values) if isinstance(node.op, ast.And) else max(values)


def _softplus(x):
    return x + math.log1p(math.exp(-x)) if x > 0 else math.log1p(math.exp(x))


def _coregflux_growth(state, metabolites=None):
    """CoRegFlux step: exchange uptake limited to concentration / (biomass * time step), with the wild-type
    growth as biomass; reactions whose genes all have a state get bounds +-softplus(continuous GPR)"""
    cobra_model = _cobra_model()
    state = {_cobra_id(g): v for g, v in state.items()}
    biomass = _wild_type()
    with cobra_model:
        for met_id, concentration in (metabolites or {}).items():
            rxn = next(r for r in cobra_model.metabolites.get_by_id(_cobra_id(met_id)).reactions if r.boundary)
            uptake = concentration / (biomass * TIME_STEP)
            uptake = 0 if uptake < TOLERANCE else uptake
            if rxn.lower_bound == 0 or abs(rxn.lower_bound) >= uptake:
                rxn.lower_bound = -uptake
        for rxn in cobra_model.reactions:
            if rxn.gpr.body is None or not {g.id for g in rxn.genes} <= set(state):
                continue
            level = _gpr_level(rxn.gpr.body, state)
            if level > TOLERANCE:
                lb = -_softplus(abs(level)) if abs(rxn.lower_bound) > TOLERANCE else rxn.lower_bound
                ub = _softplus(level) if rxn.upper_bound > TOLERANCE else rxn.upper_bound
                rxn.bounds = (lb, ub)
        return cobra_model.slim_optimize()


class TestPROMValidation:
    """Validation tests for PROM implementation."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        return _integrated_model()

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
        """PROM with interaction probabilities of 0.5: a knockout limits the reactions its targets switch off
        to half of their maximum flux; regulators without metabolic targets leave growth unchanged"""
        from mewpy.germ.analysis import PROM

        prom = PROM(integrated_model).build()

        # Create some sample probabilities (target, regulator): probability
        # Probability of 1.0 means no effect, < 1.0 means reduced flux
        regulators = list(integrated_model.regulators.keys())[:3]
        targets = list(integrated_model.targets.keys())[:5]
        probabilities = {(target, regulator): 0.5 for target in targets for regulator in regulators}

        result = prom.optimize(initial_state=probabilities, regulators=regulators)

        growth = {reg: result.solutions[f"ko_{reg}"].objective_value for reg in regulators}
        for regulator in regulators:
            expected = _prom_growth(probabilities, regulator, integrated_model)
            assert growth[regulator] == pytest.approx(expected, abs=1e-6), regulator
        # b0113 (PdhR) represses the pyruvate dehydrogenase genes b0114 and b0115
        assert "b0113" in regulators
        assert growth["b0113"] < _wild_type() - 1e-3

    def test_prom_single_regulator_ko(self, integrated_model):
        """Without interaction probabilities (1.0) a regulator knockout does not change growth"""
        from mewpy.germ.analysis import PROM
        from mewpy.solvers.solution import Status

        prom = PROM(integrated_model).build()

        regulator = list(integrated_model.regulators.keys())[0]
        result = prom.optimize(regulators=[regulator])

        assert list(result.solutions) == [f"ko_{regulator}"]
        sol = result.solutions[f"ko_{regulator}"]
        assert sol.status == Status.OPTIMAL
        assert sol.objective_value == pytest.approx(_wild_type(), abs=1e-6)

    def test_prom_multiple_regulator_ko(self, integrated_model):
        """One knockout solution per regulator; without probabilities none changes growth"""
        from mewpy.germ.analysis import PROM

        prom = PROM(integrated_model).build()

        regulators = list(integrated_model.regulators.keys())[:3]
        result = prom.optimize(regulators=regulators)

        assert set(result.solutions) == {f"ko_{reg}" for reg in regulators}
        for sol in result.solutions.values():
            assert sol.objective_value == pytest.approx(_wild_type(), abs=1e-6)

    def test_prom_knockout_lowers_growth(self):
        """A regulator knockout with interaction probability 0 restricts its targets' reactions."""
        from mewpy.germ.analysis import PROM
        from mewpy.germ.models import RegulatoryExtension

        data = Path(__file__).parent / "data"
        model = RegulatoryExtension.from_sbml(
            str(data / "e_coli_core.xml"),
            str(data / "e_coli_core_trn.csv"),
            regulatory_format="csv",
            sep=",",
            id_col=0,
            rule_col=2,
            aliases_cols=[1],
            header=0,
        )
        probabilities = {}
        for _, interaction in model.yield_interactions():
            for regulator in interaction.yield_regulators():
                probabilities[(interaction.target.id, regulator.id)] = 0.0

        arca = next(regulator for regulator in model.regulators if regulator.endswith("b4401"))
        wild_type = model.simulator.simulate().objective_value
        result = PROM(model).build().optimize(initial_state=probabilities, regulators=[arca])
        growth = next(iter(result.solutions.values())).objective_value
        assert growth < wild_type - 1e-3

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
        return _integrated_model()

    def test_coregflux_basic_functionality(self, integrated_model):
        """CoRegFlux builds; without gene states or metabolites it is FBA"""
        from mewpy.germ.analysis import CoRegFlux

        # Create CoRegFlux instance
        coregflux = CoRegFlux(integrated_model)

        # Test that it builds without errors
        coregflux.build()

        assert coregflux.synchronized
        assert coregflux.optimize().objective_value == pytest.approx(_wild_type(), abs=1e-6)
        assert _wild_type() == pytest.approx(AEROBIC_GROWTH, abs=5e-5)

    def test_coregflux_with_gene_state(self, integrated_model):
        """Gene states bound the reactions of those genes by the softplus of their continuous GPR"""
        from mewpy.germ.analysis import CoRegFlux
        from mewpy.solvers.solution import Status

        coregflux = CoRegFlux(integrated_model).build()

        # Create gene state (all genes active)
        # model.genes is a list of gene IDs
        genes = integrated_model.genes[:10]
        initial_state = {gene: 1.0 for gene in genes}

        result = coregflux.optimize(initial_state=initial_state)

        assert result.status == Status.OPTIMAL
        expected = _coregflux_growth(initial_state)
        assert expected < _wild_type() - 1e-3
        assert result.objective_value == pytest.approx(expected, abs=1e-6)

    def test_coregflux_gene_state_lowers_growth(self, integrated_model):
        """Gene expression bounds the reactions: low predicted expression must lower growth below FBA."""
        from mewpy.germ.analysis import CoRegFlux

        fba = integrated_model.simulator.simulate().objective_value
        state = {gene: 10.0 for gene in integrated_model.genes}
        result = CoRegFlux(integrated_model).build().optimize(initial_state=state)
        assert 0 < result.objective_value < fba - 1e-3

    def test_coregflux_dynamic_simulation(self, integrated_model):
        """One solution per time step, each the CoRegFlux step of that step's gene state"""
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

        result = coregflux.optimize(initial_state=initial_states, time_steps=time_steps)

        solutions = list(result.solutions.values())
        assert len(solutions) == len(time_steps)
        growth = [sol.objective_value for sol in solutions]
        for state, value in zip(initial_states, growth):
            assert value == pytest.approx(_coregflux_growth(state), abs=1e-6)
        # lower expression, tighter bounds
        assert growth[0] > growth[1] > growth[2]

    def test_coregflux_gene_expression_prediction(self, integrated_model):
        """Predicted expression of a target: linear regression on its regulators' influence (least squares)"""
        pytest.importorskip("sklearn", reason="sklearn is required for CoRegFlux gene expression prediction")
        from mewpy.germ.analysis import predict_gene_expression

        rng = np.random.default_rng(0)
        targets = list(integrated_model.targets.keys())[:10]
        regulators = list(integrated_model.regulators.keys())[:5]

        # Influence matrix (regulators x samples)
        n_samples = 20
        influence = pd.DataFrame(rng.standard_normal((len(regulators), n_samples)), index=regulators)

        # Expression matrix (targets x samples)
        expression = pd.DataFrame(rng.standard_normal((len(targets), n_samples)), index=targets)

        # Experiments (regulators x test_conditions)
        n_experiments = 5
        experiments = pd.DataFrame(rng.standard_normal((len(regulators), n_experiments)), index=regulators)

        predictions = predict_gene_expression(integrated_model, influence, expression, experiments)

        target_regulators = {t: [r.id for r in integrated_model.get_target(t).yield_regulators()] for t in targets}
        predictable = {t for t, regs in target_regulators.items() if regs and set(regs) <= set(regulators)}
        assert predictable
        assert set(predictions.index) == predictable
        assert list(predictions.columns) == list(experiments.columns)
        for target in predictable:
            x = influence.loc[target_regulators[target]].T.to_numpy()
            coefficients, *_ = np.linalg.lstsq(
                np.column_stack([x, np.ones(n_samples)]), expression.loc[target].to_numpy()
            )
            for column in experiments.columns:
                point = np.append(experiments.loc[target_regulators[target], column].to_numpy(), 1.0)
                assert predictions.loc[target, column] == pytest.approx(point @ coefficients, abs=1e-8)

    def test_coregflux_with_metabolites(self, integrated_model):
        """Metabolite concentrations bound uptake by concentration / (biomass * time step): extracellular
        acetate, acetaldehyde, 2-oxoglutarate and ethanol become carbon sources besides glucose"""
        from mewpy.germ.analysis import CoRegFlux

        coregflux = CoRegFlux(integrated_model).build()

        # Create gene state
        genes = list(integrated_model.targets.keys())[:10]
        initial_state = {gene: 1.0 for gene in genes}

        # Create metabolite concentrations using external metabolites that have exchange reactions
        # External metabolites typically end with '_e'
        external_mets = [m for m in integrated_model.metabolites if m.endswith("_e")][:5]
        metabolites = {met_id: 1.0 for met_id in external_mets}

        result = coregflux.optimize(initial_state=initial_state, metabolites=metabolites)

        expected = _coregflux_growth(initial_state, metabolites)
        assert result.objective_value == pytest.approx(expected, abs=1e-6)
        assert expected > _coregflux_growth(initial_state) + 1e-3


class TestPROMvsCoRegFlux:
    """Compare PROM and CoRegFlux results."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        return _integrated_model()

    def test_compare_with_fba(self, integrated_model):
        """FBA is the textbook growth; PROM without probabilities equals it; CoRegFlux gene states lower it"""
        from mewpy.germ.analysis import PROM, CoRegFlux

        # Get FBA baseline
        fba_result = integrated_model.simulator.simulate()
        fba_obj = fba_result.objective_value

        # Test PROM (with single regulator)
        prom = PROM(integrated_model).build()
        regulators = list(integrated_model.regulators.keys())[:1]
        prom_result = prom.optimize(regulators=regulators)
        prom_obj = list(prom_result.solutions.values())[0].objective_value

        # Test CoRegFlux
        coregflux = CoRegFlux(integrated_model).build()
        genes = list(integrated_model.targets.keys())[:10]
        initial_state = {gene: 1.0 for gene in genes}
        coregflux_obj = coregflux.optimize(initial_state=initial_state).objective_value

        assert fba_obj == pytest.approx(_wild_type(), abs=1e-6)
        assert fba_obj == pytest.approx(AEROBIC_GROWTH, abs=5e-5)
        assert prom_obj == pytest.approx(fba_obj, abs=1e-6)
        assert coregflux_obj == pytest.approx(_coregflux_growth(initial_state), abs=1e-6)
        assert coregflux_obj < fba_obj - 1e-3


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
