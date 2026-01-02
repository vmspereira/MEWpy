"""
Comprehensive CoRegFlux Example

This example demonstrates:
1. Loading an integrated metabolic-regulatory model
2. Predicting gene expression using linear regression from influence scores
3. Performing steady-state CoRegFlux simulations
4. Performing dynamic CoRegFlux simulations with metabolite tracking
5. Comparing CoRegFlux results with FBA baseline

CoRegFlux integrates transcriptional regulatory networks and gene expression
to improve phenotype prediction using continuous gene states.

Reference: Trébulle et al. (2017) https://doi.org/10.1186/s12918-017-0507-0
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

from mewpy.germ.analysis import CoRegFlux, predict_gene_expression
from mewpy.germ.models import RegulatoryExtension
from mewpy.omics import ExpressionSet


def load_ecoli_model():
    """Load E. coli core model with regulatory network."""
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    model_path = path.joinpath('models', 'germ', 'e_coli_core.xml')
    reg_path = path.joinpath('models', 'germ', 'e_coli_core_trn.csv')

    model = RegulatoryExtension.from_sbml(
        str(model_path),
        str(reg_path),
        regulatory_format='csv',
        sep=',',
        flavor='reframed'
    )

    return model


def create_mock_data(model, n_samples=50, n_experiments=10):
    """
    Create mock expression and influence data for demonstration.

    In real applications:
    - expression: Training gene expression data (genes x samples)
    - influence: Regulator influence scores from CoRegNet (regulators x samples)
    - experiments: Test condition influence scores (regulators x experiments)
    """
    # Get genes and regulators
    genes = list(model.targets.keys())[:20]
    regulators = list(model.regulators.keys())[:10]

    # Create random data
    expression = pd.DataFrame(
        np.random.randn(len(genes), n_samples),
        index=genes
    )

    influence = pd.DataFrame(
        np.random.randn(len(regulators), n_samples),
        index=regulators
    )

    experiments = pd.DataFrame(
        np.random.randn(len(regulators), n_experiments),
        index=regulators
    )

    return expression, influence, experiments


def example_1_basic_coregflux():
    """Example 1: Basic CoRegFlux steady-state simulation."""
    print("=" * 80)
    print("Example 1: Basic CoRegFlux Steady-State Simulation")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()
    print(f"Loaded model with {len(model.reactions)} reactions, {len(model.genes)} genes")
    print(f"Regulatory network: {len(model.regulators)} regulators, {len(model.targets)} targets")

    # Create CoRegFlux instance
    coregflux = CoRegFlux(model)
    coregflux.build()

    print(f"\nCoRegFlux synchronized: {coregflux.synchronized}")

    # Create initial gene state (all genes fully active)
    genes = list(model.targets.keys())[:15]
    initial_state = {gene: 1.0 for gene in genes}

    print(f"\nRunning CoRegFlux with {len(initial_state)} genes at full expression...")

    # Run CoRegFlux
    result = coregflux.optimize(initial_state=initial_state)

    print(f"  Status: {result.status}")
    print(f"  Objective value: {result.objective_value:.4f}")

    # Compare with FBA baseline
    fba_result = model.simulator.simulate()
    print(f"\nFBA baseline objective: {fba_result.objective_value:.4f}")
    print(f"CoRegFlux objective: {result.objective_value:.4f}")


def example_2_gene_expression_prediction():
    """Example 2: Predict gene expression using linear regression."""
    print("\n" + "=" * 80)
    print("Example 2: Gene Expression Prediction")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()

    # Create mock data
    expression, influence, experiments = create_mock_data(model, n_samples=50, n_experiments=5)

    print("\nTraining data:")
    print(f"  Expression: {expression.shape[0]} genes, {expression.shape[1]} samples")
    print(f"  Influence: {influence.shape[0]} regulators, {influence.shape[1]} samples")
    print("\nTest data:")
    print(f"  Experiments: {experiments.shape[0]} regulators, {experiments.shape[1]} conditions")

    # Predict gene expression for experiments
    print("\nPredicting gene expression using linear regression...")
    predictions = predict_gene_expression(model, influence, expression, experiments)

    print(f"\nPredicted expression for {predictions.shape[0]} genes in {predictions.shape[1]} experiments")
    print(f"\nExpression range: [{predictions.min().min():.2f}, {predictions.max().max():.2f}]")

    # Use predictions in CoRegFlux
    print("\nRunning CoRegFlux with predicted expression...")
    coregflux = CoRegFlux(model).build()

    # Use first experiment
    initial_state = predictions.iloc[:, 0].to_dict()
    result = coregflux.optimize(initial_state=initial_state)

    print(f"  Status: {result.status}")
    print(f"  Objective value: {result.objective_value:.4f}")


def example_3_reduced_gene_expression():
    """Example 3: CoRegFlux with varying gene expression levels."""
    print("\n" + "=" * 80)
    print("Example 3: CoRegFlux with Reduced Gene Expression")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()
    coregflux = CoRegFlux(model).build()

    # Test different expression levels
    genes = list(model.targets.keys())[:15]
    expression_levels = [1.0, 0.8, 0.5, 0.3, 0.1]

    print(f"Testing {len(expression_levels)} different expression levels...")
    print(f"Number of genes: {len(genes)}\n")

    results = []
    for level in expression_levels:
        initial_state = {gene: level for gene in genes}
        result = coregflux.optimize(initial_state=initial_state)
        results.append(result.objective_value)

        print(f"Expression level {level:.1f}: objective = {result.objective_value:.4f}")

    # Analyze trend
    print("\nTrend: As gene expression decreases, growth typically decreases")
    print(f"Reduction from 1.0 to 0.1: {(1 - results[-1]/results[0])*100:.1f}%")


def example_4_dynamic_simulation():
    """Example 4: Dynamic CoRegFlux simulation over time."""
    print("\n" + "=" * 80)
    print("Example 4: Dynamic CoRegFlux Simulation")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()
    coregflux = CoRegFlux(model).build()

    # Create time-varying gene states
    genes = list(model.targets.keys())[:15]
    n_steps = 5

    # Gene expression decreases over time (simulating stress response)
    initial_states = [
        {gene: 1.0 - (0.1 * i) for gene in genes}
        for i in range(n_steps)
    ]

    # Time steps
    time_steps = [0.1 * (i + 1) for i in range(n_steps)]

    print("Running dynamic simulation:")
    print(f"  Number of time steps: {n_steps}")
    print(f"  Time points: {time_steps}")
    print("  Gene expression pattern: decreasing from 1.0 to 0.6")

    # Run dynamic simulation
    result = coregflux.optimize(
        initial_state=initial_states,
        time_steps=time_steps
    )

    print("\nDynamic simulation completed!")
    print(f"Number of solutions: {len(result.solutions)}")

    # Show results for each time point
    print("\nTime-course results:")
    for time_key, sol in result.solutions.items():
        print(f"  {time_key}: objective = {sol.objective_value:.4f}, status = {sol.status}")


def example_5_metabolites_and_biomass():
    """Example 5: CoRegFlux with metabolite concentrations and biomass tracking."""
    print("\n" + "=" * 80)
    print("Example 5: CoRegFlux with Metabolites and Biomass")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()
    coregflux = CoRegFlux(model).build()

    # Create gene state
    genes = list(model.targets.keys())[:15]
    initial_state = {gene: 0.8 for gene in genes}

    # Get external metabolites (those with exchange reactions)
    external_mets = [m for m in model.metabolites if m.endswith('_e')][:5]
    metabolites = {met_id: 10.0 for met_id in external_mets}

    print(f"Tracking {len(metabolites)} external metabolites:")
    for met in external_mets:
        print(f"  {met}: 10.0 mM")

    # Initial biomass
    biomass = 0.1

    print(f"\nInitial biomass: {biomass:.2f}")
    print("Time step: 0.1")

    # Run CoRegFlux
    result = coregflux.optimize(
        initial_state=initial_state,
        metabolites=metabolites,
        biomass=biomass,
        time_steps=0.1
    )

    print("\nResults:")
    print(f"  Status: {result.status}")
    print(f"  Objective value: {result.objective_value:.4f}")

    # Check updated metabolites and biomass (if available in result)
    if hasattr(result, 'metabolites'):
        print(f"  Updated metabolites tracked: {len(result.metabolites)} metabolites")

    if hasattr(result, 'biomass'):
        print(f"  Updated biomass: {result.biomass:.4f}")


def example_6_soft_plus_parameter():
    """Example 6: Effect of soft_plus parameter on constraints."""
    print("\n" + "=" * 80)
    print("Example 6: Soft Plus Parameter Effect")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()
    coregflux = CoRegFlux(model).build()

    # Create gene state
    genes = list(model.targets.keys())[:15]
    initial_state = {gene: 0.7 for gene in genes}

    # Test different soft_plus values
    soft_plus_values = [0, 1, 2, 5]

    print("Testing different soft_plus parameter values...")
    print("soft_plus controls the smoothness of reaction bound constraints\n")

    for sp in soft_plus_values:
        result = coregflux.optimize(
            initial_state=initial_state,
            soft_plus=sp
        )

        print(f"soft_plus={sp}: objective = {result.objective_value:.4f}")

    print("\nHigher soft_plus values typically lead to smoother constraint transitions")


def main():
    """Run all examples."""
    print("\n")
    print("*" * 80)
    print("CoRegFlux - Integration of Transcriptional Regulatory Networks")
    print("and Gene Expression with Metabolic Models")
    print("Comprehensive Examples")
    print("*" * 80)

    try:
        example_1_basic_coregflux()
        example_2_gene_expression_prediction()
        example_3_reduced_gene_expression()
        example_4_dynamic_simulation()
        example_5_metabolites_and_biomass()
        example_6_soft_plus_parameter()

        print("\n" + "=" * 80)
        print("All examples completed successfully!")
        print("=" * 80)
    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
