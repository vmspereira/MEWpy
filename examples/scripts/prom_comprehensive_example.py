"""
Comprehensive PROM (Probabilistic Regulation of Metabolism) Example

This example demonstrates:
1. Loading an integrated metabolic-regulatory model
2. Computing target-regulator interaction probabilities from gene expression data
3. Performing single and multiple regulator knockout simulations with PROM
4. Comparing PROM results with FBA baseline

PROM uses probabilistic constraints to predict the effect of transcriptional
regulator perturbations on metabolic fluxes.

Reference: Chandrasekaran & Price (2010) https://doi.org/10.1073/pnas.1005139107
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

from mewpy.germ.analysis import PROM, target_regulator_interaction_probability
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


def create_mock_expression_data(model, n_samples=50):
    """
    Create mock gene expression data for demonstration.

    In real applications, use actual microarray or RNA-seq data.
    """
    # Get all genes and regulators
    genes = list(model.targets.keys())[:20]

    # Create random expression matrix (genes x samples)
    expression = pd.DataFrame(
        np.random.randn(len(genes), n_samples),
        index=genes
    )

    return expression


def example_1_basic_prom():
    """Example 1: Basic PROM setup and single regulator knockout."""
    print("=" * 80)
    print("Example 1: Basic PROM with Single Regulator Knockout")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()
    print(f"Loaded model with {len(model.reactions)} reactions, {len(model.genes)} genes")
    print(f"Regulatory network: {len(model.regulators)} regulators, {len(model.targets)} targets")

    # Create PROM instance
    prom = PROM(model)
    prom.build()

    print(f"\nPROM method: {prom.method}")
    print(f"Synchronized with model: {prom.synchronized}")

    # Get first regulator for knockout
    regulators = list(model.regulators.keys())
    if regulators:
        regulator = regulators[0]
        print(f"\nTesting knockout of regulator: {regulator}")

        # Run PROM with default probabilities (1.0 = no effect)
        result = prom.optimize(regulators=[regulator])

        # Get the solution
        sol = result.solutions[f"ko_{regulator}"]
        print(f"  Status: {sol.status}")
        print(f"  Objective value: {sol.objective_value:.4f}")

        # Compare with FBA baseline
        fba_result = model.simulator.simulate()
        print(f"\nFBA baseline objective: {fba_result.objective_value:.4f}")
        print(f"PROM knockout objective: {sol.objective_value:.4f}")
        print(f"Growth reduction: {(1 - sol.objective_value/fba_result.objective_value)*100:.2f}%")
    else:
        print("No regulators found in model")


def example_2_probabilities_from_expression():
    """Example 2: Calculate interaction probabilities from gene expression."""
    print("\n" + "=" * 80)
    print("Example 2: Computing Target-Regulator Interaction Probabilities")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()

    # Create mock expression data (in real use, load actual data)
    expression = create_mock_expression_data(model, n_samples=50)

    # Quantile preprocessing
    print(f"\nExpression data: {expression.shape[0]} genes, {expression.shape[1]} samples")

    # Create binary expression (threshold at median)
    binary_expression = (expression > expression.median(axis=1).values.reshape(-1, 1)).astype(int)

    # Calculate probabilities
    print("Calculating target-regulator interaction probabilities...")
    probabilities, missed = target_regulator_interaction_probability(
        model, expression, binary_expression
    )

    print(f"\nCalculated {len(probabilities)} interaction probabilities")
    print(f"Missed {sum(missed.values())} interactions (no significant correlation)")

    # Show some example probabilities
    print("\nExample interaction probabilities:")
    for i, ((target, regulator), prob) in enumerate(list(probabilities.items())[:5]):
        print(f"  P({target}=1 | {regulator}=0) = {prob:.3f}")

    # Run PROM with these probabilities
    print("\nRunning PROM with calculated probabilities...")
    prom = PROM(model).build()

    # Test first 3 regulators
    test_regulators = list(model.regulators.keys())[:3]
    result = prom.optimize(initial_state=probabilities, regulators=test_regulators)

    print(f"\nTested {len(test_regulators)} regulator knockouts:")
    for reg in test_regulators:
        sol = result.solutions[f"ko_{reg}"]
        print(f"  {reg}: objective = {sol.objective_value:.4f}, status = {sol.status}")


def example_3_multiple_knockouts():
    """Example 3: Multiple regulator knockouts with custom probabilities."""
    print("\n" + "=" * 80)
    print("Example 3: Multiple Regulator Knockouts")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()
    prom = PROM(model).build()

    # Create custom probabilities (reduced flux through some reactions)
    probabilities = {}

    # Get regulators and targets
    regulators = list(model.regulators.keys())[:5]
    targets = list(model.targets.keys())[:10]

    # Set probabilities: lower values = stronger regulatory effect
    for target in targets:
        for regulator in regulators:
            # Random probability between 0.3 and 1.0
            probabilities[(target, regulator)] = np.random.uniform(0.3, 1.0)

    print(f"Testing {len(regulators)} regulators with custom probabilities")
    print(f"Total {len(probabilities)} target-regulator interactions defined")

    # Run PROM for all regulators
    result = prom.optimize(initial_state=probabilities, regulators=regulators)

    # Analyze results
    print("\nRegulator knockout results:")
    objectives = []
    for regulator in regulators:
        sol = result.solutions[f"ko_{regulator}"]
        objectives.append(sol.objective_value)
        print(f"  {regulator}: {sol.objective_value:.4f}")

    # Find most and least impactful regulators
    max_idx = np.argmax(objectives)
    min_idx = np.argmin(objectives)

    print(f"\nMost impactful regulator (lowest growth): {regulators[min_idx]} ({objectives[min_idx]:.4f})")
    print(f"Least impactful regulator (highest growth): {regulators[max_idx]} ({objectives[max_idx]:.4f})")


def example_4_prom_fva_workflow():
    """Example 4: PROM workflow showing FVA computation."""
    print("\n" + "=" * 80)
    print("Example 4: PROM Workflow with FVA")
    print("=" * 80)

    # Load model
    model = load_ecoli_model()
    prom = PROM(model).build()

    print("PROM workflow:")
    print("1. Compute wild-type FBA solution")
    print("2. Compute maximum rates using FVA at 99% of wild-type growth")
    print("3. For each regulator knockout:")
    print("   - Identify affected target genes")
    print("   - Apply probabilistic constraints based on interaction probabilities")
    print("   - Solve FBA with modified constraints")

    # Get a regulator
    regulators = list(model.regulators.keys())[:1]
    if not regulators:
        print("\nNo regulators available")
        return

    regulator = regulators[0]
    print(f"\nDemonstrating with regulator: {regulator}")

    # Get targets of this regulator
    reg_obj = model.regulators[regulator]
    targets = list(reg_obj.yield_targets())
    print(f"Number of targets regulated by {regulator}: {len(targets)}")

    if targets:
        print(f"First 5 targets: {[t.id for t in targets[:5]]}")

    # Run with probabilities
    probabilities = {(t.id, regulator): 0.5 for t in targets}

    print("\nRunning PROM with P(target|regulator KO) = 0.5...")
    result = prom.optimize(initial_state=probabilities, regulators=[regulator])

    sol = result.solutions[f"ko_{regulator}"]
    print(f"Result: objective = {sol.objective_value:.4f}")


def main():
    """Run all examples."""
    print("\n")
    print("*" * 80)
    print("PROM - Probabilistic Regulation of Metabolism")
    print("Comprehensive Examples")
    print("*" * 80)

    try:
        example_1_basic_prom()
        example_2_probabilities_from_expression()
        example_3_multiple_knockouts()
        example_4_prom_fva_workflow()

        print("\n" + "=" * 80)
        print("All examples completed successfully!")
        print("=" * 80)
    except Exception as e:
        print(f"\nError running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
