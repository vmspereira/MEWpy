"""Test PROM workflow from GERM_Models_analysis notebook"""

from pathlib import Path
from mewpy.omics import ExpressionSet
from mewpy.io import Engines, Reader, read_model
from mewpy.germ.analysis.prom import target_regulator_interaction_probability

print("=" * 80)
print("TESTING PROM WORKFLOW FROM NOTEBOOK")
print("=" * 80)
print()

# Setup paths
path = Path("examples/models/germ")
inj661_gene_expression_path = path.joinpath("iNJ661_gene_expression.csv")

# Step 1: Load expression data
print("Step 1: Load expression data")
print("-" * 80)
try:
    expression = ExpressionSet.from_csv(file_path=inj661_gene_expression_path, sep=";", index_col=0, header=None)
    print(f"✓ Expression data loaded: {expression.shape()}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback

    traceback.print_exc()
    exit(1)
print()

# Step 2: Run quantile pipeline
print("Step 2: Run quantile pipeline")
print("-" * 80)
try:
    quantile_expression, binary_expression = expression.quantile_pipeline()
    print(f"✓ Quantile pipeline completed")
    print(f"  Quantile expression type: {type(quantile_expression)}")
    print(f"  Quantile expression shape: {quantile_expression.shape}")
    print(f"  Binary expression type: {type(binary_expression)}")
    print(f"  Binary expression shape: {binary_expression.shape}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback

    traceback.print_exc()
    exit(1)
print()

# Step 3: Load model
print("Step 3: Load iNJ661 model")
print("-" * 80)
try:
    inj661_gem_reader = Reader(Engines.MetabolicSBML, path.joinpath("iNJ661.xml"))
    inj661_trn_reader = Reader(
        Engines.TargetRegulatorRegulatoryCSV,
        path.joinpath("iNJ661_trn.csv"),
        sep=";",
        target_col=0,
        regulator_col=1,
        header=None,
    )

    model = read_model(inj661_gem_reader, inj661_trn_reader)
    BIOMASS_ID = "biomass_Mtb_9_60atp_test_NOF"
    model.objective = {BIOMASS_ID: 1}

    print(f"✓ Model loaded: {model.id}")
    print(f"  Reactions: {len(model.reactions)}")
    print(f"  Regulators: {len(model.regulators) if hasattr(model, 'regulators') else 'N/A'}")
    print(f"  Targets: {len(model.targets) if hasattr(model, 'targets') else 'N/A'}")
    print(f"  Interactions: {len(model.interactions) if hasattr(model, 'interactions') else 'N/A'}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback

    traceback.print_exc()
    exit(1)
print()

# Step 4: Compute target-regulator interaction probabilities
print("Step 4: Compute target-regulator interaction probabilities")
print("-" * 80)
try:
    initial_state, _ = target_regulator_interaction_probability(
        model, expression=quantile_expression, binary_expression=binary_expression
    )

    print(f"✓ PROM initial state computed")
    print(f"  Type: {type(initial_state)}")
    print(f"  Number of interactions: {len(initial_state)}")
    if initial_state:
        first_items = list(initial_state.items())[:3]
        print(f"  First 3 interactions:")
        for key, value in first_items:
            print(f"    {key}: {value:.4f}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback

    traceback.print_exc()
    exit(1)
print()

# Step 5: Run PROM optimization
print("Step 5: Run PROM optimization")
print("-" * 80)
try:
    from mewpy.germ.analysis import PROM

    prom = PROM(model).build()
    solution = prom.optimize(initial_state=initial_state)

    print(f"✓ PROM optimization completed")
    print(f"  Type: {type(solution)}")
    print(f"  Number of KO solutions: {len(solution.solutions)}")

    # Show first solution
    if solution.solutions:
        first_ko = list(solution.solutions.keys())[0]
        first_sol = solution.solutions[first_ko]
        print(f"  First KO ({first_ko}): {first_sol.objective_value:.6f}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback

    traceback.print_exc()
    exit(1)
print()

print("=" * 80)
print("PROM WORKFLOW TEST COMPLETE")
print("=" * 80)
