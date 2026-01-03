"""Test CoRegFlux workflow from GERM_Models_analysis notebook"""

from pathlib import Path
from mewpy.omics import ExpressionSet
from mewpy.io import Engines, Reader, read_model
from mewpy.germ.analysis.coregflux import predict_gene_expression

print("=" * 80)
print("TESTING COREGFLUX WORKFLOW FROM NOTEBOOK")
print("=" * 80)
print()

# Setup paths
path = Path("examples/models/germ")

# Step 1: Load expression data files
print("Step 1: Load expression data files")
print("-" * 80)
try:
    expression = ExpressionSet.from_csv(
        path.joinpath("iMM904_gene_expression.csv"), sep=";", index_col=0, header=0
    ).dataframe
    print(f"✓ Expression data loaded: {expression.shape}")

    influence = ExpressionSet.from_csv(
        path.joinpath("iMM904_influence.csv"), sep=";", index_col=0, header=0
    ).dataframe
    print(f"✓ Influence data loaded: {influence.shape}")

    experiments = ExpressionSet.from_csv(
        path.joinpath("iMM904_experiments.csv"), sep=";", index_col=0, header=0
    ).dataframe
    print(f"✓ Experiments data loaded: {experiments.shape}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback

    traceback.print_exc()
    exit(1)
print()

# Step 2: Load model
print("Step 2: Load iMM904 model")
print("-" * 80)
try:
    imm904_gem_reader = Reader(Engines.MetabolicSBML, path.joinpath("iMM904.xml"))
    imm904_trn_reader = Reader(
        Engines.CoExpressionRegulatoryCSV,
        path.joinpath("iMM904_trn.csv"),
        sep=",",
        target_col=2,
        co_activating_col=3,
        co_repressing_col=4,
        header=0,
    )

    model = read_model(imm904_gem_reader, imm904_trn_reader)
    print(f"✓ Model loaded: {model.id}")
    print(f"  Reactions: {len(model.reactions)}")
    print(
        f"  Regulators: {len(model.regulators) if hasattr(model, 'regulators') else 'N/A'}"
    )
    print(f"  Targets: {len(model.targets) if hasattr(model, 'targets') else 'N/A'}")
    print(
        f"  Interactions: {len(model.interactions) if hasattr(model, 'interactions') else 'N/A'}"
    )
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback

    traceback.print_exc()
    exit(1)
print()

# Step 3: Predict gene expression
print("Step 3: Predict gene expression using CoRegFlux")
print("-" * 80)
try:
    gene_expression_prediction = predict_gene_expression(
        model=model, influence=influence, expression=expression, experiments=experiments
    )

    print(f"✓ Gene expression prediction completed")
    print(f"  Type: {type(gene_expression_prediction)}")
    if hasattr(gene_expression_prediction, "shape"):
        print(f"  Shape: {gene_expression_prediction.shape}")
    print(f"  Length: {len(gene_expression_prediction)}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback

    traceback.print_exc()
    exit(1)
print()

print("=" * 80)
print("COREGFLUX WORKFLOW TEST COMPLETE")
print("=" * 80)
