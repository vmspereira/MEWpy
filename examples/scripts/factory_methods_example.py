"""
Example: Using RegulatoryExtension Factory Methods

This script demonstrates the convenient factory methods for creating
RegulatoryExtension models from files or existing models.

Note: By default, factory methods use 'reframed' as it is more lightweight
than COBRApy. You can explicitly use 'cobra' by passing flavor='cobra'.
"""

print("=" * 80)
print("REGULATORYEXTENSION FACTORY METHODS EXAMPLE")
print("=" * 80)
print("Note: Uses reframed by default (lightweight)")
print("=" * 80)

# =============================================================================
# Method 1: from_sbml() - The Simplest Way
# =============================================================================
print("\n1. Using from_sbml() - Load from files")
print("-" * 80)

from mewpy.germ.models import RegulatoryExtension

# Load metabolic model only (no regulatory network)
from pathlib import Path
examples_dir = Path(__file__).parent.parent
model_path = examples_dir / 'models' / 'germ' / 'e_coli_core.xml'
regulatory_path = examples_dir / 'models' / 'germ' / 'e_coli_core_trn.csv'

model_metabolic_only = RegulatoryExtension.from_sbml(
    str(model_path)
    # Uses reframed by default (lightweight)
    # Can specify flavor='cobra' if needed
)

print(f"Created metabolic-only model:")
print(f"  - ID: {model_metabolic_only.id}")
print(f"  - Reactions: {len(model_metabolic_only.reactions)}")
print(f"  - Genes: {len(model_metabolic_only.genes)}")
print(f"  - Has regulatory network: {model_metabolic_only.has_regulatory_network()}")

# Load with regulatory network (uses reframed by default)
model_integrated = RegulatoryExtension.from_sbml(
    str(model_path),
    str(regulatory_path),
    regulatory_format='csv',
    sep=','  # Additional CSV parameters passed to reader
    # flavor='reframed' is the default (lightweight)
    # Use flavor='cobra' if you specifically need COBRApy
)

print(f"\nCreated integrated model:")
print(f"  - ID: {model_integrated.id}")
print(f"  - Reactions: {len(model_integrated.reactions)}")
print(f"  - Genes: {len(model_integrated.genes)}")
print(f"  - Has regulatory network: {model_integrated.has_regulatory_network()}")
print(f"  - Interactions: {len(list(model_integrated.yield_interactions()))}")  # Counts (id, interaction) tuples

# =============================================================================
# Method 2: from_model() - From COBRApy/reframed Model
# =============================================================================
print("\n2. Using from_model() - From existing COBRApy model")
print("-" * 80)

import cobra

# Load a COBRApy model first (maybe you already have one)
cobra_model = cobra.io.read_sbml_model(str(model_path))
print(f"Loaded COBRApy model: {cobra_model.id} ({len(cobra_model.reactions)} reactions)")

# Create RegulatoryExtension from it
model_from_cobra = RegulatoryExtension.from_model(
    cobra_model,
    str(regulatory_path),
    regulatory_format='csv',
    sep=','
)

print(f"Created RegulatoryExtension from COBRApy model:")
print(f"  - Reactions: {len(model_from_cobra.reactions)}")
print(f"  - Interactions: {len(list(model_from_cobra.yield_interactions()))}")  # Counts (id, interaction) tuples

# =============================================================================
# Method 3: Using the Models in Analysis
# =============================================================================
print("\n3. Using factory-created models in analysis")
print("-" * 80)

from mewpy.germ.analysis import FBA, RFBA, SRFBA

# FBA on metabolic-only model
fba = FBA(model_metabolic_only)
fba_solution = fba.optimize()
fba_obj = fba_solution.objective_value if fba_solution.objective_value is not None else fba_solution.fobj
print(f"FBA (metabolic only): {fba_obj:.6f}")

# RFBA on integrated model
rfba = RFBA(model_integrated)
rfba_solution = rfba.optimize()
rfba_obj = rfba_solution.objective_value if rfba_solution.objective_value is not None else rfba_solution.fobj
print(f"RFBA (with regulation): {rfba_obj:.6f}")

# SRFBA on integrated model
srfba = SRFBA(model_integrated)
srfba_solution = srfba.optimize()
srfba_obj = srfba_solution.objective_value if srfba_solution.objective_value is not None else srfba_solution.fobj
print(f"SRFBA (steady-state): {srfba_obj:.6f}")

# =============================================================================
# Comparison with Traditional Method
# =============================================================================
print("\n4. Comparison: Factory vs Traditional")
print("-" * 80)

# Traditional method (more verbose)
from mewpy.io import Reader, Engines, read_model
from mewpy.simulation import get_simulator

metabolic_reader = Reader(Engines.MetabolicSBML, str(model_path))
regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                           str(regulatory_path),
                           sep=',')
traditional_model = read_model(metabolic_reader, regulatory_reader)

# Factory method (one line) - uses reframed by default
factory_model = RegulatoryExtension.from_sbml(
    str(model_path),
    str(regulatory_path),
    sep=','
    # Uses reframed by default (lightweight)
)

print(f"Traditional method - Type: {type(traditional_model).__name__}")
print(f"Factory method - Type: {type(factory_model).__name__}")
print(f"\nBoth work with analysis methods:")

rfba_traditional = RFBA(traditional_model)
rfba_factory = RFBA(factory_model)

print(f"  - Traditional: {rfba_traditional.optimize().objective_value:.6f}")
print(f"  - Factory: {rfba_factory.optimize().objective_value:.6f}")

# =============================================================================
# Advanced: Using with Different Backends
# =============================================================================
print("\n5. Advanced: Using with reframed backend")
print("-" * 80)

try:
    # Use reframed instead of cobra
    model_reframed = RegulatoryExtension.from_sbml(
        str(model_path),
        str(regulatory_path),
        flavor='reframed',  # Different backend
        sep=','
    )

    print(f"Created model with reframed backend:")
    print(f"  - Reactions: {len(model_reframed.reactions)}")

    # Still works with analysis methods
    rfba_reframed = RFBA(model_reframed)
    solution = rfba_reframed.optimize()
    if solution and solution.objective_value is not None:
        print(f"  - RFBA result: {solution.objective_value:.6f}")
    else:
        print(f"  - RFBA result: {solution.objective_value}")

except ImportError:
    print("  ⊘ reframed not installed (optional)")
except Exception as e:
    print(f"  ⊘ reframed test skipped: {e}")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)

print("""
Factory methods provide a convenient way to create RegulatoryExtension models:

✓ from_sbml() - Load from SBML + optional CSV regulatory network
  - Simplest method for most use cases
  - Uses 'reframed' by default (lightweight)
  - Supports both 'reframed' (default) and 'cobra' flavors

✓ from_model() - Wrap existing COBRApy/reframed model
  - Useful when you already have a model object
  - Can add regulatory network from file

✓ from_json() - Load complete model from JSON
  - For serialized integrated models

Benefits:
- Less boilerplate code (1-2 lines instead of 5-6)
- Clear, readable API
- Lightweight by default (reframed preferred)
- Flexible arguments for different file formats
- Works seamlessly with all analysis methods

Compare:
  # Old way (6 lines)
  metabolic_reader = Reader(Engines.MetabolicSBML, 'model.xml')
  regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, 'reg.csv', sep=',')
  legacy_model = read_model(metabolic_reader, regulatory_reader)

  # New way (1 line, uses reframed by default)
  model = RegulatoryExtension.from_sbml('model.xml', 'reg.csv', sep=',')

  # Or explicitly use COBRApy if needed
  model = RegulatoryExtension.from_sbml('model.xml', 'reg.csv', sep=',', flavor='cobra')
""")

print("=" * 80)
