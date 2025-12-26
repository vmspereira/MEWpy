# GERM Migration Guide

## Overview

This guide helps you migrate from the old GERM implementation to the new clean architecture introduced in MEWpy v1.0+.

**Key Changes:**
- New `RegulatoryExtension` class replaces internal metabolic storage
- Simplified factory methods for model creation
- Reframed preferred as default backend (lightweight)
- Consistent API for `yield_*` methods
- Scientific correctness improvements

---

## Quick Migration

### Before (Old Way)
```python
from mewpy.io import Reader, Engines, read_model
from mewpy.simulation import get_simulator

# 6+ lines to create integrated model
metabolic_reader = Reader(Engines.MetabolicSBML, 'model.xml')
regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, 'regulatory.csv', sep=',')
model = read_model(metabolic_reader, regulatory_reader)
```

### After (New Way)
```python
from mewpy.germ.models import RegulatoryExtension

# 1 line to create integrated model (uses reframed by default)
model = RegulatoryExtension.from_sbml('model.xml', 'regulatory.csv', sep=',')

# Or explicitly use COBRApy if needed
model = RegulatoryExtension.from_sbml('model.xml', 'regulatory.csv', sep=',', flavor='cobra')
```

---

## Detailed Migration Steps

### 1. Update Imports

**Before:**
```python
from mewpy.io import Reader, Engines, read_model
from mewpy.germ.models import MetabolicModel, RegulatoryModel
```

**After:**
```python
from mewpy.germ.models import RegulatoryExtension
# That's it! Much simpler.
```

### 2. Update Model Creation

#### Option A: From SBML Files

**Before:**
```python
metabolic_reader = Reader(Engines.MetabolicSBML, 'ecoli_core.xml')
regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, 'ecoli_trn.csv', sep=',')
model = read_model(metabolic_reader, regulatory_reader)
```

**After:**
```python
model = RegulatoryExtension.from_sbml(
    'ecoli_core.xml',
    'ecoli_trn.csv',
    regulatory_format='csv',
    sep=','
)
```

#### Option B: From Existing COBRApy Model

**Before:**
```python
import cobra
from mewpy.simulation import get_simulator

cobra_model = cobra.io.read_sbml_model('model.xml')
simulator = get_simulator(cobra_model)
# Then manually integrate with regulatory network...
```

**After:**
```python
import cobra
from mewpy.germ.models import RegulatoryExtension

cobra_model = cobra.io.read_sbml_model('model.xml')
model = RegulatoryExtension.from_model(
    cobra_model,
    'regulatory.csv',
    regulatory_format='csv',
    sep=','
)
```

### 3. Update Analysis Code

Analysis methods work the same way, but initialization is simpler:

**Before:**
```python
from mewpy.germ.analysis import RFBA

rfba = RFBA(model, solver='cplex', build=True, attach=True)
solution = rfba.optimize()
```

**After:**
```python
from mewpy.germ.analysis import RFBA

rfba = RFBA(model, solver='cplex', build=True)  # attach is unused now
solution = rfba.optimize()
```

**Note:** The `attach` parameter is kept for backwards compatibility but is not used in the new architecture.

### 4. Update yield_interactions() Usage

The API changed to return tuples for consistency:

**Before:**
```python
for interaction in model.yield_interactions():
    print(interaction.target.id)
```

**After:**
```python
for int_id, interaction in model.yield_interactions():
    print(interaction.target.id)

# Or if you don't need the ID:
for _, interaction in model.yield_interactions():
    print(interaction.target.id)
```

---

## Breaking Changes

### 1. yield_interactions() Return Type

**What changed:** `yield_interactions()` now returns `(id, interaction)` tuples instead of just `interaction` objects.

**Why:** API consistency with `yield_regulators()` and `yield_targets()` which already returned tuples.

**Migration:**
```python
# Old code (breaks):
for interaction in model.yield_interactions():
    process(interaction)

# New code (fixed):
for _, interaction in model.yield_interactions():
    process(interaction)
```

### 2. Default Backend Changed to Reframed

**What changed:** Factory methods now default to `flavor='reframed'` instead of `flavor='cobra'`.

**Why:** Reframed is more lightweight and faster.

**Migration:** If you specifically need COBRApy, explicitly set it:
```python
model = RegulatoryExtension.from_sbml('model.xml', flavor='cobra')
```

### 3. Deprecated Classes

The following classes are deprecated but still work for backwards compatibility:

- `MetabolicModel` → Use `RegulatoryExtension.from_sbml()`
- `SimulatorBasedMetabolicModel` → Use `RegulatoryExtension.from_model()`

**Migration:** Update to use factory methods as shown above.

---

## Backwards Compatibility

### Legacy Models Still Work

Old code using `read_model()` continues to work:

```python
from mewpy.io import Reader, Engines, read_model

metabolic_reader = Reader(Engines.MetabolicSBML, 'model.xml')
regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, 'regulatory.csv', sep=',')
legacy_model = read_model(metabolic_reader, regulatory_reader)

# Analysis methods work with both old and new models
from mewpy.germ.analysis import RFBA
rfba = RFBA(legacy_model)  # Still works!
solution = rfba.optimize()
```

The new architecture maintains full backwards compatibility through helper methods in the analysis base class.

---

## Scientific Correctness Improvements

Several scientific issues were fixed in the new implementation:

### 1. pFBA Infeasible Handling (Fixed)

**Before:** Infeasible pFBA solutions were masked as "optimal" with zero fluxes.

**After:** Infeasible status is correctly propagated, allowing proper debugging.

**Impact:** You may now see `Status.INFEASIBLE` where you saw zero solutions before. This is correct - it means your constraints are contradictory.

### 2. RFBA Dynamic State Update (Documented)

**Before:** State update logic was undocumented.

**After:** Clearly documented as a heuristic (regulator active if |flux| > tolerance).

**Impact:** You can now override `_update_state_from_solution()` if you need custom state update logic.

### 3. SRFBA Exception Handling (Logged)

**Before:** Constraint building failures were silently ignored.

**After:** Failures are logged as warnings.

**Impact:** You'll see warnings if GPR linearization or interaction constraints fail, helping you debug issues.

---

## Performance Improvements

### 1. No Data Duplication

**Before:** Metabolic data was duplicated between external simulator and GERM.

**After:** Single source of truth in external simulator, GERM only stores regulatory network.

**Impact:** Lower memory usage, especially for large models.

### 2. Cached GPR Parsing

**Before:** GPRs parsed on every call.

**After:** Parsed GPRs cached in `RegulatoryExtension`.

**Impact:** Faster repeated access to GPR rules.

---

## Examples

### Complete Example: Migrating RFBA Analysis

**Before (Old Code):**
```python
from mewpy.io import Reader, Engines, read_model
from mewpy.germ.analysis import RFBA

# Create model (6 lines)
metabolic_reader = Reader(Engines.MetabolicSBML, 'ecoli_core.xml')
regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, 'ecoli_trn.csv', sep=',')
model = read_model(metabolic_reader, regulatory_reader)

# Run RFBA
rfba = RFBA(model, build=True, attach=True)
solution = rfba.optimize()

# Access results
print(f"Objective: {solution.objective_value}")
```

**After (New Code):**
```python
from mewpy.germ.models import RegulatoryExtension
from mewpy.germ.analysis import RFBA

# Create model (1 line - uses reframed by default)
model = RegulatoryExtension.from_sbml('ecoli_core.xml', 'ecoli_trn.csv', sep=',')

# Run RFBA (attach is unused but kept for compatibility)
rfba = RFBA(model, build=True)
solution = rfba.optimize()

# Access results (same as before)
obj_val = solution.objective_value if solution.objective_value is not None else solution.fobj
print(f"Objective: {obj_val}")
```

---

## Troubleshooting

### Issue: "cannot unpack non-iterable Interaction object"

**Cause:** Code expects old `yield_interactions()` format.

**Fix:** Update loop to unpack tuples:
```python
# Change this:
for interaction in model.yield_interactions():
    ...

# To this:
for _, interaction in model.yield_interactions():
    ...
```

### Issue: "ModuleNotFoundError: No module named 'reframed'"

**Cause:** Reframed is now the default but not installed.

**Fix:** Either install reframed:
```bash
pip install reframed
```

Or explicitly use COBRApy:
```python
model = RegulatoryExtension.from_sbml('model.xml', flavor='cobra')
```

### Issue: "Infeasible solution where I had results before"

**Cause:** pFBA now correctly reports infeasible status instead of masking it.

**Fix:** Check your model constraints - infeasibility indicates contradictory constraints that need fixing.

---

## Testing Your Migration

After migrating, verify your code works:

```python
# Test model creation
model = RegulatoryExtension.from_sbml('model.xml', 'regulatory.csv', sep=',')
assert model.has_regulatory_network()
assert len(model.reactions) > 0

# Test analysis
from mewpy.germ.analysis import RFBA
rfba = RFBA(model)
solution = rfba.optimize()
assert solution.objective_value > 0  # or solution.fobj > 0

# Test yield_interactions()
interactions = list(model.yield_interactions())
assert len(interactions) > 0
assert isinstance(interactions[0], tuple)
print("✓ Migration successful!")
```

---

## Need Help?

- **Documentation:** See `docs/germ.md` for full API reference
- **Examples:** Check `examples/scripts/factory_methods_example.py`
- **Issues:** Report problems at https://github.com/vmspereira/MEWpy/issues

---

## Summary Checklist

- [ ] Updated imports to use `RegulatoryExtension`
- [ ] Replaced `read_model()` with factory methods
- [ ] Updated `yield_interactions()` loops to unpack tuples
- [ ] Tested with your models and analysis workflows
- [ ] Verified infeasible solutions are handled correctly
- [ ] Updated any custom analysis methods if needed

**Estimated migration time:** 15-30 minutes for typical projects.
