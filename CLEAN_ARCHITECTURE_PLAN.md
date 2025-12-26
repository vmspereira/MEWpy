# Clean GERM Architecture - No Backwards Compatibility

**Goal:** Focus on the best implementation for integrating regulatory networks with cobrapy/reframed models

---

## Core Principles

1. **Single Source of Truth:** Metabolic data lives only in cobrapy/reframed models
2. **Clean Separation:** Regulatory network is separate from metabolic model
3. **Pure Delegation:** All metabolic operations delegate to external simulator
4. **No Legacy Support:** Remove all backwards compatibility code
5. **Simple & Clean:** Minimal complexity, maximum clarity

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    User Entry Point                          │
│  cobra.Model OR reframed.CBModel (metabolic model)          │
└────────────────────┬────────────────────────────────────────┘
                     │
                     ▼
           ┌──────────────────┐
           │    Simulator     │  (MEWpy's unified interface)
           └────────┬─────────┘
                    │
                    ▼
      ┌─────────────────────────────┐
      │  RegulatoryExtension        │
      │  - Wraps simulator          │
      │  - Stores regulatory only   │
      │  - No metabolic duplication │
      └─────────────┬───────────────┘
                    │
                    ▼
      ┌─────────────────────────────┐
      │  Analysis Methods           │
      │  - RFBA, SRFBA, PROM        │
      │  - All use RegulatoryExt    │
      │  - No legacy paths          │
      └─────────────────────────────┘
```

---

## Key Changes Required

### 1. Analysis Methods (RFBA, SRFBA, PROM, CoRegFlux)

**Current Issues:**
- Conditional logic checking for legacy models
- Dual code paths (extension vs legacy)
- Complex isinstance() checks
- Backwards compatibility bloat

**Clean Solution:**
```python
class RFBA(FBA):
    """RFBA using RegulatoryExtension only."""

    def __init__(self, model: RegulatoryExtension, solver=None, build=False, attach=False):
        """Only accepts RegulatoryExtension."""
        super().__init__(model=model, solver=solver, build=build, attach=attach)

    def build(self):
        """Simplified build - no backwards compatibility."""
        simulator = self.model.simulator
        self._solver = solver_instance(simulator)
        self._linear_objective = dict(self.model.objective)  # Always string keys
        self._minimize = False
        self._synchronized = True
        return self

    def decode_constraints(self, state):
        """Simplified - only RegulatoryExtension path."""
        constraints = {}
        for rxn_id in self.model.reactions:
            gpr = self.model.get_parsed_gpr(rxn_id)
            if not gpr.is_none and not gpr.evaluate(values=state):
                constraints[rxn_id] = (0.0, 0.0)
        return constraints

    def decode_regulatory_state(self, state):
        """Simplified - only regulatory network path."""
        if not self.model.has_regulatory_network():
            return {}  # No regulatory network

        result = {}
        for interaction in self.model.yield_interactions():
            target_value = interaction.solve(state)
            result[interaction.target.id] = target_value
        return result
```

**Changes:**
- ❌ Remove `self._extension` variable
- ❌ Remove legacy model checks
- ❌ Remove dual code paths
- ✅ Always assume RegulatoryExtension
- ✅ Always use string keys for objectives
- ✅ Simplified logic flow

### 2. RegulatoryExtension Simplifications

**Current:**
- Has `yield_reactions()` etc. for backwards compatibility
- Has complex checks for legacy models
- Some redundant properties

**Clean Version:**
```python
class RegulatoryExtension:
    """
    Wraps a Simulator and adds regulatory network capabilities.

    This is the ONLY way to integrate regulatory networks with metabolic models.
    No legacy GERM models supported.
    """

    def __init__(self, simulator: Simulator, regulatory_network: RegulatoryModel = None):
        """
        :param simulator: Simulator wrapping cobra/reframed model
        :param regulatory_network: Optional regulatory network
        """
        self._simulator = simulator
        self._regulators = {}
        self._targets = {}
        self._interactions = {}
        self._gpr_cache = {}

        if regulatory_network:
            self._load_regulatory_network(regulatory_network)

    # Metabolic delegation (all delegate to simulator)
    @property
    def reactions(self): return self._simulator.reactions

    @property
    def genes(self): return self._simulator.genes

    @property
    def metabolites(self): return self._simulator.metabolites

    @property
    def objective(self): return self._simulator.objective

    @property
    def simulator(self): return self._simulator

    def get_reaction(self, rxn_id): return self._simulator.get_reaction(rxn_id)
    def get_gene(self, gene_id): return self._simulator.get_gene(gene_id)
    def get_metabolite(self, met_id): return self._simulator.get_metabolite(met_id)

    # Regulatory network (stored internally)
    def add_regulator(self, regulator): ...
    def add_target(self, target): ...
    def add_interaction(self, interaction): ...
    def yield_interactions(self): ...
    def has_regulatory_network(self): return len(self._interactions) > 0

    # GPR caching (performance)
    def get_parsed_gpr(self, rxn_id):
        if rxn_id not in self._gpr_cache:
            gpr_str = self.get_gpr(rxn_id)
            self._gpr_cache[rxn_id] = parse_expression(gpr_str)
        return self._gpr_cache[rxn_id]
```

**Simplifications:**
- ❌ Remove `yield_reactions()`, `yield_metabolites()`, `yield_genes()` (not needed)
- ❌ Remove `is_metabolic()`, `is_regulatory()` (always True, always check has_regulatory_network())
- ✅ Pure delegation pattern
- ✅ Focused API

### 3. Factory Functions

**Current:**
- Multiple factory functions
- Some support legacy models

**Clean Version:**
```python
def from_cobra_model(cobra_model, regulatory_network=None):
    """
    Create RegulatoryExtension from COBRApy model.

    :param cobra_model: cobra.Model instance
    :param regulatory_network: Optional RegulatoryModel
    :return: RegulatoryExtension instance
    """
    simulator = get_simulator(cobra_model)
    return RegulatoryExtension(simulator, regulatory_network)

def from_reframed_model(reframed_model, regulatory_network=None):
    """
    Create RegulatoryExtension from reframed model.

    :param reframed_model: reframed.CBModel instance
    :param regulatory_network: Optional RegulatoryModel
    :return: RegulatoryExtension instance
    """
    simulator = get_simulator(reframed_model)
    return RegulatoryExtension(simulator, regulatory_network)

def load_integrated_model(metabolic_path, regulatory_path=None, backend='cobra'):
    """
    Load metabolic model and optionally add regulatory network.

    :param metabolic_path: Path to SBML file
    :param regulatory_path: Optional path to regulatory network file
    :param backend: 'cobra' or 'reframed'
    :return: RegulatoryExtension instance
    """
    if backend == 'cobra':
        import cobra
        cobra_model = cobra.io.read_sbml_model(metabolic_path)
        simulator = get_simulator(cobra_model)
    elif backend == 'reframed':
        from reframed.io.sbml import load_cbmodel
        ref_model = load_cbmodel(metabolic_path)
        simulator = get_simulator(ref_model)
    else:
        raise ValueError(f"Unknown backend: {backend}")

    regulatory_network = None
    if regulatory_path:
        regulatory_network = RegulatoryModel.from_file(regulatory_path)

    return RegulatoryExtension(simulator, regulatory_network)
```

**Changes:**
- ❌ Remove `from_cobra_model_with_regulation` (too wordy)
- ❌ Remove legacy factory functions
- ✅ Simple, clear names
- ✅ Only create RegulatoryExtension

### 4. Remove Files

**To Delete:**
- `src/mewpy/germ/models/metabolic.py` (legacy metabolic model)
- `src/mewpy/germ/models/simulator_model.py` (legacy simulator wrapper)
- Any backwards compatibility shims

**To Keep:**
- `src/mewpy/germ/models/regulatory_extension.py` (THE model)
- `src/mewpy/germ/models/regulatory.py` (for pure regulatory networks)
- `src/mewpy/germ/models/unified_factory.py` (clean factory functions)

---

## Example Usage (Clean)

### Example 1: Simple FBA with COBRApy
```python
import cobra
from mewpy.simulation import get_simulator
from mewpy.germ.models import RegulatoryExtension
from mewpy.germ.analysis import FBA

# Load metabolic model
cobra_model = cobra.io.load_model('textbook')
simulator = get_simulator(cobra_model)

# Create extension (no regulatory network)
model = RegulatoryExtension(simulator)

# Run FBA
fba = FBA(model)
solution = fba.optimize()
print(f"Growth rate: {solution.objective_value}")
```

### Example 2: RFBA with Regulatory Network
```python
import cobra
from mewpy.simulation import get_simulator
from mewpy.germ.models import RegulatoryExtension, RegulatoryModel
from mewpy.germ.analysis import RFBA

# Load metabolic model
cobra_model = cobra.io.read_sbml_model('ecoli_core.xml')
simulator = get_simulator(cobra_model)

# Load regulatory network
regulatory = RegulatoryModel.from_file('regulatory.csv')

# Create integrated model
model = RegulatoryExtension(simulator, regulatory)

# Run RFBA
rfba = RFBA(model)
solution = rfba.optimize()
print(f"Growth rate with regulation: {solution.objective_value}")
```

### Example 3: Using Factory Functions
```python
from mewpy.germ.models import load_integrated_model
from mewpy.germ.analysis import SRFBA

# Load everything in one call
model = load_integrated_model(
    metabolic_path='ecoli_core.xml',
    regulatory_path='regulatory.csv',
    backend='cobra'
)

# Run SRFBA
srfba = SRFBA(model)
solution = srfba.optimize()
```

---

## Implementation Plan

### Phase 1: Analysis Methods Cleanup (Priority: HIGH)

1. **RFBA** (`src/mewpy/germ/analysis/rfba.py`)
   - Remove `self._extension` variable
   - Remove legacy model checks in `__init__`
   - Simplify `build()` - always use RegulatoryExtension
   - Simplify `decode_constraints()` - only RegulatoryExtension path
   - Simplify `decode_regulatory_state()` - only regulatory network path
   - Change type hints to only accept RegulatoryExtension

2. **SRFBA** (`src/mewpy/germ/analysis/srfba.py`)
   - Same changes as RFBA
   - Remove dual code paths
   - Always fetch GPRs from simulator via extension

3. **PROM** (`src/mewpy/germ/analysis/prom.py`)
   - Same changes as RFBA
   - Simplify _max_rates, _optimize_ko

4. **CoRegFlux** (`src/mewpy/germ/analysis/coregflux.py`)
   - Same changes as RFBA
   - Simplify next_state

5. **FBA** (`src/mewpy/germ/analysis/fba.py`)
   - Already mostly clean
   - Remove legacy objective handling

### Phase 2: RegulatoryExtension Cleanup (Priority: MEDIUM)

1. Remove `yield_reactions()`, `yield_metabolites()`, `yield_genes()` if not needed
2. Remove `is_metabolic()`, `is_regulatory()` if not needed
3. Keep only essential methods

### Phase 3: Factory Functions (Priority: LOW)

1. Rename to cleaner names
2. Remove legacy support
3. Update documentation

### Phase 4: Delete Legacy Code (Priority: LOW)

1. Delete `metabolic.py`
2. Delete `simulator_model.py`
3. Update imports throughout

---

## Benefits of Clean Architecture

### 1. Simpler Code
- ❌ No dual code paths
- ❌ No isinstance() checks everywhere
- ❌ No `self._extension` variable
- ✅ Single, clear flow

### 2. Easier to Understand
- New developers see ONE way to do things
- No confusion about legacy vs new
- Clear separation of concerns

### 3. Easier to Maintain
- Less code to maintain
- No backwards compatibility burden
- Can evolve freely

### 4. Better Performance
- No redundant checks
- Simpler code paths
- More opportunities for optimization

### 5. Cleaner Documentation
- Document ONE way
- No "legacy vs new" sections
- Clear examples

---

## Migration Guide for Users

### Old Way (Legacy - NO LONGER SUPPORTED)
```python
from mewpy.io import read_model, Reader, Engines

# Legacy method
metabolic_reader = Reader(Engines.MetabolicSBML, 'model.xml')
regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, 'regulatory.csv', ...)
model = read_model(metabolic_reader, regulatory_reader)

# Use with analysis
rfba = RFBA(model)
```

### New Way (Clean Architecture)
```python
from mewpy.germ.models import load_integrated_model

# New method
model = load_integrated_model('model.xml', 'regulatory.csv', backend='cobra')

# Use with analysis
rfba = RFBA(model)
```

**Migration is simple:**
1. Replace `read_model()` with `load_integrated_model()`
2. That's it! Same API for analysis methods

---

## Timeline

- **Week 1:** Clean up analysis methods (RFBA, SRFBA, PROM, CoRegFlux)
- **Week 2:** Simplify RegulatoryExtension
- **Week 3:** Update factory functions and documentation
- **Week 4:** Delete legacy code and final testing

---

## Decision: Proceed?

**Recommendation:** YES - Clean up the code and remove backwards compatibility

**Rationale:**
- User explicitly doesn't care about backwards compatibility
- Cleaner code is easier to maintain and extend
- Performance and memory benefits
- Simpler for new users to understand

**Next Steps:**
1. Start with Phase 1 (analysis methods cleanup)
2. Test thoroughly after each cleanup
3. Update examples and documentation
4. Delete legacy code last

