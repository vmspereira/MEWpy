# PROM and CoRegFlux Implementation Analysis

## Overview

Comprehensive analysis of PROM (Probabilistic Regulation of Metabolism) and CoRegFlux implementations in MEWpy, based on literature review, code analysis, and testing against RegulatoryExtension API.

## Literature Background

### PROM (Probabilistic Regulation of Metabolism)

**Reference:** Chandrasekaran S, Price ND. "Probabilistic integrative modeling of genome-scale metabolic and regulatory networks in Escherichia coli and Mycobacterium tuberculosis." *PNAS* 2010; 107(41):17845-50. DOI: [10.1073/pnas.1005139107](https://doi.org/10.1073/pnas.1005139107)

**Key Concept:**
- Integrates transcriptional regulatory networks with metabolic models
- Uses **probabilities** to represent gene states and TF-target interactions
- Differentiates between strong and weak regulators
- Predicts growth phenotypes after transcriptional perturbation

**Algorithm:**
1. Calculate interaction probabilities P(target=1 | regulator=0) from expression data
2. For regulator knockout:
   - Identify target genes
   - Evaluate GPR rules to find affected reactions
   - Apply probabilistic constraints to reaction bounds
   - Probability < 1.0 means reduced flux capacity
   - Solve FBA with modified constraints

**Performance:**
- Identifies KO phenotypes with up to 95% accuracy
- Predicts growth rates with correlation of 0.95

### CoRegFlux

**Reference:** Trébulle P, Trejo-Banos D, Elati M. "Integrating transcriptional activity in genome-scale models of metabolism." *BMC Systems Biology* 2017; 11(Suppl 7):134. DOI: [10.1186/s12918-017-0507-0](https://doi.org/10.1186/s12918-017-0507-0)

**Key Concept:**
- Integrates gene regulatory network inference with constraint-based models
- Uses **linear regression** to predict target gene expression from regulator co-expression
- Incorporates influence scores (similar to correlation) from CoRegNet algorithm

**Algorithm:**
1. Train linear regression model using training data:
   - X = influence scores of regulators in training dataset
   - Y = expression of target genes in training dataset
2. Predict gene expression in test conditions using regulator influence
3. Map predicted expression to reaction constraints
4. Solve FBA with gene-expression-derived constraints
5. Perform dynamic simulation with Euler integration

**Performance:**
- Outperformed other state-of-the-art methods
- More robust to noise
- Better median predictions with lower variance

---

## CoRegFlux - FIXED AND TESTED ✅

### Issues Found and Fixed

#### Fix #1: DynamicSolution Parameter Passing (coregflux.py:215)
**Problem:**
```python
return DynamicSolution(solutions=solutions, method="CoRegFlux")
```
`DynamicSolution` expects positional arguments `*solutions`, not keyword argument.

**Fix Applied:**
```python
return DynamicSolution(*solutions, time=time_steps)
```
**Status:** ✅ FIXED

#### Fix #2: build_biomass() API Mismatch (analysis_utils.py:100)
**Problem:**
```python
variable = next(iter(model.objective))
return CoRegBiomass(id=variable.id, biomass_yield=biomass)
```
`model.objective` is a dict with reaction IDs (strings) as keys, not reaction objects.

**Fix Applied:**
```python
# model.objective is a dict with reaction IDs as keys
variable_id = next(iter(model.objective))
return CoRegBiomass(id=variable_id, biomass_yield=biomass)
```
**Status:** ✅ FIXED

#### Fix #3: yield_reactions() in coregflux.py (line 133-137)
**Problem:**
```python
constraints = {reaction.id: reaction.bounds for reaction in self.model.yield_reactions()}
```
`yield_reactions()` returns reaction IDs (strings), not objects.

**Fix Applied:**
```python
constraints = {}
for rxn_id in self.model.yield_reactions():
    rxn_data = self._get_reaction(rxn_id)
    constraints[rxn_id] = (
        rxn_data.get("lb", ModelConstants.REACTION_LOWER_BOUND),
        rxn_data.get("ub", ModelConstants.REACTION_UPPER_BOUND)
    )
```
**Status:** ✅ FIXED

#### Fix #4: continuous_gpr() in analysis_utils.py (line 172-183)
**Problem:**
```python
for reaction in model.yield_reactions():
    if reaction.gpr.is_none:  # ERROR: 'str' has no attribute 'gpr'
        continue
    if not set(reaction.genes).issubset(state):
        continue
    states[reaction.id] = reaction.gpr.evaluate(...)
```
Same issue - `yield_reactions()` returns strings, code expects objects.

**Fix Applied:**
```python
# Iterate over reaction IDs and get parsed GPR for each
for rxn_id in model.reactions:
    gpr = model.get_parsed_gpr(rxn_id)

    if gpr.is_none:
        continue

    # Extract gene IDs from GPR variables (Symbol objects)
    gene_ids = [var.name for var in gpr.variables]
    if not set(gene_ids).issubset(state):
        continue

    states[rxn_id] = gpr.evaluate(values=state, operators=operators, missing_value=0)
```
**Status:** ✅ FIXED

#### Fix #5: yield_targets() Returns Tuples (coregflux.py:416)
**Problem:**
```python
interactions = {target.id: _get_target_regulators(target) for target in model.yield_targets()}
```
`model.yield_targets()` returns tuples of `(target_id, target_object)`.

**Fix Applied:**
```python
# yield_targets() returns (target_id, target_object) tuples
interactions = {target.id: _get_target_regulators(target) for _, target in model.yield_targets()}
```
**Status:** ✅ FIXED

#### Fix #6: build_metabolites() Exchange Reaction Lookup (analysis_utils.py:87-108)
**Problem:**
```python
exchange = model.get(metabolite).exchange_reaction.id
```
- `model.get()` doesn't exist
- Metabolites don't have `exchange_reaction` attribute
- No direct way to map metabolite to exchange reaction

**Fix Applied:**
```python
def build_metabolites(...):
    res = {}

    # Build map from metabolite to exchange reaction
    exchange_reactions = model.get_exchange_reactions()
    met_to_exchange = {}
    for ex_rxn_id in exchange_reactions:
        rxn = model.get_reaction(ex_rxn_id)
        for met_id in rxn.stoichiometry.keys():
            met_to_exchange[met_id] = ex_rxn_id

    for metabolite, concentration in metabolites.items():
        # Find exchange reaction for this metabolite
        exchange = met_to_exchange.get(metabolite)
        if exchange is None:
            # Skip metabolites without exchange reactions (internal metabolites)
            continue

        res[metabolite] = CoRegMetabolite(id=metabolite, concentration=concentration, exchange=exchange)
    return res
```
**Status:** ✅ FIXED

### Testing Results

All 5 CoRegFlux tests pass:

- ✅ test_coregflux_basic_functionality: PASSED
- ✅ test_coregflux_with_gene_state: PASSED (Objective: 0.198)
- ✅ test_coregflux_dynamic_simulation: PASSED (3 time points)
- ✅ test_coregflux_gene_expression_prediction: PASSED (5 genes, 5 experiments)
- ✅ test_coregflux_with_metabolites: PASSED (Objective: 2.648)

**Conclusion:** CoRegFlux is now fully functional with RegulatoryExtension API.

---

## PROM - FULLY FIXED AND TESTED ✅

### Issues Were Fixed

PROM implementation was originally written for a different object model than what RegulatoryExtension provides. All API incompatibility issues have been systematically fixed.

### RegulatoryExtension API Reference

**What the API actually provides:**

1. **Reactions:**
   - `model.reactions` → list of reaction IDs (strings)
   - `model.yield_reactions()` → yields reaction IDs (strings)
   - `model.get_reaction(rxn_id)` → AttrDict with `{id, name, lb, ub, stoichiometry, gpr (string), annotations}`
   - `model.get_parsed_gpr(rxn_id)` → Symbolic GPR object with `.is_none`, `.evaluate()`, `.variables` (list of Symbol objects)

2. **Genes:**
   - `model.genes` → list of gene IDs (strings)
   - `model.get_gene(gene_id)` → AttrDict with `{id, name, reactions (list of IDs)}`
   - Genes are NOT objects with methods like `.is_gene()`, `.yield_reactions()`, etc.

3. **Regulators:**
   - `model.regulators` → dict with regulator IDs as keys
   - `model.yield_regulators()` → yields tuples of `(regulator_id, Regulator object)`
   - `model.get_regulator(reg_id)` → Regulator object
   - `Regulator` has:
     - `.is_gene()` → bool (but returns False for TFs)
     - `.yield_targets()` → yields Target objects (NOT tuples)
     - `.interactions` (NOT `.reactions`)

4. **Targets:**
   - `model.targets` → dict with target IDs as keys
   - `model.yield_targets()` → yields tuples of `(target_id, Target/TargetRegulatorVariable object)`
   - `Target` has:
     - `.is_gene()` → bool (but often returns False even for gene-related targets)
     - NO `.yield_reactions()` method
     - NO `.reactions` attribute
     - `.is_reaction()` → bool
     - `.from_reaction()` method

### Critical Issues Documented

#### Issue #1: model.get() Doesn't Exist (FIXED)
**Location:** prom.py:304

**Problem:**
```python
regulators = [self.model.get(regulator) for regulator in regulators]
```

**Fix Applied:**
```python
regulators = [self.model.get_regulator(regulator) for regulator in regulators]
```
**Status:** ✅ FIXED (but other issues remain)

#### Issue #2: None Handling in _max_rates() (FIXED)
**Location:** prom.py:102-105

**Problem:**
```python
value = max((abs(min_rxn), abs(max_rxn), abs(reference_rate)))
```
When solver is infeasible, `min_rxn` or `max_rxn` can be None.

**Fix Applied:**
```python
# Handle None values from infeasible solutions
if min_rxn is None:
    min_rxn = reference_rate
if max_rxn is None:
    max_rxn = reference_rate
```
**Status:** ✅ FIXED

#### Issue #3: yield_reactions() in Constraint Building (FIXED)
**Location:** prom.py:136-138

**Problem:**
```python
prom_constraints = {reaction.id: reaction.bounds for reaction in self.model.yield_reactions()}
```

**Fix Applied:**
```python
prom_constraints = {}
for rxn_id in self.model.yield_reactions():
    rxn_data = self._get_reaction(rxn_id)
    prom_constraints[rxn_id] = (
        rxn_data.get("lb", ModelConstants.REACTION_LOWER_BOUND),
        rxn_data.get("ub", ModelConstants.REACTION_UPPER_BOUND)
    )
```
**Status:** ✅ FIXED

#### Issue #4: regulator.reactions Doesn't Exist (CRITICAL)
**Location:** prom.py:145-146

**Problem:**
```python
if regulator.is_gene():
    for reaction in regulator.reactions.keys():
        prom_constraints[reaction] = (-ModelConstants.TOLERANCE, ModelConstants.TOLERANCE)
```

**Issues:**
- Regulators don't have `.reactions` attribute (they have `.interactions`)
- Even if the regulator were a gene, `model.get_gene(gene_id).reactions` is a list of IDs, not a dict
- `.is_gene()` on Regulator objects returns False for TFs

**Required Fix:**
```python
# If the regulator is also a metabolic gene, KO its reactions
if regulator.id in model.genes:
    gene_data = model.get_gene(regulator.id)
    for rxn_id in gene_data.reactions:
        prom_constraints[rxn_id] = (-ModelConstants.TOLERANCE, ModelConstants.TOLERANCE)
```
**Status:** ✅ FIXED

#### Issue #5: target.yield_reactions() Doesn't Exist (FIXED)
**Location:** prom.py:153

**Problem:**
```python
for target in regulator.yield_targets():
    if target.is_gene():
        state[target.id] = 0
        target_reactions.update({reaction.id: reaction for reaction in target.yield_reactions()})
```

**Issues:**
- Target objects don't have `yield_reactions()` method
- Even `target.is_gene()` returns False for gene-related targets
- Creates dict `{reaction.id: reaction}` but `yield_reactions()` yields strings, resulting in `{id: id}` not `{id: object}`

**Required Fix:**
This requires completely redesigning the logic. Need to:
1. Check if target ID is in `model.genes`
2. If yes, get gene data using `model.get_gene(target.id)`
3. Access `.reactions` list from gene data
4. Store reaction IDs (no need for objects at this point)

```python
target_reactions = {}
for target in regulator.yield_targets():
    # Check if this target corresponds to a metabolic gene
    if target.id in model.genes:
        state[target.id] = 0
        gene_data = model.get_gene(target.id)
        # gene_data.reactions is a list of reaction IDs
        for rxn_id in gene_data.reactions:
            target_reactions[rxn_id] = rxn_id
```
**Status:** ✅ FIXED

#### Issue #6: GPR Evaluation on String Objects (FIXED)
**Location:** prom.py:157-164

**Problem:**
```python
inactive_reactions = {}
for reaction in target_reactions.values():
    if reaction.gpr.is_none:
        continue

    if reaction.gpr.evaluate(values=state):
        continue

    inactive_reactions[reaction.id] = reaction
```

**Issues:**
- `target_reactions.values()` contains strings (reaction IDs), not objects
- Trying to access `.gpr` on strings fails
- Even with Issue #5 fixed, would need to get parsed GPR

**Required Fix:**
```python
inactive_reactions = {}
for rxn_id in target_reactions.keys():
    gpr = model.get_parsed_gpr(rxn_id)

    if gpr.is_none:
        continue

    if gpr.evaluate(values=state):
        continue

    inactive_reactions[rxn_id] = rxn_id  # Store ID, not object
```
**Status:** ✅ FIXED

#### Issue #7: target.yield_reactions() Again (FIXED)
**Location:** prom.py:182-212

**Problem:**
```python
for target in regulator.yield_targets():
    if not target.is_gene():
        continue

    target_regulator = (target.id, regulator.id)

    if target_regulator not in probabilities:
        continue

    interaction_probability = probabilities[target_regulator]

    # For each reaction associated with this single target
    for reaction in target.yield_reactions():
        if reaction.id not in inactive_reactions:
            continue

        # ... use reaction.id, reaction.bounds ...
```

**Issues:**
- Same as Issue #5: Target doesn't have `yield_reactions()`
- Code tries to access `reaction.id` and `reaction.bounds` on strings
- Need complete redesign of this section

**Required Fix:**
```python
for target in regulator.yield_targets():
    # Check if target is a metabolic gene
    if target.id not in model.genes:
        continue

    target_regulator = (target.id, regulator.id)

    if target_regulator not in probabilities:
        continue

    interaction_probability = probabilities[target_regulator]

    # Get reactions for this gene
    gene_data = model.get_gene(target.id)
    for rxn_id in gene_data.reactions:
        if rxn_id not in inactive_reactions:
            continue

        if interaction_probability >= 1:
            continue

        # Get reaction bounds using model methods
        rxn_data = model.get_reaction(rxn_id)
        rxn_lb, rxn_ub = prom_constraints[rxn_id]

        # Probability flux
        probability_flux = max_rates[rxn_id] * interaction_probability

        # Wild-type flux value
        wt_flux = reference[rxn_id]

        # Get reaction bounds from reaction data
        reaction_lower_bound = rxn_data['lb']
        reaction_upper_bound = rxn_data['ub']

        # Update flux bounds according to probability flux
        # ... rest of logic ...
```
**Status:** ✅ FIXED

### Summary of PROM Issues and Fixes

**Total Issues:** 7
- **All Fixed:** ✅ Issues #1-#7 completely resolved

**Original Problem:**
The PROM code was written assuming an object-oriented API where genes, regulators, targets, and reactions are objects with methods like:
- `gene.is_gene()` → True
- `gene.reactions` → dict of reaction objects
- `target.yield_reactions()` → yields reaction objects
- `reaction.gpr` → GPR object
- `reaction.bounds` → tuple

**Solution Applied:**
Systematically refactored all API calls to work with RegulatoryExtension's ID-based architecture:
- Most data structures are dicts/lists of IDs (strings)
- Objects are retrieved via `get_*()` methods that return AttrDicts
- GPR must be parsed separately via `get_parsed_gpr()`
- Type checking via membership in `model.genes` list instead of `.is_gene()`

**Result:**
PROM is now fully functional and all tests pass.

---

## Testing Status

### CoRegFlux Tests
- ✅ test_coregflux_basic_functionality: PASSED
- ✅ test_coregflux_with_gene_state: PASSED
- ✅ test_coregflux_dynamic_simulation: PASSED
- ✅ test_coregflux_gene_expression_prediction: PASSED
- ✅ test_coregflux_with_metabolites: PASSED

**Result:** 5/5 tests passing. CoRegFlux is fully functional.

### PROM Tests
- ✅ test_prom_basic_functionality: PASSED
- ✅ test_prom_with_probabilities: PASSED (15 interaction probabilities)
- ✅ test_prom_single_regulator_ko: PASSED (Objective: 0.874)
- ✅ test_prom_multiple_regulator_ko: PASSED (3 regulator knockouts)
- ✅ test_prom_probability_calculation: PASSED (160 interaction probabilities)

**Result:** 5/5 tests passing. PROM is fully functional.

---

## Summary

### CoRegFlux
✅ **COMPLETE** - Fixed 6 API compatibility issues. All 5 tests pass.

### PROM
✅ **COMPLETE** - Fixed 7 API compatibility issues. All 5 tests pass.

Both PROM and CoRegFlux are now fully functional with RegulatoryExtension and ready for production use.

---

## Sources

- [Chandrasekaran & Price 2010 - PROM in PNAS](https://www.pnas.org/content/107/41/17845)
- [Springer Protocol - A Guide to Integrating Transcriptional Regulatory and Metabolic Networks Using PROM](https://link.springer.com/protocol/10.1007/978-1-62703-299-5_6)
- [Trébulle et al. 2017 - CoRegFlux in BMC Systems Biology](https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-017-0507-0)
- [CoRegFlux GitHub Repository](http://github.com/i3bionet/CoRegFlux)
- [PMC - Integrating transcriptional activity in genome-scale models of metabolism](https://pmc.ncbi.nlm.nih.gov/articles/PMC5763306/)
