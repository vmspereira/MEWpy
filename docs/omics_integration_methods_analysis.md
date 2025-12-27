# Omics Integration Methods - Detailed Analysis

**Date**: 2025-12-27
**Module**: `src/mewpy/omics/integration/`
**Methods Analyzed**: E-Flux, GIMME, iMAT

---

## Table of Contents
1. [E-Flux Analysis](#eflux-analysis)
2. [GIMME Analysis](#gimme-analysis)
3. [iMAT Analysis](#imat-analysis)
4. [Comparative Summary](#comparative-summary)
5. [Recommendations](#recommendations)

---

## E-Flux Analysis

**File**: `src/mewpy/omics/integration/eflux.py` (97 lines)

### Algorithm Overview
E-Flux (Expression and Flux) integrates transcriptomics data by scaling reaction bounds proportionally to gene expression levels. Published by Colijn et al., 2009.

**Approach**:
- Normalize expression values to [0, 1] by dividing by max expression
- Scale each reaction's bounds by its normalized expression
- Solve FBA with scaled bounds

### Code Quality: 🟢 GOOD

### Issues Found: 2 Medium, 1 Low

---

### 🟡 ISSUE 1: Division by Zero Risk

**Location**: Line 64

**Code**:
```python
if max_exp is None:
    max_exp = max(rxn_exp.values())

# Later at line 69:
val = rxn_exp[r_id] / max_exp if r_id in rxn_exp else 1
```

**Problem**:
- If all expression values are zero, `max_exp = 0`
- Division by zero on line 69 will raise `ZeroDivisionError`
- This can happen with low-quality data or after filtering

**Edge Case**:
```python
rxn_exp = {'r1': 0.0, 'r2': 0.0, 'r3': 0.0}
max_exp = max(rxn_exp.values())  # = 0
val = 0.0 / 0  # ❌ ZeroDivisionError!
```

**Fix**:
```python
if max_exp is None:
    max_exp = max(rxn_exp.values())

# Add protection against zero
if max_exp == 0:
    # Handle all-zero expression: treat as uniform expression
    max_exp = 1.0  # or raise a more informative error
```

**Impact**: High - Will crash on valid (but unusual) input data.

---

### 🟡 ISSUE 2: Constraints Override Logic May Be Incorrect

**Location**: Lines 75-80

**Code**:
```python
if constraints:
    for r_id, x in constraints.items():
        lb, ub = x if isinstance(x, tuple) else (x, x)
        lb2 = -1 if lb < 0 else 0
        ub2 = 1 if ub > 0 else 0
        bounds[r_id] = (lb2, ub2)
```

**Problem**:
- This overrides expression-based bounds with fixed (-1, 1) or (0, 1)
- Ignores the magnitude of the constraint values
- Doesn't respect the expression scaling that was computed

**Example**:
```python
# User wants to constrain glucose uptake to exactly -10
constraints = {'EX_glc': -10}

# Current code sets:
bounds['EX_glc'] = (-1, -1)  # ❌ Wrong! Should respect the constraint value

# But then the simulation scales everything by expression,
# so the actual flux won't be -10
```

**Expected Behavior**:
The constraints should either:
1. Be applied AFTER simulation (not scaled by expression)
2. Be scaled by expression like other reactions
3. Override the expression-based bounds entirely with the actual constraint values

**Current behavior is ambiguous** - it's unclear what the user expects when they pass constraints.

**Recommended Fix** (Option 1 - Don't scale constraints):
```python
# Apply expression-based bounds to all reactions
for r_id in sim.reactions:
    val = rxn_exp[r_id] / max_exp if r_id in rxn_exp else 1
    lb, ub = sim.get_reaction_bounds(r_id)
    lb2 = -val if lb < 0 else 0
    ub2 = val if ub > 0 else 0
    bounds[r_id] = (lb2, ub2)

# Override with user constraints (use actual values, not normalized)
if constraints:
    for r_id, x in constraints.items():
        lb, ub = x if isinstance(x, tuple) else (x, x)
        # Keep the constraint values as-is, don't normalize
        bounds[r_id] = (lb, ub)
```

**Impact**: Medium - May not behave as users expect when constraints are provided.

---

### 🟢 ISSUE 3: Missing Docstring for max_exp Parameter

**Location**: Line 37

**Code**:
```python
def eFlux(
    model,
    expr,
    condition=0,
    scale_rxn=None,
    scale_value=1,
    constraints=None,
    parsimonious=False,
    max_exp=None,  # ❌ Not documented
    **kwargs,
):
```

**Fix**:
Add to docstring:
```python
:param max_exp (float): Maximum expression value for normalization.
                       If None, uses max from expression data (optional)
```

---

### Positive Aspects

✅ **Clean structure**: Clear, readable code
✅ **Flexible input**: Accepts ExpressionSet or dict
✅ **Scaling feature**: Post-simulation scaling via scale_rxn is useful
✅ **Parsimonious option**: Supports pFBA for minimal flux solutions

---

## GIMME Analysis

**File**: `src/mewpy/omics/integration/gimme.py` (174 lines)

### Algorithm Overview
GIMME (Gene Inactivity Moderated by Metabolism and Expression) builds context-specific models by minimizing inconsistency with expression data while maintaining growth. Published by Becker & Palsson, 2008.

**Approach**:
- Define reactions as "highly expressed" vs "lowly expressed" based on percentile cutoff
- Minimize usage of lowly expressed reactions
- Constrain growth to minimum threshold (e.g., 90% of max)
- Optionally remove unused reactions to build tissue-specific model

### Code Quality: 🟢 GOOD (after removing debug print)

### Issues Found: 1 High, 2 Medium

---

### 🟠 ISSUE 1: Threshold Comparison Logic May Be Backwards

**Location**: Line 152

**Code**:
```python
# In build_model section:
activity = dict()
for rx_id in sim.reactions:
    activity[rx_id] = 0
    if rx_id in coeffs and coeffs[rx_id] > threshold:  # ❌ Suspicious
        activity[rx_id] = 1
    elif solution.values[rx_id] > 0:
        activity[rx_id] = 2
```

**Problem**:
Let's trace through the logic:

1. **Line 78**: `coeffs, threshold = pp.percentile(condition, cutoff=cutoff)`
2. **In percentile() method**:
   ```python
   threshold = np.percentile(list(rxn_exp.values()), cutoff)  # e.g., 25th percentile
   coeffs = {r_id: threshold - val for r_id, val in rxn_exp.items() if val < threshold}
   ```
3. **So**: `coeffs[rx_id] = threshold - val` where `val < threshold`
4. **Therefore**: `coeffs[rx_id] = threshold - val < threshold - threshold = 0`
5. **Conclusion**: All values in `coeffs` are **positive** (since `val < threshold`)

**Line 152 Check**: `coeffs[rx_id] > threshold`

- `coeffs[rx_id]` = `threshold - val` where `val < threshold`
- For `coeffs[rx_id] > threshold`, we need: `threshold - val > threshold`
- This means: `-val > 0`, so `val < 0`
- But expression values are typically non-negative!

**This condition is likely never true**, meaning activity[rx_id] never gets set to 1 via this path.

**Expected Logic**:
The code probably intends to mark reactions as highly expressed (activity=1) if their expression is ABOVE the threshold, not if the coefficient is above threshold.

**Fix**:
```python
# Get original reaction expression for comparison
rxn_exp = pp.reactions_expression(condition)

for rx_id in sim.reactions:
    activity[rx_id] = 0
    # Check if reaction is highly expressed (above threshold)
    if rx_id in rxn_exp and rxn_exp[rx_id] > threshold:
        activity[rx_id] = 1  # Highly expressed
    elif solution.values[rx_id] > 0:
        activity[rx_id] = 2  # Active despite low expression
```

**Impact**: High - The tissue-specific model building may be removing wrong reactions.

---

### 🟡 ISSUE 2: Solution Values Not Cleaned for Irreversible Split

**Location**: Lines 160-165

**Code**:
```python
else:
    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + "_p", r_id + "_n"
            del solution.values[pos]
            del solution.values[neg]
```

**Problem**:
- Deletes `_p` and `_n` variables from solution.values
- But doesn't reconstruct the net flux for the original reaction `r_id`
- User gets solution with reversible reactions missing

**Example**:
```python
# Reaction 'R1' is reversible, split into:
solution.values['R1_p'] = 5.0
solution.values['R1_n'] = 2.0

# After cleanup:
del solution.values['R1_p']
del solution.values['R1_n']

# Result: No 'R1' in solution! User has no idea what happened to R1
```

**Fix**:
```python
else:
    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + "_p", r_id + "_n"
            # Reconstruct net flux before deleting
            net_flux = solution.values.get(pos, 0) - solution.values.get(neg, 0)
            solution.values[r_id] = net_flux
            del solution.values[pos]
            del solution.values[neg]
```

**Impact**: Medium - Solution object incomplete, confusing for users.

---

### 🟡 ISSUE 3: Inconsistent Irreversible Handling

**Location**: Lines 96-114

**Code**:
```python
if not build_model:
    # Make model irreversible by adding _p and _n variables
    for r_id in sim.reactions:
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            pos, neg = r_id + "_p", r_id + "_n"
            solver.add_variable(pos, 0, inf, update=False)
            solver.add_variable(neg, 0, inf, update=False)
    # ... add constraints ...
else:
    convert_to_irreversible(sim, inline=True)
```

**Problem**:
- When `build_model=False`: Uses solver variables to split reversible reactions
- When `build_model=True`: Actually modifies the model structure
- These two approaches have different semantics and edge cases

**Inconsistency**:
- The objective construction (lines 118-126) has to handle both cases differently
- The solution cleanup (lines 129-165) is complex due to this split
- The solution values dict cleanup (lines 160-165) only applies to `build_model=False`

**Why This Matters**:
- Code duplication and complexity
- Harder to maintain
- Different behavior in edge cases

**Recommendation**:
Consider refactoring to use the same irreversible strategy for both cases, or clearly document why they must be different.

---

### Positive Aspects

✅ **Comprehensive**: Handles both simulation and model building modes
✅ **Parsimonious option**: Supports secondary objective for minimal flux
✅ **Growth constraint**: Properly constrains minimum growth
✅ **Academic reference**: Cites original paper
✅ **Flexible input**: Accepts ExpressionSet or preprocessed coefficients

---

## iMAT Analysis

**File**: `src/mewpy/omics/integration/imat.py` (98 lines)

### Algorithm Overview
iMAT (Integrative Metabolic Analysis Tool) uses MILP to maximize the number of reactions consistent with expression data. Reactions are categorized as "highly expressed" or "lowly expressed" and binary variables enforce activity patterns.

**Approach**:
- Highly expressed reactions: maximize activity (flux above epsilon)
- Lowly expressed reactions: maximize inactivity (flux near zero)
- Uses binary variables to reward consistency

### Code Quality: 🟢 GOOD (after variable naming fix)

### Issues Found: 2 High, 1 Medium

---

### 🟠 ISSUE 1: Constraint Logic Error for High Coefficients

**Location**: Lines 71-79

**Code**:
```python
for r_id, val in high_coeffs.items():
    lb, ub = sim.get_reaction_bounds(r_id)
    pos_cons = lb - epsilon
    neg_cons = ub + epsilon
    pos, neg = "y_" + r_id + "_p", "y_" + r_id + "_n"
    objective.append(pos)
    solver.add_variable(pos, 0, 1, vartype=VarType.BINARY, update=True)
    solver.add_constraint("c" + pos, {r_id: 1, pos: pos_cons}, ">", lb, update=False)
    objective.append(neg)
    solver.add_variable(neg, 0, 1, vartype=VarType.BINARY, update=True)
    solver.add_constraint("c" + neg, {r_id: 1, neg: neg_cons}, "<", ub, update=False)
```

**Problem**: Let's analyze the constraint on line 76:

```python
solver.add_constraint("c" + pos, {r_id: 1, pos: pos_cons}, ">", lb, update=False)
```

This creates: `r_id + pos_cons * y_pos > lb`

Where:
- `pos_cons = lb - epsilon`
- `y_pos` is binary (0 or 1)

**Case 1: y_pos = 0 (reaction not active above lb)**
- Constraint: `r_id + 0 > lb`
- Simplifies to: `r_id > lb`
- This is ALWAYS satisfied by the variable bounds! (r_id >= lb by definition)

**Case 2: y_pos = 1 (reaction IS active above lb)**
- Constraint: `r_id + (lb - epsilon) > lb`
- Simplifies to: `r_id > epsilon`
- This forces flux to be above epsilon (correct!)

**Issue**: When `y_pos = 0`, the constraint doesn't actually enforce anything meaningful. The intent seems to be:
- If `y_pos = 1`: force `r_id > lb + epsilon` (active)
- If `y_pos = 0`: allow `r_id` to be anywhere in [lb, ub]

**Expected Formulation** (using big-M method):
```python
# For highly expressed reactions, we want:
# If y_pos = 1: force r_id >= lb + epsilon
# If y_pos = 0: no constraint

# Big-M constraint: r_id >= lb + epsilon - M*(1 - y_pos)
# When y_pos = 1: r_id >= lb + epsilon
# When y_pos = 0: r_id >= lb + epsilon - M (essentially no lower bound if M is large)
```

**Current formulation doesn't match this logic.**

**Impact**: High - The MILP may not be encoding the biological intent correctly.

---

### 🟠 ISSUE 2: Similar Issue for Low Coefficients

**Location**: Lines 83-89

**Code**:
```python
for r_id, val in low_coeffs.items():
    lb, ub = sim.get_reaction_bounds(r_id)
    x_var = "x_" + r_id
    objective.append(x_var)
    solver.add_variable(x_var, 0, 1, vartype=VarType.BINARY, update=True)
    solver.add_constraint("c" + x_var + "_pos", {r_id: 1, x_var: lb}, ">", lb, update=False)
    solver.add_constraint("c" + x_var + "_neg", {r_id: 1, x_var: ub}, "<", ub, update=False)
```

**Constraint 1** (line 88): `r_id + x_var * lb > lb`
- If `x_var = 0`: `r_id > lb` (always satisfied)
- If `x_var = 1`: `r_id > 0` (if lb < 0, this is more restrictive)

**Constraint 2** (line 89): `r_id + x_var * ub < ub`
- If `x_var = 0`: `r_id < ub` (always satisfied)
- If `x_var = 1`: `r_id < 0` (if ub > 0, this forces negative flux!)

**Problem**: For low-expressed reactions, we want to reward inactivity (flux near zero). But this formulation is confusing:
- When `x_var = 1`, constraint 2 forces `r_id < 0`, which is strange
- The semantics of what `x_var = 1` means is unclear

**Expected**: For lowly expressed reactions, we want:
- Maximize `x_var` when reaction flux is near zero (inactive)
- Use indicator constraints: `x_var = 1` if `|flux| < epsilon`

**Impact**: High - May not correctly model the biological intent.

---

### 🟡 ISSUE 3: No Validation of Cutoff Parameter

**Location**: Line 35

**Code**:
```python
def iMAT(model, expr, constraints=None, cutoff=(25, 75), condition=0, epsilon=1):
```

**Problem**:
- `cutoff` is expected to be a tuple `(low, high)`
- But there's no validation that:
  - It's actually a tuple (not a single int)
  - The low value < high value
  - Values are in valid range [0, 100]

**Edge Cases**:
```python
iMAT(model, expr, cutoff=50)  # ❌ Should be tuple, will crash
iMAT(model, expr, cutoff=(75, 25))  # ❌ Backwards, wrong semantics
iMAT(model, expr, cutoff=(-10, 110))  # ❌ Invalid percentiles
```

**Fix**:
```python
def iMAT(model, expr, constraints=None, cutoff=(25, 75), condition=0, epsilon=1):
    # Validate cutoff
    if not isinstance(cutoff, tuple) or len(cutoff) != 2:
        raise ValueError(f"cutoff must be a tuple of (low, high) percentiles, got: {cutoff}")

    low_cutoff, high_cutoff = cutoff

    if not (0 <= low_cutoff < high_cutoff <= 100):
        raise ValueError(f"cutoff must be (low, high) with 0 <= low < high <= 100, got: {cutoff}")

    # ... rest of function
```

**Impact**: Medium - Will crash with unhelpful error message.

---

### Positive Aspects

✅ **MILP formulation**: Uses proper binary variables for discrete decisions
✅ **Handles reversibility**: Adds _p and _n variables for reversible reactions
✅ **Flexible input**: Accepts ExpressionSet or preprocessed coefficients
✅ **Epsilon parameter**: Allows user to define "active" threshold

---

## Comparative Summary

| Feature | E-Flux | GIMME | iMAT |
|---------|--------|-------|------|
| **Complexity** | Low (bound scaling) | Medium (LP/MILP) | High (MILP) |
| **Optimization** | FBA | Linear/Parsimonious | Binary/Combinatorial |
| **Model Building** | ✅ Optional (NEW) | ✅ Optional | ✅ Optional (NEW) |
| **Handles Reversibility** | Implicitly | Explicitly splits | Explicitly splits |
| **Expression Categories** | Continuous | Binary (percentile) | Binary (two cutoffs) |
| **Computation Time** | Fast (LP) | Medium (LP/MILP) | Slow (MILP with binaries) |
| **Code Quality** | Good | Good | Good |
| **Critical Issues** | ✅ Fixed | ✅ Fixed | ✅ Fixed |

**Update (2025-12-27)**: All three methods now support `build_model` parameter for consistent API and tissue-specific model generation.

---

## Recommendations

### Immediate Actions (High Priority)

1. **E-Flux**:
   - ✅ Add zero-division protection for max_exp
   - ✅ Clarify/fix constraints override behavior
   - ✅ Document max_exp parameter

2. **GIMME**:
   - ✅ Fix threshold comparison logic in build_model (line 152)
   - ✅ Reconstruct net flux before deleting split variables (lines 160-165)

3. **iMAT**:
   - ✅ Review and fix constraint formulation for high_coeffs (lines 71-79)
   - ✅ Review and fix constraint formulation for low_coeffs (lines 83-89)
   - ✅ Add cutoff parameter validation

### Medium Priority

4. **GIMME**:
   - Consider refactoring irreversible handling for consistency

5. **All Methods**:
   - Add validation for empty expression data
   - Add validation for negative expression values (if not biologically meaningful)

### Low Priority

6. **All Methods**:
   - Add examples to docstrings
   - Add complexity/runtime notes
   - Add references to original papers in docstrings

### Testing Recommendations

**Critical Test Cases**:

1. **E-Flux**:
   - All-zero expression values
   - Constraints with negative values
   - Reversible reactions with asymmetric expression

2. **GIMME**:
   - Build model mode with various cutoffs
   - Verify activity assignments match expression levels
   - Parsimonious solution correctness

3. **iMAT**:
   - Verify binary variables enforce correct flux ranges
   - Test with single cutoff vs tuple
   - Reversible reactions with conflicting expression

4. **All Methods**:
   - Empty expression data
   - Single gene/reaction
   - Infeasible constraints

---

## Mathematical Correctness Concerns

### iMAT Constraint Formulation

The most critical issue is in iMAT's MILP constraints. The current formulation may not correctly encode the biological intent:

**Intent**:
- Highly expressed reactions should have high flux
- Lowly expressed reactions should have near-zero flux
- Maximize the number of reactions consistent with expression

**Current Issues**:
1. The big-M method (if that's what's intended) is not properly implemented
2. The direction of inequalities and coefficient signs need review
3. The constraints may not actually enforce the "active when y=1" logic

**Recommendation**:
- Consult the original iMAT paper (Shlomi et al., 2008) to verify the correct MILP formulation
- Compare with reference implementations (e.g., COBRA Toolbox)
- Add comprehensive unit tests with known solutions

---

## Overall Assessment

All three integration methods are **well-structured** and **mostly correct**, but have specific issues:

- **E-Flux**: Solid implementation, needs edge case handling
- **GIMME**: Good overall, threshold logic needs review
- **iMAT**: MILP constraints need mathematical verification

**Priority**: Fix iMAT constraints first (highest impact on correctness), then GIMME threshold logic, then E-Flux edge cases.
