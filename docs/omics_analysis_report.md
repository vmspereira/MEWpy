# Omics Module Code Analysis Report

**Date**: 2025-12-27
**Module**: `src/mewpy/omics/`
**Total Lines**: ~852 lines
**Status**: ✅ PASSES ALL LINTING

---

## Executive Summary

The omics module implements gene expression data handling and integration algorithms (E-Flux, GIMME, iMAT) for constraint-based metabolic modeling. The code is **well-structured and clean**, passing all linting checks (flake8, black, isort) with zero issues.

**Overall Quality**: 🟢 **GOOD**

**Priority Breakdown**:
- 🔴 **CRITICAL**: 0 issues
- 🟠 **HIGH**: 3 issues (debug code, potential bugs, missing validation)
- 🟡 **MEDIUM**: 6 issues (code quality, optimization)
- 🟢 **LOW**: 4 issues (documentation, type hints)

---

## Module Structure

```
src/mewpy/omics/
├── __init__.py                (4 lines) - Clean exports
├── expression.py              (479 lines) - Main expression data handling
└── integration/
    ├── __init__.py            (0 lines) - Empty
    ├── eflux.py               (97 lines) - E-Flux algorithm
    ├── gimme.py               (174 lines) - GIMME algorithm
    └── imat.py                (98 lines) - iMAT algorithm
```

---

## 🟠 HIGH PRIORITY ISSUES

### 1. **Debug Print Statement Left in Production Code**

**Location**: `src/mewpy/omics/integration/gimme.py`, line 82

**Issue**:
```python
def GIMME(...):
    # ... setup code ...
    print(coeffs)  # ❌ Debug print statement!
    solver = solver_instance(sim)
```

**Problem**:
- Print statement in production code
- Will clutter output in automated pipelines
- Not controlled by logging framework

**Fix**:
```python
# Remove the print statement or replace with proper logging:
import logging
logger = logging.getLogger(__name__)

def GIMME(...):
    # ... setup code ...
    logger.debug(f"Coefficients: {coeffs}")  # ✅ Use logging
    solver = solver_instance(sim)
```

**Impact**: Minor annoyance, but unprofessional in production code.

---

### 2. **Potential Bug: p_values Property Check**

**Location**: `src/mewpy/omics/expression.py`, line 171

**Issue**:
```python
@property
def p_values(self):
    """Returns the numpy array of p-values."""
    if not self._p_values.all():  # ❌ WRONG CHECK!
        raise ValueError("No p-values defined.")
    else:
        return self._p_values
```

**Problem**:
- `.all()` checks if all values are truthy (non-zero)
- But `_p_values` could be `None` → will raise `AttributeError`
- Should check if `_p_values is None` instead

**Correct Check**:
```python
@property
def p_values(self):
    """Returns the numpy array of p-values."""
    if self._p_values is None:  # ✅ Correct check
        raise ValueError("No p-values defined.")
    else:
        return self._p_values
```

**Impact**: Will crash with `AttributeError: 'NoneType' object has no attribute 'all'` when p-values are not defined.

---

### 3. **Missing Input Validation in ExpressionSet Constructor**

**Location**: `src/mewpy/omics/expression.py`, lines 39-62

**Issue**:
```python
def __init__(self, identifiers: list, conditions: list, expression: np.array, p_values: np.array = None):
    # Checks expression shape matches identifiers/conditions
    n = len(identifiers)
    m = len(conditions)
    if expression.shape != (n, m):
        raise ValueError(...)

    # But doesn't check:
    # - Empty lists
    # - Duplicate identifiers
    # - Duplicate conditions
    # - p_values shape if provided
```

**Recommended Validation**:
```python
def __init__(self, identifiers: list, conditions: list, expression: np.array, p_values: np.array = None):
    # Validate non-empty
    if not identifiers or not conditions:
        raise ValueError("Identifiers and conditions cannot be empty")

    # Check for duplicates
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Duplicate identifiers found")

    if len(conditions) != len(set(conditions)):
        raise ValueError("Duplicate conditions found")

    # Check shape
    n, m = len(identifiers), len(conditions)
    if expression.shape != (n, m):
        raise ValueError(f"Expression shape {expression.shape} doesn't match ({n},{m})")

    # Validate p_values shape if provided
    if p_values is not None:
        expected_p_cols = len(list(combinations(conditions, 2)))
        if p_values.shape != (n, expected_p_cols):
            raise ValueError(f"p_values shape {p_values.shape} doesn't match expected ({n},{expected_p_cols})")

    # ... rest of initialization ...
```

**Impact**: Could lead to silent errors or confusing behavior downstream.

---

## 🟡 MEDIUM PRIORITY ISSUES

### 4. **Inefficient Repeated Computation in Preprocessing.percentile()**

**Location**: `src/mewpy/omics/expression.py`, lines 374-389

**Issue**:
```python
def percentile(self, condition=None, cutoff=25):
    if type(cutoff) is tuple:
        coef = []
        thre = []
        for cut in cutoff:
            rxn_exp = self.reactions_expression(condition)  # ❌ Computed repeatedly!
            threshold = np.percentile(list(rxn_exp.values()), cut)
            coeffs = {r_id: threshold - val for r_id, val in rxn_exp.items() if val < threshold}
            coef.append(coeffs)
            thre.append(threshold)
```

**Problem**:
- `reactions_expression(condition)` is called once per cutoff value
- This recalculates gene→reaction expression mapping every time
- Very wasteful for tuple cutoffs like `(25, 75)`

**Fix**:
```python
def percentile(self, condition=None, cutoff=25):
    # Compute reaction expression ONCE
    rxn_exp = self.reactions_expression(condition)

    if type(cutoff) is tuple:
        coef = []
        thre = []
        for cut in cutoff:
            threshold = np.percentile(list(rxn_exp.values()), cut)  # ✅ Reuse rxn_exp
            coeffs = {r_id: threshold - val for r_id, val in rxn_exp.items() if val < threshold}
            coef.append(coeffs)
            thre.append(threshold)
        coeffs = tuple(coef)
        threshold = tuple(thre)
    else:
        threshold = np.percentile(list(rxn_exp.values()), cutoff)
        coeffs = {r_id: threshold - val for r_id, val in rxn_exp.items() if val < threshold}
    return coeffs, threshold
```

---

### 5. **Use of `type()` Instead of `isinstance()`**

**Location**: `src/mewpy/omics/expression.py`, line 374

**Issue**:
```python
if type(cutoff) is tuple:  # ❌ Should use isinstance()
```

**Problem**:
- `type()` doesn't respect subclasses
- `isinstance()` is the Pythonic way

**Fix**:
```python
if isinstance(cutoff, tuple):  # ✅ Correct
```

---

### 6. **Mutation of Input Array in quantile_binarization()**

**Location**: `src/mewpy/omics/expression.py`, lines 466-479

**Issue**:
```python
def quantile_binarization(expression: np.ndarray, q: float = 0.33) -> np.ndarray:
    threshold = np.quantile(expression, q)

    threshold_mask = expression >= threshold
    expression[threshold_mask] = 1      # ❌ Mutates input!
    expression[~threshold_mask] = 0     # ❌ Mutates input!
    return expression
```

**Problem**:
- Function mutates the input array
- Caller's data is unexpectedly modified
- Side effects are undocumented

**Fix**:
```python
def quantile_binarization(expression: np.ndarray, q: float = 0.33) -> np.ndarray:
    """
    Binarizes the expression matrix using the q-th quantile threshold.

    :param expression: Expression matrix (will NOT be modified)
    :param q: Quantile to compute
    :return: NEW binarized expression matrix
    """
    threshold = np.quantile(expression, q)

    # Create a copy to avoid mutating input
    binary_expression = expression.copy()

    threshold_mask = binary_expression >= threshold
    binary_expression[threshold_mask] = 1
    binary_expression[~threshold_mask] = 0
    return binary_expression
```

---

### 7. **Redundant None Checks with Default Parameters**

**Location**: `src/mewpy/omics/expression.py`, lines 426-436

**Issue**:
```python
def knn_imputation(
    expression: np.ndarray,
    missing_values: float = None,  # Default is None
    n_neighbors: int = 5,
    weights: str = "uniform",
    metric: str = "nan_euclidean",
):
    # ... import ...

    if missing_values is None:  # ❌ Redundant - just use np.nan as default
        missing_values = np.nan

    if n_neighbors is None:     # ❌ Will never be None (default is 5)
        n_neighbors = 5

    if weights is None:         # ❌ Will never be None (default is "uniform")
        weights = "uniform"

    if metric is None:          # ❌ Will never be None
        metric = "nan_euclidean"
```

**Fix**:
```python
def knn_imputation(
    expression: np.ndarray,
    missing_values: float = np.nan,  # ✅ Use np.nan directly
    n_neighbors: int = 5,
    weights: str = "uniform",
    metric: str = "nan_euclidean",
):
    # Remove all the redundant None checks
    imputation = KNNImputer(
        missing_values=missing_values,
        n_neighbors=n_neighbors,
        weights=weights,
        metric=metric
    )
    return imputation.fit_transform(expression)
```

---

### 8. **Inconsistent Return Type in get_condition()**

**Location**: `src/mewpy/omics/expression.py`, lines 76-104

**Issue**:
```python
def get_condition(self, condition: Union[int, str] = None, **kwargs):
    # ... gets values ...

    form = kwargs.get("format", "dict")
    if form and condition is not None:  # ❌ Complex conditional logic
        if form == "list":
            return values.tolist()
        elif form == "dict":
            return dict(zip(self._identifiers, values.tolist()))
        else:
            return values  # numpy array
    else:
        return values  # numpy array
```

**Problem**:
- Returns different types (list, dict, np.array) based on `format` parameter
- Logic for when to apply format is confusing
- Default format is "dict" but sometimes returns array

**Recommendation**: Simplify logic and document return types clearly:
```python
def get_condition(self, condition: Union[int, str] = None, format: str = None):
    """
    Retrieves omics data for a specific condition.

    :param condition: Condition identifier (int index or str name).
                     If None, returns all data.
    :param format: Output format: "dict", "list", or None for numpy array
    :return: Expression values in requested format
    """
    if isinstance(condition, int):
        values = self[:, condition]
    elif isinstance(condition, str):
        values = self[:, self._condition_index[condition]]
    else:
        values = self[:, :]

    # Always apply format conversion if specified and single condition
    if format and condition is not None:
        if format == "list":
            return values.tolist()
        elif format == "dict":
            return dict(zip(self._identifiers, values.tolist()))

    return values  # Default: numpy array
```

---

### 9. **Missing Error Handling for Missing Conditions/Identifiers**

**Location**: `src/mewpy/omics/expression.py`, line 90

**Issue**:
```python
elif isinstance(condition, str):
    values = self[:, self._condition_index[condition]]  # ❌ KeyError if not found
```

**Fix**:
```python
elif isinstance(condition, str):
    if condition not in self._condition_index:
        raise ValueError(f"Unknown condition: {condition}. Available: {self._conditions}")
    values = self[:, self._condition_index[condition]]
```

---

## 🟢 LOW PRIORITY ISSUES

### 10. **Incomplete Docstrings**

Many docstrings use `[description]` placeholders:

**Examples**:
```python
# Line 167
def p_values(self):
    """Returns the numpy array of p-values.

    Raises:
        ValueError: [description]  # ❌ Not filled in
    """
```

```python
# Line 181
def p_values(self, p_values: np.array):
    """Sets p-values

    Args:
        p_values (np.array): [description]  # ❌ Not filled in

    Raises:
        ValueError: [description]  # ❌ Not filled in
    """
```

**Fix**: Replace `[description]` with actual descriptions.

---

### 11. **Missing Type Hints**

Several functions lack complete type hints:

```python
# Line 239
def apply(self, function: None):  # ❌ Should be Callable, not None
    """Apply a function to all expression values."""
```

**Fix**:
```python
from typing import Callable, Optional

def apply(self, function: Optional[Callable[[float], float]] = None):
    """Apply a function to all expression values."""
```

---

### 12. **Variable Naming: `object` vs `objective`**

**Location**: `src/mewpy/omics/integration/imat.py`, line 93

**Issue**:
```python
object = {x: 1 for x in objective}  # ❌ Typo: 'object' should be 'objective'
```

**Problem**:
- Shadows Python builtin `object`
- Inconsistent with variable name used elsewhere in code

**Fix**:
```python
objective_dict = {x: 1 for x in objective}
solution = solver.solve(objective_dict, minimize=False, constraints=constraints)
```

---

### 13. **Typo in Docstring**

**Location**: `src/mewpy/omics/expression.py`, line 333

**Issue**:
```python
# Line 333
For Order 2, thresholding ofgene expression is followed by its
#                          ^^^^^^^ Missing space: "of gene"
```

---

## Algorithm-Specific Issues

### GIMME Algorithm

**Line 152**: Threshold comparison may be incorrect:
```python
if rx_id in coeffs and coeffs[rx_id] > threshold:
```

This checks if the coefficient (which is `threshold - val`) is greater than the threshold itself. This seems backwards - coefficients are negative differences for low-expressed reactions. Review the logic.

### iMAT Algorithm

**Lines 71-79**: Variable naming is confusing:
```python
pos_cons = lb - epsilon
neg_cons = ub + epsilon
pos, neg = "y_" + r_id + "_p", "y_" + r_id + "_n"
```

The variables `pos_cons` and `neg_cons` are used in constraints but their meaning is unclear. Consider renaming to `lower_bound_offset` and `upper_bound_offset`.

---

## Positive Aspects

### ✅ Code Quality Strengths:

1. **Clean Code Structure**:
   - Well-organized module hierarchy
   - Clear separation of concerns (expression handling vs. integration algorithms)
   - Good use of classes and functions

2. **Linting**:
   - ✅ Passes flake8 (0 issues)
   - ✅ Passes black formatting
   - ✅ Passes isort

3. **Documentation**:
   - Good docstrings for main functions
   - Academic references included (GIMME algorithm)
   - Clear parameter descriptions

4. **Error Handling**:
   - Shape validation in ExpressionSet constructor
   - Import error handling with helpful messages (sklearn imports)

5. **Type Hints**:
   - Present in many function signatures
   - Helps with IDE support

6. **Pandas Integration**:
   - Good DataFrame interoperability
   - CSV file loading support

---

## Recommendations Summary

### Immediate Actions (High Priority):
1. ✅ Remove debug print statement (gimme.py:82)
2. ✅ Fix p_values property check (expression.py:171)
3. ✅ Add input validation to ExpressionSet.__init__

### Short Term (Medium Priority):
4. Optimize percentile() to avoid repeated computation
5. Replace `type()` with `isinstance()`
6. Fix quantile_binarization() mutation issue
7. Simplify get_condition() logic
8. Add error handling for missing conditions

### Long Term (Low Priority):
9. Complete docstring placeholders
10. Add missing type hints
11. Fix variable naming issues
12. Fix typos

---

## Testing Recommendations

1. **Unit Tests Needed**:
   - ExpressionSet validation (empty lists, duplicates, shape mismatches)
   - p_values property with None vs. defined values
   - get_condition() with invalid conditions
   - quantile_binarization() doesn't mutate input

2. **Integration Tests**:
   - E-Flux, GIMME, iMAT with real expression data
   - Preprocessing pipeline end-to-end

3. **Edge Cases**:
   - Single condition/single identifier
   - All zero expression values
   - Missing data handling

---

## Statistics

| Category | Count |
|----------|-------|
| Total Lines | 852 |
| Files | 6 |
| Classes | 2 (ExpressionSet, Preprocessing) |
| Functions | 10 |
| Linting Issues | 0 |
| **Total Issues Found** | **13** |

### Issues by Priority:
- 🔴 Critical: 0
- 🟠 High: 3
- 🟡 Medium: 6
- 🟢 Low: 4

---

**Overall Assessment**: The omics module is in **good shape** with clean, well-structured code. The main issues are minor bugs (p_values check, debug print) and code quality improvements (mutation, efficiency). No critical security or mathematical issues found.

**Next Steps**: Address high priority issues first (debug print, p_values bug, validation), then tackle medium priority items (efficiency, mutations).
