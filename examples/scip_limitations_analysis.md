# SCIP Solver Limitations in MEWpy

## Issues Identified

### 1. **State Machine Constraints** ✅ FIXED
**Problem**: SCIP has a strict state machine where you cannot add constraints after solving without calling `freeTransform()` first.

**Location**: Encountered in `pFBA.build()` where we:
1. Solve FBA to find optimal objective
2. Try to add biomass constraint
3. Add auxiliary variables for minimization

**Error**:
```
Exception: SCIP: method cannot be called at this time in solution process!
```

**Fix Applied**: Create a temporary solver for the initial FBA solve, then create a fresh solver for the pFBA problem (src/mewpy/germ/analysis/pfba.py:63-77)

**Status**: ✅ Fixed

---

### 2. **Repeated Problem Modifications**
**Problem**: During gene/reaction deletion analysis, the same FBA object is reused with many different constraint sets. SCIP handles temporary constraints by:
- Calling `freeTransform()` before each modification
- Changing variable bounds
- Re-solving

**Location**: Used in:
- `single_gene_deletion()` - loops through all genes
- `single_reaction_deletion()` - loops through all reactions
- `fva()` - solves for min/max of each reaction

**Impact**: Works but may be slower than CPLEX/Gurobi due to overhead of freeing and rebuilding the transform.

**Current Implementation**:
```python
# src/mewpy/solvers/pyscipopt_solver.py:119-122
try:
    self.problem.freeTransform()
except:
    pass  # Might not be transformed yet
```

**Status**: ⚠️ Working but suboptimal performance

---

### 3. **LP Solver Internal Errors**
**Problem**: "SCIP: error in LP solver!" during gene deletion analysis

**Possible Causes**:
1. **Numerical instability**: Repeated modifications may accumulate numerical errors
2. **SoPlex issues**: SCIP's default LP solver (SoPlex) may encounter edge cases
3. **Constraint conflicts**: Temporary constraints may create infeasible/unbounded problems

**Location**: Appears intermittently in `single_gene_deletion()` when applying gene knockouts

**Error**:
```
Exception: SCIP: error in LP solver!
```

**Status**: ⚠️ Intermittent, needs investigation

---

## Potential Improvements

### Option 1: Fresh Solver Per Solve (Easy) ✅ RECOMMENDED

Instead of reusing solvers with temporary constraints, create fresh solvers for each deletion:

```python
# In single_gene_deletion loop
for gene in genes:
    # Create fresh FBA instance for this gene
    gene_fba = FBA(model).build()
    solution, status = run_method_and_decode(
        method=gene_fba,
        constraints={**constraints, **gene_constraints}
    )
```

**Pros**:
- Avoids state machine issues
- Cleaner problem structure for each solve
- May be more stable

**Cons**:
- Slight overhead from rebuilding (but SCIP is fast at this)
- Uses more memory temporarily

---

### Option 2: Better Temporary Constraint Handling (Medium)

Improve how temporary constraints are applied/removed:

```python
def _apply_temporary_constraints_optimized(self, constraints):
    """
    Apply temporary constraints more efficiently for SCIP.
    Only free transform once, apply all changes, then optimize.
    """
    # Free transform once
    try:
        self.problem.freeTransform()
    except:
        pass

    # Store original bounds
    temp_constrs = []
    for var_id, bounds in constraints.items():
        if var_id in self._vars:
            orig_lb = self._cached_lower_bounds[var_id]
            orig_ub = self._cached_upper_bounds[var_id]
            temp_constrs.append((var_id, orig_lb, orig_ub))

            # Apply new bounds without calling freeTransform again
            lb, ub = bounds if isinstance(bounds, tuple) else (bounds, bounds)
            var = self._vars[var_id]
            if lb is not None:
                self.problem.chgVarLb(var, lb if lb != -inf else -self.problem.infinity())
            if ub is not None:
                self.problem.chgVarUb(var, ub if ub != inf else self.problem.infinity())

    return temp_constrs
```

**Pros**:
- More efficient than current approach
- Reduces freeTransform() calls

**Cons**:
- Requires changes to solver interface
- Still has some overhead

---

### Option 3: Problem Pooling (Hard)

Maintain a pool of solver instances to avoid rebuild overhead:

```python
class SCIPSolverPool:
    def __init__(self, base_model, pool_size=4):
        self.solvers = [PySCIPOptSolver(base_model) for _ in range(pool_size)]
        self.available = self.solvers.copy()

    def get_solver(self):
        if self.available:
            return self.available.pop()
        return PySCIPOptSolver(self.base_model)  # Create new if pool empty

    def return_solver(self, solver):
        solver.problem.freeTransform()  # Reset state
        self.available.append(solver)
```

**Pros**:
- Amortizes solver creation cost
- Good for analyses with many similar problems

**Cons**:
- Complex to implement
- Memory overhead
- May not help much given SCIP's fast problem building

---

### Option 4: SCIP Parameters Tuning (Easy) 🔧 RECOMMENDED

Add better default parameters for stability:

```python
def __init__(self, model=None):
    # ... existing code ...

    # Tuning for numerical stability
    self.problem.setParam("numerics/epsilon", 1e-9)
    self.problem.setParam("numerics/sumepsilon", 1e-6)
    self.problem.setParam("numerics/feastol", 1e-6)
    self.problem.setParam("numerics/lpfeastolfactor", 1.0)

    # Disable presolving for repeated solves (faster)
    # self.problem.setParam("presolving/maxrounds", 0)  # Optional

    # Use more stable LP solver settings
    self.problem.setParam("lp/threads", 1)  # Single-threaded for consistency
```

**Pros**:
- Easy to implement
- May fix numerical issues
- No API changes needed

**Cons**:
- May not solve all issues
- Could slightly impact performance

---

## Recommendations

### Short Term (Do Now):
1. ✅ **pFBA fix is already applied** - use fresh solver approach
2. 🔧 **Add parameter tuning** for numerical stability (Option 4)
3. 📝 **Document SCIP limitations** in code and user docs

### Medium Term (Consider):
1. **Implement Option 1** for deletion analyses - create fresh FBA per deletion
2. **Add solver selection warnings** - warn users if analyzing large models with SCIP
3. **Benchmark SCIP vs CPLEX/Gurobi** - quantify performance differences

### Long Term (Future):
1. **Optimize temporary constraints** (Option 2)
2. **Consider solver pooling** for very large analyses (Option 3)

---

## Comparison with Other Solvers

| Feature | CPLEX | Gurobi | SCIP | OptLang |
|---------|-------|--------|------|---------|
| Add constraints after solve | ✅ Easy | ✅ Easy | ⚠️ Need freeTransform() | ✅ Easy |
| Temporary constraints | ✅ Fast | ✅ Fast | ⚠️ Slower | ✅ Fast |
| Numerical stability | ✅✅ Excellent | ✅✅ Excellent | ✅ Good | Depends on backend |
| Open source | ❌ No | ❌ No | ✅ Yes | ✅ Yes |
| Performance | ✅✅ Best | ✅✅ Best | ✅ Good | Depends on backend |

---

## Testing Notes

- ✅ GERM_Models.ipynb runs successfully with SCIP
- ⚠️ GERM_Models_analysis.ipynb has intermittent LP errors in gene deletion
- ✅ pFBA works after fix
- ✅ FBA, RFBA, SRFBA work correctly

---

## Conclusion

SCIP is a viable open-source alternative to CPLEX/Gurobi for MEWpy, but requires:
1. Awareness of state machine limitations (mostly handled)
2. Parameter tuning for numerical stability
3. Potentially different coding patterns for repeated solves

For most use cases, SCIP will work fine. For large-scale gene deletion or FVA analyses, commercial solvers may be faster.
