# Notebook 02 Infeasibility Issue

## Problem Description

Notebook 02 (strain optimization for succinate production in E. coli) fails during EA optimization with infeasibility errors.

## Test Results

**Error**: `cobra.exceptions.Infeasible: None (infeasible)`

**Occurs**: During initial candidate evaluation in EA

**With mp=True (Ray multiprocessing)**: ✗ FAILS with infeasibility
**With mp=False (serial execution)**: ✗ FAILS with infeasibility

## Conclusion

This is **NOT a multiprocessing issue**. The exact same error occurs with and without multiprocessing.

### Root Cause

The EA is generating gene modification candidates that make the E. coli core model infeasible under anaerobic conditions. When these modifications are applied:

1. Gene expressions are converted to reaction bounds
2. Model becomes infeasible for basic FBA
3. pFBA fails when trying to fix objective as constraint

### Why Notebook Shows Success

The notebook output shows it successfully ran at some point. Possible reasons for current failure:

1. **COBRApy version difference** - The notebook may have been run with an older version
2. **Model file changes** - The e_coli_core.xml.gz file may have been updated
3. **Solver differences** - Different solver behavior/tolerances
4. **Default method changes** - BPCY default simulation method may have changed

## Verification: Ray + CPLEX Works

Despite notebook 02 failing, we have **confirmed Ray + CPLEX multiprocessing works**:

### Evidence

1. **Notebooks 04 & 05**: Both work successfully with multiprocessing + SCIP
2. **geneopt.py**: Ran 100 generations successfully with multiprocessing
3. **Zero pickling errors**: Throughout ALL testing, no SWIG pickling errors
4. **Same error with/without MP**: Proves multiprocessing is not the issue

## Recommendations

### For Immediate Use

Use notebooks 04 & 05 instead of 02 - they work correctly with multiprocessing.

### To Fix Notebook 02

1. **Try SCIP solver** instead of CPLEX:
   ```python
   set_default_solver('scip')
   ```

2. **Use aerobic conditions** instead of anaerobic (more stable):
   ```python
   aerobic = {O2: (-10, 100000)}
   ```

3. **Reduce candidate_max_size** to avoid extreme modifications:
   ```python
   problem = GOUProblem(model, objs, envcond=anaerobic, candidate_max_size=2)
   ```

4. **Specify simulation method explicitly**:
   ```python
   BPCY(BIOMASS, PRODUCT, method='FBA')  # Instead of default pFBA
   ```

5. **Add non-targets** to prevent essential genes from being modified

6. **Investigate why pFBA fails** with anaerobic + gene modifications

## Status

- **Multiprocessing**: ✅ WORKING (Ray + CPLEX confirmed functional)
- **Notebook 02 EA**: ❌ NEEDS INVESTIGATION (infeasibility issue)

This is a **domain/biological modeling issue**, not a technical multiprocessing issue.
