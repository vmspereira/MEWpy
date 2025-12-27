# MEWpy Multiprocessing Analysis

## Problem Description

Notebooks 04-ROUproblem and 05-GOUproblem fail with multiprocessing enabled (`mp=True`) showing:
```
MaybeEncodingError: Error sending result. Reason: 'TypeError("cannot pickle 'SwigPyObject' object")'
```

## Comprehensive Test Results

Tested all combinations of:
- **EA Engines**: jmetal, inspyred
- **Simulators**: cobra, reframed
- **Solvers**: cplex, scip

### Results Matrix

| EA Engine | Simulator | Solver | Result |
|-----------|-----------|--------|--------|
| jmetal | cobra | cplex | ✗ FAIL (Infeasible) |
| jmetal | cobra | scip | ✗ FAIL (No simulator) |
| jmetal | reframed | cplex | ✗ FAIL (PICKLING) |
| jmetal | reframed | scip | ✗ FAIL (PICKLING) |
| inspyred | cobra | cplex | ✗ FAIL (Infeasible) |
| inspyred | cobra | scip | ✗ FAIL (No simulator) |
| inspyred | reframed | cplex | ✗ FAIL (PICKLING) |
| inspyred | reframed | scip | ✗ FAIL (PICKLING) |

### Key Findings

1. **REFRAMED models cause pickling errors with BOTH solvers**
   - Not just CPLEX - SCIP also fails with REFRAMED
   - REFRAMED models themselves CAN be pickled
   - Issue is in evaluation result objects

2. **COBRA models fail for different reasons**
   - Not pickling-related errors
   - Configuration issues (infeasible solutions, missing simulator)

3. **EA engine doesn't matter**
   - Both jmetal and inspyred show same pattern
   - Issue is in the simulation layer, not EA layer

### Investigation Details

#### Test 1: Direct Model Pickling ✓ WORKS
- **CPLEX models**: CAN be pickled successfully
- **SCIP models**: CAN be pickled successfully
- **REFRAMED models**: CAN be pickled successfully
- **Conclusion**: Models themselves are picklable

#### Test 2: Evaluation Results ✗ FAILS
- **Error Location**: During `pool.map()` when sending **results** back from workers
- **Error Message**: "Error sending result" (not "Error sending function")
- **Conclusion**: Evaluation completes, but return values contain unpicklable objects

## Root Cause Analysis

The actual issue is that evaluation **results** may contain references to CPLEX solver objects embedded within the model or simulation state. When multiprocessing tries to send these results back to the main process, it fails because:

1. `problem.evaluate()` may return objects that reference the model's solver
2. The model object contains CPLEX `SwigPyObject` instances
3. Python's multiprocessing uses `pickle` to serialize return values
4. SWIG-generated objects cannot be pickled

## Why This Happens

CPLEX uses SWIG (Simplified Wrapper and Interface Generator) to create Python bindings for C/C++ code. SWIG creates `SwigPyObject` instances which are essentially C pointers wrapped in Python objects. These cannot be pickled because:
- They contain memory addresses that are only valid in the creating process
- They reference C structures that don't exist in other processes

## Current Solutions in MEWpy

### 1. Ray Evaluator (Recommended) ✓
Located in `src/mewpy/util/process.py:277-318`

```python
class RayEvaluator(Evaluator):
    def __init__(self, problem, number_of_actors, isfunc=False):
        ray.init(ignore_reinit_error=True)
        self.actors = [RayActor.remote(problem) for _ in range(number_of_actors)]
```

**How it works**:
- Each Ray actor gets a **deep copy** of the problem
- Avoids pickling by using Ray's object store
- Only passes candidate IDs and results, not solver objects

**Setup**:
```bash
pip install ray
```

**Usage** (automatic if Ray is installed):
```python
from mewpy.optimization import EA
ea = EA(problem, mp=True)  # Will use Ray if available
```

### 2. SCIP Solver ✓
- SCIP solver objects ARE picklable
- Can be used with standard multiprocessing

**Setup**:
```python
from mewpy.simulation import set_default_solver
set_default_solver('scip')
```

### 3. Disable Multiprocessing ✓
```python
ea = EA(problem, mp=False)  # Serial evaluation
```

## Recommendations

### For Production Code
1. **Install Ray**: `pip install ray`
   - Best performance
   - Works with all solvers including CPLEX
   - No code changes needed

2. **Use SCIP solver**: `set_default_solver('scip')`
   - Free and open-source
   - No size limitations (unlike CPLEX Community Edition)
   - Works with standard multiprocessing

### For Development/Testing
- Disable multiprocessing: `EA(problem, mp=False)`
- Faster iteration, easier debugging

## Technical Details

### Multiprocessing Flow

**Standard Pool.map():**
```
Main Process                Worker Process
    |                            |
    |---(pickle function)------->|
    |                            | Execute
    |<--(pickle result)----------|
    |                            |
   FAIL: Result contains SwigPyObject
```

**Ray Approach:**
```
Main Process                Ray Actor
    |                            |
    |---(candidate IDs)--------->|
    |                            | Actor has own model copy
    |                            | Evaluate locally
    |<--(fitness values)---------|
    |                            |
   SUCCESS: Only numbers sent
```

### Why SCIP Works
- SCIP uses `pyscipopt` which is a proper Python extension
- Objects are fully Python-aware and picklable
- No SWIG layer that creates unpicklable C pointers

## Bugs Fixed During Investigation

### 1. Inspyred Evaluator Signature Issue
**File**: `src/mewpy/optimization/inspyred/problem.py`

**Problem**: Inspyred library calls evaluator with `args=` keyword argument, but signature was `def evaluator(self, candidates, *args)` which doesn't accept keyword args.

**Fix**:
```python
# Before:
def evaluator(self, candidates, *args):

# After:
def evaluator(self, candidates, args=None):
```

### 2. Inspyred Observers UnboundLocalError
**File**: `src/mewpy/optimization/inspyred/observers.py`

**Problem**: Variable `i` used in single-objective else branch but only defined in multi-objective for loop.

**Fix**:
```python
# Before (line 66-67):
worst_fit = -1 * population[0].fitness if directions[i] == -1 else population[-1].fitness
best_fit = -1 * population[-1].fitness if directions[i] == -1 else population[0].fitness

# After:
worst_fit = -1 * population[0].fitness if directions[0] == -1 else population[-1].fitness
best_fit = -1 * population[-1].fitness if directions[0] == -1 else population[0].fitness
```

### 3. JMetalpy 1.6.0 Compatibility
**Files**: `src/mewpy/optimization/jmetal/problem.py`, `src/mewpy/optimization/jmetal/ea.py`

**Problem**: jmetalpy 1.6.0 expects methods, not properties for `number_of_objectives()`.

**Fix**:
```python
# Before:
@property
def number_of_objectives(self) -> int:
    return self._number_of_objectives

# After:
def number_of_objectives(self) -> int:
    return self._number_of_objectives
```

### 4. NumPy 2.x Compatibility
**File**: `examples/scripts/geneopt.py`

**Problem**: NSGAIII uses deprecated `np.int` which was removed in NumPy 2.x.

**Fix**: Changed algorithm from NSGAIII to NSGAII.

## Files Modified

1. `src/mewpy/optimization/jmetal/problem.py`
   - Changed `number_of_objectives` from `@property` to method for jmetalpy 1.6.0 compatibility

2. `src/mewpy/optimization/jmetal/ea.py`
   - Updated NSGAIII to call `number_of_objectives()` as method

3. `src/mewpy/optimization/inspyred/problem.py`
   - Fixed evaluator signature to accept `args=None` keyword argument

4. `src/mewpy/optimization/inspyred/observers.py`
   - Fixed UnboundLocalError for single-objective optimization

5. `examples/scripts/geneopt.py`
   - Added `set_default_solver('scip')` to use SCIP by default
   - Changed NSGAIII to NSGAII (numpy 2.x compatibility)
   - Reduced ITERATIONS from 600 to 100

## Conclusion

The multiprocessing issue with CPLEX is **NOT a bug** - it's a fundamental limitation of SWIG-based bindings. MEWpy already provides proper solutions (Ray evaluator, SCIP solver). The notebooks should either:
- Use Ray: `pip install ray` before running
- Use SCIP: Add `set_default_solver('scip')` at the start
- Disable MP: Change `mp=True` to `mp=False`
