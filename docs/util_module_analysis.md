# Code Quality Analysis - mewpy.util Module

**Date:** 2025-12-27
**Analyzed by:** Automated code quality review
**Module:** src/mewpy/util/

---

## Executive Summary

Comprehensive analysis of the mewpy.util module identified **50+ code quality issues** across all priority levels:
- **2 CRITICAL issues** requiring immediate attention (logic errors)
- **18+ HIGH priority issues** (bare excepts, print statements, mutable defaults)
- **20+ MEDIUM priority issues** (spelling errors, wildcard imports)
- **10+ LOW priority issues** (missing type hints, TODO comments)

---

## 🔴 CRITICAL PRIORITY ISSUES

### 1. **Incorrect Symbol Definition (Logic Error)**

**Location**: `parsing.py`, line 42

**Issue**:
```python
S_LESS_THAN_EQUAL = "=>"
```

**Problem**:
- The "less than or equal to" operator is defined as "=>" (which typically means "greater than or equal" in some languages)
- This is a logic error that breaks functionality
- Expressions using "<=" will be interpreted incorrectly

**Fix**:
```python
S_LESS_THAN_EQUAL = "<="
```

---

### 2. **Missing Return Statement in Decorator Wrapper**

**Location**: `history.py`, lines 99-114

**Issue**:
```python
def recorder(func: Callable):
    @wraps(func)
    def wrapper(self: Union["Model", "Variable"], value):
        history = self.history
        old_value = getattr(self, func.__name__)
        if old_value != value:
            history.queue_command(undo_func=func, undo_args=(self, old_value), func=func, args=(self, value))
        func(self, value)
    return wrapper
```

**Problem**:
- The wrapper doesn't explicitly return anything
- Will return `None` instead of what the wrapped function might return
- While type hint says `-> None`, this may break expected behavior

**Fix**:
```python
@wraps(func)
def wrapper(self: Union["Model", "Variable"], value):
    history = self.history
    old_value = getattr(self, func.__name__)
    if old_value != value:
        history.queue_command(undo_func=func, undo_args=(self, old_value), func=func, args=(self, value))
    return func(self, value)  # Add explicit return
```

---

## 🟠 HIGH PRIORITY ISSUES

### 3. **Bare Except Clauses (5 instances)**

**Location**: `request.py`, lines 32, 75, 83, 87, 193

**Issue**:
```python
# Line 32
except:
    tokens = entry["uniProtkbId"].split("_")
    name = tokens[0]
    print(f"No gene name for {protein} using uniProtkbId")

# Line 75
except:
    print("No comments")

# Lines 83, 87, 193
except:
    pass
```

**Problem**:
- Bare `except:` catches ALL exceptions including `KeyboardInterrupt`, `SystemExit`
- Makes debugging extremely difficult
- Masks unexpected errors

**Fix**:
```python
# Line 32
except (KeyError, IndexError, TypeError):
    tokens = entry["uniProtkbId"].split("_")
    name = tokens[0]
    print(f"No gene name for {protein} using uniProtkbId")

# Line 75
except KeyError:
    print("No comments")

# Lines 83, 87, 193
except (KeyError, AttributeError, TypeError):
    pass
```

---

### 4. **Print Statements in Library Code (8 instances)**

**Locations**:
- `graph.py`: lines 116, 129
- `parsing.py`: lines 304, 308
- `request.py`: lines 35, 76
- `process.py`: lines 154, 275
- `utilities.py`: line 102

**Issue**:
```python
# graph.py, line 116
print(s[:5])

# request.py, line 35
print(f"No gene name for {protein} using uniProtkbId")

# process.py, line 154
print("nodaemon")
```

**Problem**:
- Library code should never use print statements
- Prevents users from controlling output
- Can't be suppressed or redirected
- No log levels or filtering

**Fix**:
```python
import logging
logger = logging.getLogger(__name__)

# Instead of print(f"No gene name for {protein} using uniProtkbId")
logger.warning(f"No gene name for {protein} using uniProtkbId")

# Instead of print(s[:5])
logger.debug(f"Current state: {s[:5]}")
```

---

### 5. **Mutable Default Arguments (3 instances)**

**Locations**:
- `parsing.py`, line 628
- `graph.py`, lines 40, 136

**Issue**:
```python
# parsing.py, line 628
def __init__(self, true_list=[], variables={}):
    self.true_list = true_list
    self.vars = variables

# graph.py, line 136
def shortest_distance(model, reaction, reactions=None, remove=[]):
    ...

# graph.py, line 40
def create_metabolic_graph(..., remove=[], ...):
    ...
```

**Problem**:
- Mutable default arguments are shared across function calls
- Modifications persist across invocations
- Common Python anti-pattern that causes subtle bugs

**Fix**:
```python
# parsing.py
def __init__(self, true_list=None, variables=None):
    self.true_list = true_list if true_list is not None else []
    self.vars = variables if variables is not None else {}

# graph.py
def shortest_distance(model, reaction, reactions=None, remove=None):
    remove = remove if remove is not None else []

def create_metabolic_graph(..., remove=None, ...):
    remove = remove if remove is not None else []
```

---

### 6. **Generic Exception Types (2 instances)**

**Locations**:
- `parsing.py`, line 96
- `request.py`, line 166

**Issue**:
```python
# parsing.py
raise Exception(f"Unrecognized constant: {type(value).__name__}")

# request.py
raise Exception("zeep library is required.")
```

**Problem**:
- Generic `Exception` is too broad
- Doesn't indicate what went wrong
- Makes error handling difficult

**Fix**:
```python
# parsing.py
raise ValueError(f"Unrecognized constant: {type(value).__name__}")

# request.py
raise ImportError("zeep library is required.")
```

---

### 7. **Broad Exception Handler**

**Location**: `request.py`, line 60

**Issue**:
```python
except Exception:
    pass
```

**Problem**:
- Catching broad `Exception` silently masks errors
- Makes debugging difficult

**Fix**:
```python
except (KeyError, AttributeError):
    # ecNumber not available
    pass
```

---

### 8. **Incomplete Docstring with Placeholder Values**

**Location**: `process.py`, lines 309-320, 333-344

**Issue**:
```python
def get_evaluator(problem, n_mp=cpu_count(), evaluator=ModelConstants.MP_EVALUATOR):
    """Retuns a multiprocessing evaluator

    Args:
        problem: a class implementing an evaluate(candidate) function
        n_mp (int, optional): The number of cpus. Defaults to cpu_count().
        evaluator (str, optional): The evaluator name: options 'ray','dask','spark'.\
            Defaults to ModelConstants.MP_EVALUATOR.

    Returns:
        [type]: [description]
    """
```

**Problem**:
- Return type `[type]` and description `[description]` are placeholder values
- Typo: "Retuns" should be "Returns"

**Fix**:
```python
    """Returns a multiprocessing evaluator

    Args:
        problem: a class implementing an evaluate(candidate) function
        n_mp (int, optional): The number of cpus. Defaults to cpu_count().
        evaluator (str, optional): The evaluator name: options 'ray','dask','spark'.
            Defaults to ModelConstants.MP_EVALUATOR.

    Returns:
        Evaluator: A multiprocessing evaluator instance based on the specified evaluator type
    """
```

---

## 🟡 MEDIUM PRIORITY ISSUES

### 9. **Wildcard Import**

**Location**: `parsing.py`, line 29

**Issue**:
```python
from math import *
```

**Problem**:
- Pollutes the namespace
- Makes it unclear what symbols are available
- Can cause naming conflicts

**Fix**:
```python
from math import (sqrt, sin, cos, tan, log, exp, ceil, floor,
                  log10, log2, radians, degrees, pi, e)
```

---

### 10. **Spelling Errors in Docstrings and Comments**

**Locations and issues**:

1. **parsing.py**, lines 115, 383, 385
```python
# OLD: Operators precedence used to add parentesis when
# NEW: Operators precedence used to add parentheses when
```

2. **parsing.py**, lines 525, 741
```python
# OLD: """Defines a basic arithmetic sintax."""
# NEW: """Defines a basic arithmetic syntax."""
```

3. **request.py**, lines 101, 123, 143
```python
# OLD: def retreive(data, organism=None):
# NEW: def retrieve(data, organism=None):
```

4. **process.py**, line 121
```python
# OLD: When using COBRApy, mewmory resources are not released
# NEW: When using COBRApy, memory resources are not released
```

5. **process.py**, line 309
```python
# OLD: """Retuns a multiprocessing evaluator
# NEW: """Returns a multiprocessing evaluator
```

---

### 11. **Return Type Annotation Inconsistency**

**Location**: `history.py`, lines 65-67

**Issue**:
```python
def __call__(self, *args, **kwargs) -> None:
    return self.queue_command(*args, **kwargs)
```

**Problem**:
- Function is annotated to return `None`
- But explicitly returns the result of `queue_command()`
- Unusual for a function annotated as `-> None`

**Fix**:
Either remove the explicit return or verify if it should return a value:
```python
def __call__(self, *args, **kwargs) -> None:
    self.queue_command(*args, **kwargs)
```

---

### 12. **TODO Comments Without Sufficient Context**

**Location**: `parsing.py`, line 88

**Issue**:
```python
# TODO(odashi): Support other symbols for the imaginary unit than j.
```

**Problem**:
- TODO without creation date or context
- Unclear if this is still relevant

**Fix**:
```python
# TODO(odashi): Support other symbols for the imaginary unit than j.
# Current implementation only supports 'j' for complex numbers
# Consider using a configurable symbol or supporting both 'j' and 'i'
```

---

### 13. **Missing Error Messages in Silent Failures**

**Location**: `request.py`, lines 83, 87

**Issue**:
```python
except:
    pass
```

**Problem**:
- Silent failures without logging make debugging very difficult

**Fix**:
```python
except (KeyError, AttributeError):
    logger.debug("Optional field not available")
```

---

### 14. **Unused Parameter**

**Location**: `graph.py`, line 39-40

**Issue**:
```python
def create_metabolic_graph(
    model, directed=True, carbon=True, reactions=None, remove=[],
    edges_labels=False, biomass=False, metabolites=False
):
```

**Problem**:
- Parameter `biomass` is defined but never used in the function body (lines 39-108)

**Fix**:
Remove unused parameter or implement its functionality

---

## 🟢 LOW PRIORITY ISSUES

### 15. **Missing Type Hints**

**Locations**:
- `graph.py`: All functions lack type hints
- `parsing.py`: Many functions lack type hints
- `utilities.py`: Several functions lack type hints

**Problem**:
- Reduces code readability
- Prevents IDE type checking
- Makes API usage less clear

**Impact**: Low - code works but type hints would improve development experience

---

### 16. **Inconsistent Docstring Formats**

**Issue**:
Different docstring styles are used throughout the module (reST style vs Google style vs PEP 257)

**Fix**:
Standardize on one format across the module

---

### 17. **Variable Name That Shadows Built-in**

**Location**: `parsing.py`, lines 459-464

**Issue**:
```python
l, pl = self.left.to_latex()
r, pr = self.right.to_latex()
```

**Problem**:
- Single-letter `l` (lowercase L) is hard to distinguish from `1` (one)
- Shadows the built-in `list` type if used as `l`

**Fix**:
```python
left_latex, left_prec = self.left.to_latex()
right_latex, right_prec = self.right.to_latex()
```

---

## Summary Statistics

| Category | Count | Files Affected |
|----------|-------|----------------|
| CRITICAL bugs | 2 | 2 files |
| Bare except clauses | 5 | 1 file |
| Print statements | 8 | 5 files |
| Broad exception handlers | 1 | 1 file |
| Mutable default arguments | 3 | 2 files |
| Generic Exception types | 2 | 2 files |
| Spelling errors (docstrings) | 10+ | 3 files |
| Wildcard imports | 1 | 1 file |
| Incomplete docstrings | 2 | 1 file |
| TODO comments | 1+ | 1 file |
| Missing type hints | Many | Multiple files |

---

## Files Requiring Attention

### Critical Priority:
1. **parsing.py** - Logic error in symbol definition (line 42)
2. **history.py** - Missing return in decorator wrapper (line 99-114)

### High Priority:
3. **request.py** - 5 bare except clauses, print statements, generic exceptions
4. **graph.py** - Mutable defaults, print statements
5. **parsing.py** - Mutable defaults, generic exception, wildcard import
6. **process.py** - Print statements, incomplete docstrings
7. **utilities.py** - Print statement

### Medium Priority:
8. All files - Spelling errors in docstrings

---

## Recommended Fix Order

1. **IMMEDIATE**: Fix S_LESS_THAN_EQUAL symbol in parsing.py:42
2. **IMMEDIATE**: Add return statement in history.py:99-114
3. **URGENT**: Replace all 5 bare except clauses in request.py
4. **URGENT**: Fix all mutable default arguments (3 instances)
5. **HIGH**: Replace 8 print statements with logging
6. **HIGH**: Use specific exception types instead of generic Exception
7. **MEDIUM**: Fix spelling errors in docstrings (10+ instances)
8. **MEDIUM**: Replace wildcard import in parsing.py
9. **MEDIUM**: Complete placeholder docstrings in process.py
10. **LOW**: Add type hints where missing
11. **LOW**: Improve TODO comment context

---

## Testing Recommendations

After fixes:
1. Run full test suite: `pytest tests/ -x --tb=short`
2. Run flake8: `flake8 src/mewpy/util/`
3. Run black: `black --check src/mewpy/util/`
4. Run isort: `isort --check-only src/mewpy/util/`
5. Verify no regressions in functionality
