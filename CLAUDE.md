# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MEWpy is a Metabolic Engineering Workbench in Python for strain design optimization. It provides methods to explore constraint-based models (CBM) including:
- Simulating single organisms with steady-state metabolic models (GECKO, ETFL) and kinetic models
- Evolutionary computation-based strain design optimization (KO/OU of reactions, genes, enzymes)
- Omics data integration (eFlux, GIMME, iMAT)
- Regulatory network integration
- Microbial community modeling (SteadyCOM, SMETANA)

Supports REFRAMED and COBRApy simulation environments. Optimization relies on inspyred or jMetalPy packages.

## Development Commands

### Environment Setup
```bash
# Install in development mode with all dev dependencies
pip install -e ".[dev]"

# Install with test dependencies
pip install -e ".[test]"

# Install with optional SCIP solver
pip install -e ".[solvers]"
```

### Testing
```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_a_simulator.py

# Run tests with coverage
pytest --cov=mewpy --cov-report=html

# Run tests using tox (multiple Python versions)
tox
```

### Linting and Formatting
```bash
# Check code style with flake8
flake8 src/mewpy

# Format code with black
black src/mewpy

# Sort imports with isort
isort src/mewpy

# Run all linting (flake8 uses .flake8 config)
# Max line length: 125 for flake8, 120 for black/isort
# Black and flake8 configured to be compatible (E203, W503 ignored)
```

### Building
```bash
# Build package
python -m build

# Install from source
pip install .
```

## Architecture Overview

### Package Structure

The codebase is organized under `src/mewpy/` with these main modules:

**Core Framework:**
- `germ/` - GEneric Representation of Models (metabolic, regulatory, integrated models)
- `simulation/` - Phenotype simulation framework with adapters for COBRA/REFRAMED/GERM
- `optimization/` - Evolutionary algorithms for strain design (inspyred/jmetal backends)
- `problems/` - Optimization problem definitions (KO/OU for reactions, genes, enzymes)
- `solvers/` - LP solver interfaces (CPLEX, Gurobi, SCIP, OptLang) and ODE solvers
- `io/` - Model I/O operations (SBML, JSON, CSV) using builder/director pattern

**Domain-Specific:**
- `model/` - Specialized model types (GECKO, kinetic, SMoment)
- `cobra/` - COBRApy integration utilities
- `com/` - Community modeling (SteadyCOM, SMETANA)
- `omics/` - Omics data integration (eFlux, GIMME, iMAT)
- `util/` - Utilities (parsing, constants, history management)
- `visualization/` - Plotting and visualization

### Key Design Patterns

#### 1. Simulator Pattern (Adapter + Strategy)
Location: `src/mewpy/simulation/`

Abstracts different metabolic modeling platforms:
- `get_simulator()` factory function selects appropriate simulator
- `Simulator` base class with implementations for COBRA, REFRAMED, GERM, Hybrid, Kinetic
- `SimulationMethod` enum defines FBA variants (FBA, pFBA, MOMA, lMOMA, ROOM)

```python
from mewpy import get_simulator
from mewpy.simulation import SimulationMethod

simul = get_simulator(model)  # Auto-detects model type
result = simul.simulate(method=SimulationMethod.FBA, constraints={...})
```

#### 2. GERM - Generic Model Representation
Location: `src/mewpy/germ/`

Uses metaclass programming for dynamic model composition:
- `MetaModel` metaclass creates composite model classes on-the-fly
- Supports metabolic-only, regulatory-only, or integrated models
- Polymorphic constructors: `Model.from_types(['metabolic', 'regulatory'])`
- Type checkers: `model.is_metabolic()`, `model.is_regulatory()`

**Sub-modules:**
- `models/` - Base model classes with dynamic typing
- `variables/` - Variable types (genes, reactions, metabolites, regulators)
- `algebra/` - Expression trees and symbolic computation
- `lp/` - Linear programming problem construction
- `analysis/` - Analysis methods (FBA, pFBA, RFBA, SRFBA, PROM, CoRegFlux)

#### 3. Optimization Framework (Strategy + Template Method)
Location: `src/mewpy/optimization/`

- `AbstractEA` defines template for evolutionary algorithms
- `EA()` factory selects engine (inspyred or jmetal)
- `evaluation/` contains pluggable fitness functions (BPCY, WYIELD, TargetFlux)
- Separate implementations in `inspyred/` and `jmetal/` subdirectories

#### 4. Problem Hierarchy (Template Method)
Location: `src/mewpy/problems/`

- `AbstractProblem` base class with template methods:
  - `generator()` - Create random solutions
  - `encode()` / `decode()` - EA representation conversion
  - `solution_to_constraints()` - Map to metabolic constraints
- Concrete implementations:
  - `RKOProblem` / `ROUProblem` - Reaction knockout/over-under expression
  - `GKOProblem` / `GOUProblem` - Gene-based optimization
  - `GeckoKOProblem` / `GeckoOUProblem` - Enzyme-constrained models
  - `ETFLGKOProblem` / `ETFLGOUProblem` - ETFL models
  - `KineticKOProblem` / `KineticOUProblem` - Kinetic models
  - `CommunityKOProblem` - Community optimization
  - `OptORFProblem` / `OptRAMProblem` - Regulatory optimization

### Component Interaction Flow

**Typical Optimization Workflow:**
```
Model (COBRA/REFRAMED/GERM)
    ↓
get_simulator() [Factory selects appropriate adapter]
    ↓
Problem (defines solution space + evaluation functions)
    ↓
EA() [Factory selects inspyred/jmetal engine]
    ↓
Evolutionary Loop:
  - generator() creates candidates
  - decode() → solution_to_constraints()
  - simulate() on modified model
  - EvaluationFunction.get_fitness()
    ↓
Solutions (values, fitness, constraints)
```

**GERM Integrated Analysis:**
```
MetabolicModel + RegulatoryModel
    ↓
Model.from_types(['metabolic', 'regulatory'])
    ↓
MetabolicRegulatoryModel (dynamic class created by metaclass)
    ↓
Analysis methods: RFBA, SRFBA, PROM, CoRegFlux
    ↓
LinearProblem construction
    ↓
Solver (CPLEX/Gurobi/SCIP)
    ↓
Solution with fluxes + regulatory states
```

### I/O Architecture
Location: `src/mewpy/io/`

Uses Builder + Director pattern:
- `Reader` / `Writer` classes for model serialization
- `Director` orchestrates multiple readers/writers
- `engines/` contains format-specific implementations

**Supported Formats:**
- SBML (metabolic FBC + regulatory QUAL plugins)
- JSON (GERM native format)
- CSV (regulatory networks)
- Direct COBRA/REFRAMED model import

```python
from mewpy.io import read_sbml, read_model, write_model

model = read_sbml('model.xml')  # Load from SBML
write_model(model, 'output.json')  # Save to JSON
```

## Important Notes

### Solver Requirements
MEWpy requires at least one LP solver:
- CPLEX (commercial, preferred for large problems)
- Gurobi (commercial)
- SCIP (open-source, install via `pip install pyscipopt`)

The default solver is configured globally with automatic fallback.

### Model Compatibility
- GERM provides unified interface across model types
- `get_simulator()` automatically detects and wraps COBRA/REFRAMED models
- Models notify attached simulators of changes (observer pattern)
- History management available via `util/history.py` for undo/redo

### Evaluation Functions
Custom fitness functions should extend `EvaluationFunction`:
- Override `get_fitness()` method
- Return tuple of (fitness_values, is_maximization)
- Can use multiple objectives for multi-objective optimization

### Test Organization
Tests follow alphabetical ordering to manage dependencies:
- `test_a_simulator.py` - Simulation framework tests (run first)
- `test_b_problem.py` - Problem definition tests
- `test_c_optimization.py` - Optimization tests
- `test_d_models.py` - Model tests
- `test_e_germ_problem.py` - GERM-specific tests
- `test_f_omics.py` - Omics integration tests
- `test_g_com.py` - Community modeling tests
- `test_h_kin.py` - Kinetic modeling tests

### Code Style
- Line length: 120 (black/isort), 125 (flake8)
- Black formatting is authoritative
- Flake8 configured to be compatible with black (E203, W503 ignored)
- Star imports allowed in `__init__.py` for API exposure (F401, F403, F405)
- Bare except allowed (E722) - review carefully when modifying
