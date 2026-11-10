## Options

MEWpy options are set as arguments of a problem or of `EA`, or globally in the `EAConstants` and `ModelConstants` classes of _mewpy.util.constants_. Set global options before creating the problem and the EA.

The examples below use the E. coli core model, `tests/data/e_coli_core.xml` in the MEWpy repository:

```python
from reframed.io.sbml import load_cbmodel
from mewpy.optimization.evaluation import BPCY, TargetFlux

model = load_cbmodel('tests/data/e_coli_core.xml')
objectives = [BPCY('R_Biomass_Ecoli_core', 'R_EX_succ_e'), TargetFlux('R_EX_succ_e')]
```

**Number of processors for parallel solutions evaluation**

By default, `EA` evaluates solutions in parallel (`mp=True`) using half of the CPUs reported by `multiprocessing.cpu_count()`, and at least one. To set the number of parallel processes, change the `NUM_CPUS` constant in `mewpy.util.constants`. A value of 0 or less uses the default.

```python
from mewpy.util.constants import EAConstants
# uses 32 parallel processes
EAConstants.NUM_CPUS = 32
```

`ModelConstants.MP_EVALUATOR` selects how solutions are distributed: `'ray'` (the default), `'mp'` (a process pool), `'nodaemon'` (a process pool whose workers can start their own processes), `'dask'` or `'spark'`. The `'dask'` evaluator does not use `NUM_CPUS`. If the package an evaluator needs is not installed, MEWpy uses `'mp'`. To evaluate solutions in a single process, run `EA(problem, mp=False)`.

```python
from mewpy.util.constants import ModelConstants
ModelConstants.MP_EVALUATOR = 'mp'
```

**Over and under expression folds.**

Over and under expression problems choose the fold of each modification from a list of levels. The default, `EAConstants.LEVELS`, is `[0, 1/32, 1/16, 1/8, 1/4, 1/2, 2, 4, 8, 16, 32]`, where 0 is a deletion. `OptRamProblem` uses `EAConstants.OPTRAM_LEVELS`, which replaces 0 with 0.001, and `MediumProblem` uses the levels from 0 to 10 in steps of 0.1. The levels may be defined in the problem definition, for example, in a genes' under expression problem:

```python
# Define the regulation folds for under expression or deletion.
# 0 for deletion.
levels = [1/8, 1/4, 1/2, 0]

from mewpy.problems import GOUProblem
problem = GOUProblem(model, objectives, levels=levels)
```

**Number of modifications**

The minimum and the maximum number of modifications may be defined directly in the problem definition. For example, to allow a maximum of 6 gene deletions:

```python
from mewpy.problems import GKOProblem
problem = GKOProblem(model, objectives, candidate_max_size=6)
```

Likewise, the minimum number of modifications may be explicitly defined:

```python
from mewpy.problems import GKOProblem
problem = GKOProblem(model, objectives, candidate_min_size=4, candidate_max_size=6)
```

The default minimum and maximum numbers of modifications are 1 and 10 (`EAConstants.MIN_SOLUTION_SIZE` and `EAConstants.MAX_SOLUTION_SIZE`). `MediumProblem` sets the default maximum to the number of uptake reactions. When the minimum and the maximum are equal, the EA keeps all solutions at that size, though the simplification of solutions described below may then remove modifications.

**Optimization algorithm**

MEWpy resorts to the jMetalPy and inspyred packages to evolve modification solutions. If both packages are installed, MEWpy uses jMetalPy by default, running the Non-dominated Sorting Genetic Algorithm (NSGA-II) for multi-objective problems and a Genetic Algorithm (GA) for single objective problems. To use inspyred instead, add the following to your script:

```python
from mewpy.optimization import set_default_engine
set_default_engine('inspyred')
```

MEWpy also allows to define the optimization algorithm when configuring the EA for the addressed optimization problem:

```python
from mewpy.optimization import EA
from mewpy.problems import RKOProblem

problem = RKOProblem(model, objectives)
ea = EA(problem, max_generations=100, algorithm='NSGAIII')
```

The algorithm names are `'GA'`, `'SA'` (Simulated Annealing), `'NSGAII'`, `'SPEA2'` and `'NSGAIII'`. inspyred provides GA, SA and NSGA-II; `'SPEA2'` and `'NSGAIII'` always run on jMetalPy. An unknown name runs the preferred algorithm, NSGA-II unless changed with `set_preferred_EA`, which also switches the default engine if it does not provide that algorithm:

```python
from mewpy.optimization import set_preferred_EA
set_preferred_EA('SPEA2')
```

The number of objectives limits the choice:

* A single objective problem runs SA when `algorithm='SA'` and GA otherwise.
* On a multi-objective problem, inspyred always runs NSGA-II. jMetalPy runs NSGA-II, SPEA2 or NSGA-III as requested and, if given GA or SA, runs SPEA2 for two objectives and NSGA-III for three or more.

To solve a multi-objective problem with a single objective algorithm, combine the objectives into one with the `AggregatedSum` evaluation function, for example `AggregatedSum(objectives, tradeoffs=[0.5, 0.5])`.

**Population size, number of generations and seed**

`EA` takes the population size (default 100), the number of generations (default 100) and a seed for the random number generators, for reproducible runs.

```python
ea = EA(problem, max_generations=50, population_size=50, seed=42)
```

`mewpy.optimization.settings.set_default_population_size` changes the default population size.

**Seeding an EA with an initial population.**

The EAs may be seeded with a list of solutions, to guide the optimization or to give it a push start.

For deletion problems, the initial population is a list of solutions, each a list of modification targets. For example, for a reaction knockout problem:

```python
from mewpy.problems import RKOProblem
problem = RKOProblem(model, objectives)

init_pop = [['R_PFK', 'R_ACKr'],
            ['R_PGI', 'R_LDH_D', 'R_PYK']]

ea = EA(problem, initial_population=init_pop)
```

Each `'R_...'` is a reaction in the modification target list. Unless a list is given with the `target` parameter, the problem builds it from the model; for reaction knockouts it excludes essential, exchange and transport reactions, and the reactions in the `non_target` parameter. The list is available from the problem instance:

```python
problem.target_list
```

For over/under expression problems, the initial population is a list of dictionaries:

```python
from mewpy.problems import ROUProblem
problem = ROUProblem(model, objectives)

init_pop = [{'R_PFK': 2, 'R_ACKr': 0},
            {'R_PGI': 8, 'R_LDH_D': 0.5, 'R_PYK': 0}]

ea = EA(problem, initial_population=init_pop)
```

where each item is of the form `modification_target: fold_level`, and each fold level is one of the problem's `levels`. In reaction and gene problems such as these, identifiers may be given with or without their SBML prefix (`'PFK'` or `'R_PFK'`). A solution with an identifier that is not in the target list, or with a level that is not in `levels`, is skipped with a warning by jMetalPy and raises a `ValueError` in inspyred.

**Simplification of solutions**

By default, MEWpy simplifies the final set of solutions by removing each modification whose removal leaves all the objective values unchanged (within 1e-6). This step, which evaluates each solution once more plus once per modification, may be skipped by setting the simplify flag to false when running the EA:

```python
ea.run(simplify=False)
```

**Other global options**

* `EAConstants.KILL_DUMP = True` makes an EA interrupted with Ctrl+C write its current population to a `mewpy-dump-<date>-<time>.json` file in the working directory.
* MEWpy logs each distinct evaluation error once. `EAConstants.DEBUG = True` also issues a warning for every solution that cannot be evaluated.
