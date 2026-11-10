# Installation

MEWpy requires Python 3.9 or later and at least one linear programming solver.

**From PyPI**

```bash
pip install mewpy
```

**From the GitHub repository**

```bash
git clone https://github.com/vmspereira/MEWpy.git
cd MEWpy
pip install -e .
```

**Optional dependencies**

```bash
pip install "mewpy[solvers]"   # the SCIP solver (pyscipopt)
pip install "mewpy[germ]"      # CoRegFlux, PROM and expression preprocessing (scikit-learn, scipy)
pip install "mewpy[examples]"  # packages some example notebooks use
pip install "mewpy[docs]"      # building this documentation
pip install "mewpy[test]"      # running the test suite
pip install "mewpy[dev]"       # development tools
```

**Solvers**

MEWpy uses the first available of these solvers, in this order. `mewpy.solvers.set_default_solver()` selects another one.

- [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) (commercial, recommended for large problems). `pip install cplex` installs the Community Edition, which is limited to 1000 variables and 1000 constraints.
- [Gurobi](https://www.gurobi.com) (commercial), through `gurobipy`.
- [SCIP](https://www.scipopt.org/) (open source), with `pip install pyscipopt`.
- [GLPK](https://www.gnu.org/software/glpk/) (open source), through optlang, which COBRApy installs.
