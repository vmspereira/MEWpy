[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/mewpy.svg)](https://badge.fury.io/py/mewpy)
[![Documentation Status](https://readthedocs.org/projects/mewpy/badge/?version=latest)](https://mewpy.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/BioSystemsUM/mewpy.svg?branch=master)](https://travis-ci.org/BioSystemsUM/mewpy)

# MEWpy

MEWpy is an integrated Metabolic Engineering Workbench for strain design optimization. It offers methods to explore different classes of constraint-based models (CBM) for:

- **Metabolic Modeling**: Simulate steady-state metabolic models with different formulations (GECKO, ETFL) and kinetic models
- **Strain Design Optimization**: Evolutionary computation-based optimization by knocking out (KO) or over/under expressing (OU) reactions, genes, or enzymes
- **Omics Integration**: Integration of transcriptomics and expression data (eFlux, GIMME, iMAT)
- **Regulatory Networks**: Integrated metabolic-regulatory modeling (PROM, SRFBA, RFBA, CoRegFlux)
- **Community Modeling**: Microbial community simulation and optimization (SteadyCOM, SMETANA)

MEWpy supports [REFRAMED](https://github.com/cdanielmachado/reframed) and [COBRApy](https://opencobra.github.io/cobrapy/) simulation environments. The optimization engine relies on either [inspyred](https://github.com/aarongarrett/inspyred) or [jMetalPy](https://github.com/jMetal/jMetalPy) packages.

## Key Features

### Metabolic and Integrated Modeling
- **Unified Simulator API**: Seamless integration with COBRApy and REFRAMED models
- **GERM Framework**: Generic representation supporting metabolic, regulatory, and integrated models
- **Analysis Methods**: FBA, pFBA, FVA, MOMA, lMOMA, ROOM, gene/reaction deletions
- **Regulatory Integration**: PROM, SRFBA, RFBA, CoRegFlux for integrated metabolic-regulatory analysis
- **Expression Data Integration**: eFlux, GIMME, iMAT with built-in tools for processing omics data
- **Enzyme-Constrained Models**: Support for GECKO and ETFL formulations

### Strain Design Optimization
- **Multi-Objective Optimization**: Evolutionary algorithms for strain design using inspyred or jMetalPy
- **Flexible Design Strategies**: Reaction/gene/enzyme knockouts (KO) and over/under expression (OU)
- **Custom Objectives**: BPCY, product yield, biomass-product coupled yield, and custom fitness functions
- **Problem Types**: Support for metabolic, kinetic, GECKO, ETFL, and community optimization problems

### Community Modeling
- **SteadyCOM**: Simulate microbial communities at steady-state
- **SMETANA**: Species METabolic interaction ANAlysis for studying metabolic interactions
- **Community Optimization**: Strain design for microbial consortia
- **Cross-feeding Analysis**: Identify and optimize metabolic dependencies

### Kinetic Modeling
- **ODE-based Simulation**: Support for kinetic metabolic models
- **Parameter Optimization**: Fit kinetic parameters to experimental data
- **Dynamic Strain Design**: Knockout optimization for kinetic models
- **Hybrid Approaches**: Integration of kinetic and constraint-based models

### Solver Support
- **Commercial Solvers**: CPLEX (recommended for large problems), Gurobi
- **Open-Source Solvers**: SCIP (`pip install pyscipopt`), GLPK
- **Automatic Fallback**: Seamlessly switches between available solvers

## Quick Start

```python
from mewpy.io import read_sbml
from mewpy.simulation import get_simulator, SimulationMethod

# Load a model
model = read_sbml('path/to/model.xml')

# Create simulator
simulator = get_simulator(model)

# Run FBA
result = simulator.simulate(method=SimulationMethod.FBA)
print(f"Objective value: {result.objective_value}")

# Perform FVA
from mewpy.germ.analysis import fva
fva_result = fva(model)
```

For integrated metabolic-regulatory models:

```python
from mewpy.io import read_model, Reader, Engines
from mewpy.germ.analysis import SRFBA, PROM

# Load integrated model
gem_reader = Reader(Engines.MetabolicSBML, 'metabolic_model.xml')
trn_reader = Reader(Engines.BooleanRegulatoryCSV, 'regulatory_network.csv',
                    sep=',', id_col=0, rule_col=2, header=0)
model = read_model(gem_reader, trn_reader)

# Run SRFBA
srfba = SRFBA(model).build()
solution = srfba.optimize()
```

For strain design optimization:

```python
from mewpy.optimization import EA
from mewpy.problems import GKOProblem
from mewpy.simulation import get_simulator

# Load model
model = read_sbml('path/to/model.xml')
simulator = get_simulator(model)

# Define optimization problem
problem = GKOProblem(
    model=model,
    fevaluation=['BPCY'],  # Biomass-Product Coupled Yield
    product='EX_succ_e',    # Target product (e.g., succinate)
    max_candidate_size=10   # Maximum number of gene knockouts
)

# Run evolutionary algorithm
ea = EA(problem, max_generations=50, mp=True)
final_pop = ea.run()

# Display results
ea.dataframe()
```

## Examples

Examples are provided as [jupyter notebooks](https://github.com/vmspereira/MEWpy/tree/master/examples) and as [python scripts](https://github.com/vmspereira/MEWpy/tree/master/examples/scripts).

## Documentation

The package documentation is available at [mewpy.readthedocs.io](https://mewpy.readthedocs.io).

## Installation

### From PyPI

```bash
pip install mewpy
```

### From GitHub

```bash
git clone https://github.com/vmspereira/MEWpy.git
cd MEWpy
pip install -e .
```

### Development Installation

```bash
pip install -e ".[dev]"     # Install with development dependencies
pip install -e ".[test]"    # Install with testing dependencies
pip install -e ".[solvers]" # Install with SCIP solver support
```

### Solver Requirements

MEWpy requires at least one compatible linear programming solver. Supported solvers:

- **[CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio)** (commercial, recommended for large problems)
- **[Gurobi](http://www.gurobi.com)** (commercial)
- **[SCIP](https://www.scipopt.org/)** (open-source, install via `pip install pyscipopt`)
- **[GLPK](https://www.gnu.org/software/glpk/)** (open-source)

For open-source solver installation:

```bash
pip install pyscipopt  # Installs SCIP solver
```

## Cite

If you use MEWpy in your research, please cite:

Vítor Pereira, Fernando Cruz, Miguel Rocha, MEWpy: a computational strain optimization workbench in Python, Bioinformatics, 2021; [https://doi.org/10.1093/bioinformatics/btab013](https://doi.org/10.1093/bioinformatics/btab013)

## Credits and License

Developed by Vítor Pereira and Centre of Biological Engineering, University of Minho. MEWpy received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement number 814408 (2019-2023).

MEWpy is currently maintained by Vítor Pereira.

Released under the GNU General Public License v3.0.
