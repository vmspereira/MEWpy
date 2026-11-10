[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI version](https://badge.fury.io/py/mewpy.svg)](https://badge.fury.io/py/mewpy)
[![Documentation Status](https://readthedocs.org/projects/mewpy/badge/?version=latest)](https://mewpy.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/vmspereira/MEWpy/actions/workflows/ci.yml/badge.svg)](https://github.com/vmspereira/MEWpy/actions/workflows/ci.yml)

# MEWpy

MEWpy is an integrated Metabolic Engineering Workbench for strain design optimization. It offers methods to explore different classes of constraint-based models (CBM) for:

- **Metabolic Modeling**: Simulate steady-state metabolic models with different formulations (GECKO, ETFL) and kinetic models, and batch cultures with dynamic FBA
- **Strain Design Optimization**: Evolutionary computation-based optimization by knocking out (KO) or over/under expressing (OU) reactions, genes, or enzymes
- **Omics Integration**: Integration of transcriptomics and expression data (eFlux, GIMME, iMAT)
- **Regulatory Networks**: Integrated metabolic-regulatory modeling (PROM, SRFBA, RFBA, dynamic RFBA, CoRegFlux)
- **Community Modeling**: Microbial community simulation and optimization (SteadyCom, SMETANA, community dynamic FBA)

MEWpy supports [REFRAMED](https://github.com/cdanielmachado/reframed) and [COBRApy](https://opencobra.github.io/cobrapy/) simulation environments. The optimization engine relies on either [inspyred](https://github.com/aarongarrett/inspyred) or [jMetalPy](https://github.com/jMetal/jMetalPy) packages.

## Key Features

### Metabolic and Integrated Modeling
- **Unified Simulator API**: The same code runs on COBRApy and REFRAMED models; reaction, metabolite and gene ids can be given with or without the SBML prefix (`R_`, `M_`, `G_`)
- **GERM Framework**: Generic representation supporting metabolic, regulatory, and integrated models
- **Analysis Methods**: FBA, pFBA, FVA, MOMA, lMOMA, ROOM, gene/reaction deletions
- **Dynamic FBA**: Batch cultures of a single organism (`DynamicFBA`) or of a community on a shared medium (`CommunityDynamicFBA`)
- **Regulatory Integration**: PROM, SRFBA, RFBA, dynamic RFBA (`DynamicRFBA`), CoRegFlux for integrated metabolic-regulatory analysis
- **Expression Data Integration**: eFlux, GIMME, iMAT with built-in tools for processing omics data
- **Enzyme-Constrained Models**: Support for GECKO and ETFL formulations

### Strain Design Optimization
- **Multi-Objective Optimization**: Evolutionary algorithms for strain design using inspyred or jMetalPy
- **Flexible Design Strategies**: Reaction/gene/enzyme knockouts (KO) and over/under expression (OU)
- **Objectives**: BPCY, WYIELD, target fluxes, batch titer, yield and productivity (dynamic FBA), and custom fitness functions
- **Problem Types**: Support for metabolic, kinetic, GECKO, ETFL, and community optimization problems

### Community Modeling
- **SteadyCom**: Simulate microbial communities at steady-state
- **Community dynamic FBA**: Members with their own biomass growing on a shared medium
- **SMETANA**: Species METabolic interaction ANAlysis for studying metabolic interactions
- **Community Optimization**: Strain design for microbial consortia
- **Cross-feeding Analysis**: Identify and optimize metabolic dependencies

### Kinetic Modeling
- **ODE-based Simulation**: Support for kinetic metabolic models
- **Strain Design**: Knockout and over/under expression problems for kinetic models
- **Hybrid Approaches**: Integration of kinetic and constraint-based models

### Solver Support
- **Commercial Solvers**: CPLEX (recommended for large problems), Gurobi
- **Open-Source Solvers**: SCIP (`pip install pyscipopt`), GLPK
- **Default Solver**: The first available of CPLEX, Gurobi, SCIP and GLPK (through optlang); `set_default_solver()` picks another

## Quick Start

```python
from mewpy.io import read_sbml
from mewpy.simulation import get_simulator, SimulationMethod

# Load a model (metabolic only)
model = read_sbml('path/to/model.xml', regulatory=False)

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
from cobra.io import read_sbml_model
from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems import GKOProblem

# Load model (COBRApy or REFRAMED)
model = read_sbml_model('path/to/e_coli_core.xml')
BIOMASS, PRODUCT = 'Biomass_Ecoli_core', 'EX_succ_e'

# Define optimization problem: gene knockouts that make anaerobic E. coli secrete succinate
problem = GKOProblem(
    model,
    [BPCY(BIOMASS, PRODUCT), WYIELD(BIOMASS, PRODUCT)],   # objectives
    envcond={'EX_o2_e': (0, 1000)},                       # medium: no oxygen
    candidate_max_size=6,                                 # maximum number of gene knockouts
)

# Run evolutionary algorithm
ea = EA(problem, max_generations=50, mp=True)
final_pop = ea.run()

# Display results
ea.dataframe()
```

## Examples

Examples are provided as [jupyter notebooks](https://github.com/vmspereira/MEWpy/tree/master/examples) and as [python scripts](https://github.com/vmspereira/MEWpy/tree/master/examples/scripts). The notebook `14-case-study.ipynb` follows a strain design from the model to the design, checking each step against measured data. Some notebooks need the packages of `pip install "mewpy[examples]"`.

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
pip install -e ".[germ]"    # Install the dependencies of CoRegFlux, PROM and expression preprocessing
pip install -e ".[examples]" # Install the packages some example notebooks use
```

MEWpy requires Python 3.9 or later.

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
