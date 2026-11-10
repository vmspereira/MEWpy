# Optimization Problems

An optimization problem combines a model, a list of objectives and problem options. The objectives are described in [Evaluation Functions](evaluation.md), including the batch titer, yield and productivity objectives, which simulate each candidate by dynamic FBA.

`mewpy.problems` exports these problem classes:

- `RKOProblem` and `ROUProblem`: reaction deletions, and reaction over/under expression;
- `GKOProblem` and `GOUProblem`: gene deletions, and gene over/under expression;
- `MediumProblem`: the uptake reactions of the medium and their uptake levels;
- `CofactorSwapProblem`: reactions whose NAD(H) cofactors are swapped for NADP(H);
- `GeckoKOProblem` and `GeckoOUProblem`: protein deletions, and protein over/under expression, in GECKO and sMOMENT models;
- `ETFLGKOProblem` and `ETFLGOUProblem`: gene deletions, and gene over/under expression, in ETFL models;
- `KineticKOProblem` and `KineticOUProblem`: changes to the maximum rate parameters of kinetic models;
- `CommunityKOProblem`: organisms removed from a community;
- `OptRamProblem`: metabolic genes and transcription factors of an OptRAM regulatory network, which `load_optram` reads into an `OptRAMRegModel`;
- `OptORFProblem`: genes and regulators of an integrated metabolic-regulatory GERM model, which `load_optorf` can build.

The constraint-based problems accept these options:

- `envcond`: the environmental conditions, a dictionary of reaction bounds;
- `constraints`: constraints added to every simulation;
- `candidate_min_size` and `candidate_max_size`: the minimum and maximum number of modifications in a solution (by default 1 and 10);
- `target`: the modification targets. Without it, each problem builds a default list. `OptRamProblem` and `OptORFProblem` always build their own;
- `non_target`: targets removed from the default list.

The constraint-based over/under expression problems also accept `levels`, the expression folds (by default `[0, 1/32, 1/16, 1/8, 1/4, 1/2, 2, 4, 8, 16, 32]`, where 0 is a deletion), `reference`, the reference flux values, and `twostep`, which computes the reference flux values of each solution after applying its deletions (by default `False`). [Optimization Tasks](optimization.md) describes how levels become constraints.

## Reaction Constraint Problems

The examples in this section use the E. coli model `iJO1366SL.xml`, in `examples/models/ec/` in the repository. Its ids encode some characters as `_DASH_`, `_LPAREN_` and `_RPAREN_`. Reaction, metabolite and gene ids can be given with or without the `R_`, `M_` or `G_` prefix: `EX_tyr_DASH_L_LPAREN_e_RPAREN_` and `R_EX_tyr_DASH_L_LPAREN_e_RPAREN_` name the same reaction.

```python
# load the model
from reframed.io.sbml import load_cbmodel
model = load_cbmodel('iJO1366SL.xml', flavor='cobra')

# Define the target
PRODUCT_ID = 'R_EX_tyr_DASH_L_LPAREN_e_RPAREN_'
BIOMASS_ID = 'R_Ec_biomass_iJO1366_core_53p95M'

# the model file defines no objective
model.set_objective({BIOMASS_ID: 1})

# environmental conditions
envcond = {'R_EX_o2_LPAREN_e_RPAREN_': (-9.66, 100000.0),
           'R_EX_glc_LPAREN_e_RPAREN_': (-12.5, 100000.0)}

# Optimization objectives
from mewpy.optimization.evaluation import BPCY, WYIELD
evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method='lMOMA')
evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)
```

The default targets of a reaction deletion problem are the reactions that are not essential, exchange or transport reactions.

```python
# build a new reaction deletion problem instance
from mewpy.problems import RKOProblem
problem = RKOProblem(model,
                     fevaluation=[evaluator_1, evaluator_2],
                     envcond=envcond)
```

The default targets of a reaction over/under expression problem are all the reactions.

```python
# build a new reaction over/under expression problem instance
from mewpy.problems import ROUProblem
problem = ROUProblem(model,
                     fevaluation=[evaluator_1, evaluator_2],
                     envcond=envcond)
```

```python
# run the optimization
from mewpy.optimization import EA
ea = EA(problem, max_generations=100, visualizer=True)
final_pop = ea.run()
```

`visualizer=True` plots the objective values during the run and requires a graphical environment.

When the optimization is concluded, duplicated solutions are removed and each solution is simplified by removing the modifications that do not change its objective values. The simplification may be skipped:

```python
final_pop = ea.run(simplify=False)
```

`MediumProblem` searches the uptake reactions of the medium. A solution sets the uptake level of each reaction it contains (by default from 0 to 10, in steps of 0.1) and closes the other uptake reactions. By default, a solution may contain all the uptake reactions. `CofactorSwapProblem` adds, for each cytosolic reaction that uses NAD and NADH, a copy that uses NADP and NADPH instead. A solution closes the original reactions and opens their copies.

## Gene Constraint Problems

Optimizations of genes' expression are run by setting and running the intended problem. Gene deletion optimization problems are defined as a GKOProblem while gene over- or under expression optimization use the GOUProblem class. The default targets of a gene deletion problem are the genes that are not essential or associated with transport reactions. The examples use the model, conditions and objectives of the previous section.

```python
# build a new problem instance
from mewpy.problems import GKOProblem
problem = GKOProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], envcond=envcond)
```

```python
# build a new problem instance
from mewpy.problems import GOUProblem
problem = GOUProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], envcond=envcond)
```

## Enzymatic Constraints Problems

MEWpy enables strain optimization using Genome-scale models enhanced with enzymatic (kcat) parameters and enzyme mass constraints:

* MEWpy supports [GECKO](https://doi.org/10.15252/msb.20167411) models, from the original COBRApy based implementation, but also implemented over REFRAMED package.

* MEWpy also supports sMOMENT and GECKO like models obtained from [AutoPACMEN](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3329-9). The repository includes AutoPACMEN models of iJO1366 in `examples/models/autopacmen/`.

The optimization API is common to the one defined for the optimization of metabolic constraints (Reactions and Genes). MEWpy automatically selects the phenotype simulator for the loaded model. The targets are protein ids, such as `P00044`, and a solution constrains the protein draw reactions (`draw_prot_P00044`; the `prot_prefix` option sets another prefix). The default targets of `GeckoKOProblem` are the proteins that are not essential, and those of `GeckoOUProblem` are all the proteins.

`GeckoModel('single-pool')` loads the ecYeastGEM model with a single protein pool, which ships with MEWpy in `src/mewpy/model/data/` (`'multi-pool'` loads the model with a pool for each measured protein). In this model, the glucose uptake is the positive flux of `r_1714_REV`.

```python
# load the model
from mewpy.model.gecko import GeckoModel
model = GeckoModel('single-pool')

# Define the target
PRODUCT_ID = 'r_1913'
BIOMASS_ID = 'r_2111'

# environmental conditions
envcond = {'r_1714_REV': (0, 12.5)}

# Optimization objectives
from mewpy.optimization.evaluation import BPCY, WYIELD
evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method="lMOMA")
evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)

# build a new problem instance for enzymatic OU
from mewpy.problems import GeckoOUProblem
problem = GeckoOUProblem(model, fevaluation=[
                         evaluator_1, evaluator_2], envcond=envcond)

# run the optimization
from mewpy.optimization import EA
ea = EA(problem, max_generations=100, visualizer=True)
final_pop = ea.run()
```

ETFL models, built with the `etfl` package, are optimized with `ETFLGKOProblem` and `ETFLGOUProblem`, which take the same options as `GKOProblem` and `GOUProblem`. In `ETFLGOUProblem`, a modification of a gene with an associated enzyme changes the bounds of the gene's translation reaction, and a modification of another gene acts through the GPRs of its reactions. In `ETFLGKOProblem`, the deletion of a gene with an associated enzyme closes the gene's translation reaction. `etfl` is not a MEWpy dependency.

## Kinetic Problems

`KineticKOProblem` and `KineticOUProblem` modify the parameters of a kinetic model. A solution multiplies each parameter it contains by a factor, and a factor of 0 is a deletion. The default targets of `KineticKOProblem` are the parameters whose name contains `rmax` or `vmax`, and those of `KineticOUProblem` are the parameters whose name contains `max` (case-insensitive). The objectives must be kinetic evaluation functions, such as `TargetFlux`.

The example uses the carbon metabolism model `chassagnole2002.xml`, in `examples/models/kinetic/` in the repository, and maximizes the serine synthesis rate, `vsersynth`. `t_points` sets the integration time span, and `timeout=0` integrates without a time limit.

```python
from mewpy.io.sbml import load_ODEModel
from mewpy.problems import KineticOUProblem
from mewpy.optimization.evaluation import TargetFlux
from mewpy.optimization import EA

model = load_ODEModel('chassagnole2002.xml')
problem = KineticOUProblem(model, [TargetFlux('vsersynth')], t_points=[0, 1e9], timeout=0)
ea = EA(problem, max_generations=20)
final_pop = ea.run()
```

## Community Problems

`CommunityKOProblem` merges a list of models into a community model and searches for the organisms to remove from the community. Models with a repeated id are discarded, so each model needs its own id. Removing an organism closes its uptake reactions from the shared compartment and its biomass reaction. The exchange reactions of the community model are unbounded. The problem does not apply `envcond` to the community model, so the medium is set on the problem's simulator.

The example builds a community of three strains of the E. coli core model `e_coli_core.xml`, in `tests/data/` in the repository, on the medium defined by the exchange bounds of that model.

```python
from reframed.io.sbml import load_cbmodel
from mewpy.problems import CommunityKOProblem
from mewpy.optimization.evaluation import TargetFlux
from mewpy.simulation import Environment

models = []
for strain in ['ec1', 'ec2', 'ec3']:
    model = load_cbmodel('e_coli_core.xml')
    model.id = strain
    models.append(model)

problem = CommunityKOProblem(models, [TargetFlux('EX_ac_e')], candidate_max_size=2)
problem.simulator.set_environmental_conditions(Environment.from_model(models[0]))
print(problem.target_list)
```

```
['ec1', 'ec2', 'ec3']
```

## Regulatory Constraints

MEWpy implements computational strain design optimization with regulatory constraints. Two methods are available, [OptRAM](https://doi.org/10.1371/journal.pcbi.1006835) and [OptORF](https://doi.org/10.1186/1752-0509-4-53).

### OptRAM Example

The example uses the yeast model `yeast_7.6-optram.xml` and the regulatory network files `regnet.csv`, `mgene.csv` and `TFnames.csv`, in `examples/models/optram/` in the repository (also in `tests/data/`). The genes of the regulatory network are matched to the model's genes with or without the `G_` prefix. The targets are the network's genes and transcription factors. By default, the expression levels are those of the other over/under expression problems, with 0.001 in place of 0.

```python
from mewpy.problems.optram import OptRamProblem, load_optram

# regulatory matrix Genes x TFs
matrix_file = 'regnet.csv'
# csv file mapping genes names entries in the regulatory matrix
gene_file = 'mgene.csv'
# csv with TFs expression
tf_file = 'TFnames.csv'

BIOMASS_ID = 'r_2111'
PRODUCT_ID = 'r_1913'  # TYR
GLC = 'r_1714'

# build the regulatory network
regnet = load_optram(gene_file, tf_file, matrix_file)

# load the model
from cobra.io import read_sbml_model
model = read_sbml_model('yeast_7.6-optram.xml')

# define the optimization objectives
from mewpy.optimization.evaluation import BPCY, WYIELD

evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID, method="lMOMA")
evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)

# environmental conditions
envcond = {GLC: (-12.5, 10000)}

# instantiate the problem
problem = OptRamProblem(model, [evaluator_1, evaluator_2],
                        regnet, candidate_min_size=1, candidate_max_size=6, envcond=envcond)

from mewpy.optimization import EA
ea = EA(problem, max_generations=100, mp=True)
final_pop = ea.run()
```

### OptORF Example

`OptORFProblem` deletes genes and regulators of an integrated GERM model: its targets are the genes, regulators and regulatory targets that are not reactions or metabolites. A solution is decoded as RFBA decodes a regulatory state (`RFBA.regulatory_constraints`): the regulators settle from the initial state and the deletions, the deletions take precedence over the decoded states, and the reactions whose GPR is false are closed. The model must be a GERM model (see [GERM models and analysis](germ.md)); COBRApy and REFRAMED models raise a `ValueError`. `load_optorf(regulatory_model, metabolic_model, config)` builds one from a regulatory and a metabolic model, given as files or readers.

The example uses the E. coli model `iJR904.xml` and its regulatory network `iMC1010.csv`, in `examples/models/germ/` in the repository.

```python
from mewpy.io import Reader, Engines, read_model
from mewpy.optimization import EA, BPCY, WYIELD
from mewpy.problems import OptORFProblem

# load a GERM model. Consult the documentation (mewpy.germ) for more details
metabolic_reader = Reader(Engines.MetabolicSBML, 'iJR904.xml')
regulatory_reader = Reader(Engines.BooleanRegulatoryCSV, 'iMC1010.csv',
                           sep=',', id_col=0, rule_col=4, aliases_cols=[1, 2, 3], header=0)
model = read_model(metabolic_reader, regulatory_reader)

BIOMASS_ID = 'BiomassEcoli'
GLC = 'EX_glc_DASH_D_e'
PRODUCT_ID = 'EX_succ_e'

# OptORF can be used with an initial state of the regulatory network.
initial_state = {
    'Stringent': 0.0,
    'high-NAD': 0.0,
    'AGDC': 0.0,
}

model.objective = {BIOMASS_ID: 1}
model.get(GLC).bounds = (-18.5, 0.0)
model.get(BIOMASS_ID).lower_bound = 0.1

evaluator_1 = BPCY(BIOMASS_ID, PRODUCT_ID)
evaluator_2 = WYIELD(BIOMASS_ID, PRODUCT_ID)
problem = OptORFProblem(model, [evaluator_1, evaluator_2], initial_state=initial_state, candidate_max_size=6)

ea = EA(problem, max_generations=10, mp=True)
final_pop = ea.run()
```
