## Evaluation Functions

All evaluation functions/optimization objectives are defined in `mewpy.optimization.evaluation` module. The examples use the *E. coli* core model, `tests/data/e_coli_core.xml` in the repository, with succinate as the product:

```python
biomass_reaction_id = 'Biomass_Ecoli_core'
product_reaction_id = 'EX_succ_e'
substrate_reaction_id = 'EX_glc__D_e'
```



**Biomass-Product Coupled Yield** (BPCY):

The maximization of the Biomass-Product Coupled Yield is one of the most commonly used objectives in Computational Strain Optimization. 

![](<https://latex.codecogs.com/svg.latex?BPCY=Product\times%20Growth>) 

```python
from mewpy.optimization.evaluation import BPCY
fevaluation = BPCY(biomass_reaction_id,product_reaction_id)
```

By default, BPCY uses the pFBA fluxes of each candidate. The `method` argument selects another phenotype simulation method, such as *lMOMA*.

```python
fevaluation = BPCY(biomass_reaction_id,product_reaction_id,method ='lMOMA')
```

Also, the BPCY computation may account for a carbon source or substrate consumption:

![](<https://latex.codecogs.com/svg.latex?BPCY=\frac{Product\times%20Growth}{Substrate}>) 


```python
from mewpy.optimization.evaluation import BPCY
fevaluation = BPCY(biomass_reaction_id,product_reaction_id,
                   uptake=substrate_reaction_id)
```

A candidate without substrate uptake scores the worst fitness.



**Weighted Yield** (WYIELD):

BPCY, on its own, has some limitations. Although the BPCY score of a mutated solution may be high, the flux of the target reaction may vary widely at the maximal growth. To guide the EA to more robust solutions, WYIELD scores the range of the product flux: an FVA of the product, with the growth at 99% or more of the candidate's FBA growth (`obj_frac=0.99`), weighted as

![](<https://latex.codecogs.com/svg.latex?WYIELD=\alpha\times\text{FVA}_{max}(Product)+(1-\alpha)\times\text{FVA}_{min}(Product)>)

```python
from mewpy.optimization.evaluation import WYIELD
fevaluation = WYIELD(biomass_reaction_id,product_reaction_id)
```

The trade-off parameter is by default set to 0.3. However it may be altered by adding a new parameter when instantiating the class, for example,  `alpha=0.5`.

```python
fevaluation = WYIELD(biomass_reaction_id,product_reaction_id,alpha=0.5)
```

A candidate scores the worst fitness unless its growth times `obj_frac` is above a minimum, explicitly defined, `min_biomass_value=<some_value>`, or a fraction of the wild type growth, `min_biomass_per=0.1` (the default), that is 10%. With `scale=True`, the score is divided by the growth times `obj_frac`.



**BPCY with FVA**

MEWpy also includes an objective function that combines BPCY and WYIELD, whose formulation is:


![](<https://latex.codecogs.com/svg.latex?BPCY_{FVA}=\frac{Product\times%20Growth}{Substrate}\times\left(1-\log\left|\frac{\text{FVA}_{max}-\text{FVA}_{min}}{\text{FVA}_{max}+\text{FVA}_{min}}\right|\right)>)


```python
from mewpy.optimization.evaluation import BPCY_FVA
fevaluation = BPCY_FVA(biomass_reaction_id,product_reaction_id,uptake=substrate_reaction_id)
```

As in BPCY, the substrate is optional and fluxes may be obtained using different phenotype simulation methods. The FVA keeps the growth at 99% or more of the candidate's (`obj_frac=0.99`). A candidate whose maximal product flux is 0 scores the worst fitness, and one whose minimal and maximal product fluxes have the same magnitude scores its BPCY.

This objective function is based on a proposal from "*OptRAM*: *In-silico* strain design via integrative regulatory-metabolic network modeling" where the additional factor to BPCY favors solutions with a smaller gap between the product minimum and maximum FVA.



**Target flux**

The most straightforward objective function is the flux of a target reaction, e.g. the production of a target product, by default maximized. Given the biomass reaction, it may require a minimum growth, explicitly defined (`min_biomass_value`) or a fraction of the wild type growth (`min_biomass_per`); a candidate that grows less scores the worst fitness. The phenotype simulation method is pFBA by default and is set with `method`. The same objective function minimizes a target reaction flux with `maximize=False`.

```python
from mewpy.optimization.evaluation import TargetFlux
fevaluation = TargetFlux(product_reaction_id)
fevaluation = TargetFlux(product_reaction_id, biomass=biomass_reaction_id, min_biomass_per=0.1)
```

`TargetFlux`, `CandidateSize` and `ModificationType` also serve kinetic optimization problems (`KineticKOProblem`, `KineticOUProblem`), where `TargetFlux` reads a kinetic rate.

 

**Number of modifications**

Although problems definition allows for setting a maximum number of modifications, and solutions from the final population may be automatically simplified to remove unnecessary modifications, the minimization of the number of perturbations can be set as an additional optimization objective. `maximize=True` maximizes it instead.

```python
from mewpy.optimization.evaluation import CandidateSize
fevaluation =  CandidateSize()
```



**Modification type**

This objective favors deletions, then under expressions, then over expressions. Each modification scores 5 (deletion), 2 (under expression) or 0 (over expression), values that `penalizations` may change, and a candidate scores the mean, maximized. As such, a solution that encompasses more deletions is considered better than one with many over expressions.

```python
from mewpy.optimization.evaluation import ModificationType
fevaluation =  ModificationType()
```



**Molecular weight**

Minimizes, over a set of reactions, the sum of each flux magnitude times the molecular weight of the reaction's products (of its substrates when it has no products), in g/gDW/h.

```python
from mewpy.optimization.evaluation import MolecularWeight
fevaluation =  MolecularWeight([product_reaction_id])
```



**Number of reactions above a threshold**

Counts the reactions of a list whose pFBA flux is above `threshold` (default 0.1), maximized by default.

```python
from mewpy.optimization.evaluation import CNRFA
fevaluation = CNRFA(['EX_succ_e', 'EX_ac_e', 'EX_etoh_e'], threshold=0.1)
```



**Distance to a flux distribution**

Minimizes the Euclidean distance between the candidate's pFBA fluxes and given fluxes, over the reactions of the given fluxes.

```python
from mewpy.optimization.evaluation import FluxDistance
fevaluation = FluxDistance({'EX_succ_e': 5.0, 'EX_ac_e': 0.0})
```



**Target flux with fixed constraints**

`TargetFluxWithConstraints(reaction, constraints)` runs an FBA with `constraints` on the candidate's simulator and scores the reaction's flux, minimized by default. The candidate's simulator is passed by `KcatOptProblem`, whose candidates are simulators with modified kcat values; in other problems this objective scores the worst fitness, 1000.



**Batch titer, yield and productivity**

These objectives simulate a batch culture of each candidate by dynamic FBA (`mewpy.simulation.DynamicFBA`, see [Dynamic FBA](https://mewpy.readthedocs.io/en/latest/simulation.html#dynamic-fba)), with the candidate's modifications as constraints, and score the end of the batch, as in DySScO ([Zhuang et al., 2013](https://doi.org/10.1186/1472-6750-13-8)):

- `BatchTiter`: the product concentration at the end of the batch (mM);
- `BatchYield`: the product formed per substrate consumed (mol/mol);
- `BatchProductivity`: the product formed per batch time, `t_end` (mM/h), also when growth stops earlier.

A `BatchCulture` describes the batch, and objectives that share it simulate each candidate once. The product and the substrate must be among its tracked compounds. Leave the tracked compounds out of the problem's `envcond`: a simulation ignores constraints on medium reactions, so their uptake would follow the medium instead of `uptake`. Here is a batch with DySScO's initial conditions and glucose uptake (0.01 gDW/L of biomass, 20 mM glucose, vmax 20 mmol/gDW/h and Km 1 mM, no oxygen), without its succinate inhibition of the uptake and with a 24 h batch instead of 50 h, for succinate in E. coli core:

```python
import pandas as pd
from cobra.io import read_sbml_model
from mewpy.optimization.evaluation import BatchCulture, BatchTiter, BatchYield, BatchProductivity
from mewpy.problems import RKOProblem
from mewpy.simulation import MichaelisMenten

model = read_sbml_model('tests/data/e_coli_core.xml')
culture = BatchCulture(initial_biomass=0.01,
                       initial_concentrations={'EX_glc__D_e': 20, 'EX_succ_e': 0},
                       uptake={'EX_glc__D_e': MichaelisMenten(vmax=20, km=1)},
                       t_end=24)
objectives = [BatchTiter(culture, 'EX_succ_e'),
              BatchYield(culture, 'EX_succ_e', 'EX_glc__D_e'),
              BatchProductivity(culture, 'EX_succ_e')]
problem = RKOProblem(model, objectives, envcond={'EX_o2_e': (0, 1000)}, candidate_max_size=4)

knockouts = {'wild type': {},
             'ACALD LDH_D': {'ACALD': 0, 'LDH_D': 0},
             'ACALD LDH_D PTAr': {'ACALD': 0, 'LDH_D': 0, 'PTAr': 0}}
pd.DataFrame({name: problem.evaluate_solution(ko, decode=False) for name, ko in knockouts.items()},
             index=['titer (mM)', 'yield (mol/mol)', 'productivity (mM/h)']).T.round(3)
```

```
                  titer (mM)  yield (mol/mol)  productivity (mM/h)
wild type              0.000            0.000                0.000
ACALD LDH_D           17.306            0.877                0.721
ACALD LDH_D PTAr      18.299            0.946                0.762
```

Without acetaldehyde dehydrogenase and lactate dehydrogenase, succinate is a main fermentation product. Removing phosphotransacetylase as well raises the titer and the yield, but the culture is slower: growth stops at 20.8 h instead of 12.1 h. The productivity divides the product formed by the batch time, 24 h, for both, so it ranks them as the titer does. `culture.simulate(simulator, constraints)` gives the trajectory of a candidate, and an `EA` on the problem searches for knockouts.



**Combining two or more objectives**

The objective functions above may be combined into a linear aggregated weighed sum and used in single objective optimization algorithms, such as Genetic Algorithm or Simulated Annealing. 

![](<https://latex.codecogs.com/svg.latex?f_{agg}=\sum_{i=1}^n%20w_i\times%20f_i=w_1\times%20f_1+w_2\times%20f_2+...+w_n\times%20f_n>)

Though the sum of all weights should be equal to 1, this is not imposed as weights may also be used to introduce a normalization for each function. When not provided, the aggregated function assigns a same weight to all functions *w=1/n*.

```python
from mewpy.optimization.evaluation import BPCY, WYIELD, AggregatedSum
f1 = BPCY(biomass_reaction_id,product_reaction_id,method ='lMOMA')
f2 = WYIELD(biomass_reaction_id,product_reaction_id)

fevaluation = AggregatedSum([f1,f2],tradeoffs=[0.7,0.3])
```
