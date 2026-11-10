# Phenotype Simulation



## Loading metabolic models

Models can be loaded using REFRAMED or COBRApy. The examples on this page use the *E. coli* iML1515 model, `examples/models/ec/iML1515.xml` in the repository, with paths relative to the repository root:

```python
# using REFRAMED
from reframed.io.sbml import load_cbmodel
model = load_cbmodel('examples/models/ec/iML1515.xml')

# using COBRApy
from cobra.io import read_sbml_model
model = read_sbml_model('examples/models/ec/iML1515.xml')
```

The outputs below are from the COBRApy model. REFRAMED keeps the SBML prefixes in its identifiers (`R_PGI`, `M_glc__D_e`); the simulator accepts identifiers with or without them.

GERM models, read with `mewpy.io.read_model` from SBML, JSON, a COBRApy model or a REFRAMED model, have a simulator too. See [MEWpy phenotype simulation using GERM models](https://mewpy.readthedocs.io/en/latest/germ.html#mewpy-phenotype-simulation-using-germ-models).

A simulator object provides a common interface to the main phenotype analysis tasks, the same for COBRApy, REFRAMED, GERM, GECKO and sMOMENT models. The *get_simulator* function returns a simulator for the provided model, or a copy when it is given a simulator.

```python
# build a phenotype simulator
from mewpy.simulation import get_simulator
simul = get_simulator(model)
```

The simulator lists the metabolites, reactions, genes and compartments of the model:

```python
simul.metabolites[:10]
```

```
['octapb_c',
 'cysi__L_e',
 'dhap_c',
 'prbatp_c',
 '10fthf_c',
 'btal_c',
 '6pgg_c',
 'co2_e',
 'akg_e',
 'gsn_e']
```

```python
simul.reactions[:10]
```

```
['CYTDK2',
 'XPPT',
 'HXPRT',
 'NDPK5',
 'SHK3Dr',
 'NDPK6',
 'NDPK8',
 'DHORTS',
 'OMPDC',
 'PYNP2r']
```

```python
simul.genes[:10]
```

```
['b2551',
 'b0870',
 'b3368',
 'b2436',
 'b0008',
 'b3500',
 'b2465',
 'b0945',
 'b4467',
 'b3126']
```

```python
simul.compartments
```

```
{'c': 'cytosol', 'e': 'extracellular space', 'p': 'periplasm'}
```



`find` searches reactions, metabolites (`find_in='m'`) or genes (`find_in='g'`) with a regular expression, and returns a data frame. The search is case-sensitive.

```python
# find reactions identified with a term containing 'BIOMASS'
simul.find('BIOMASS')
# find metabolites identified with a term containing 'o2'
simul.find('o2', find_in='m')
```



A simulator may also be built with environmental conditions, which apply to every phenotype simulation it runs. In the next example, glucose uptake is limited to 10 mmol/gDW/h and oxygen uptake to 1000 mmol/gDW/h, that is, left unconstrained.

```python
# environmental conditions
envcond = {'EX_glc__D_e': (-10.0, 100000.0),
           'EX_o2_e':(-1000,1000)}

simul = get_simulator(model,envcond=envcond)
```

The conditions only need to be set once. The simulator sets them as bounds of the model's reactions, so the model object changes too: pass a copy, `get_simulator(model.copy(), envcond=envcond)`, to keep the original. Constraints that `simulate` or `FVA` receive for a reaction of the medium are ignored; change the medium with `simul.set_reaction_bounds` or build another simulator.

An `Environment` holds a medium as exchange reaction bounds and can be passed as `envcond`:

```python
from mewpy.simulation import Environment

medium = Environment.from_model(simul)   # the bounds of the exchange reactions
medium.get_compounds(model=simul)        # the compounds that can be taken up
```

```
['pi',
 'co2',
 'fe3',
 'h',
 'mn2',
 'fe2',
 'glc__D',
 'zn2',
 'mg2',
 'ca2',
 'ni2',
 'cu2',
 'sel',
 'cobalt2',
 'h2o',
 'mobd',
 'so4',
 'nh4',
 'k',
 'na1',
 'cl',
 'o2',
 'tungs',
 'slnt']
```

`Environment.from_compounds(['glc__D', 'o2'], max_uptake=10, model=simul)` builds a medium from compound identifiers, `Environment.empty(simul)` one that blocks every uptake, and `Environment.complete(simul, max_uptake=10)` one that allows each. `medium.apply(simul, inplace=False)` returns a medium as constraints, with the uptake of every exchange reaction it does not list blocked.



## Phenotype simulation

Phenotype simulations are run with the simulator's `simulate` method.

```python
# FBA 
result = simul.simulate()
# or 
result = simul.simulate(method='FBA')
result
```

```
============================================================
Simulation Result
============================================================
Status:              OPTIMAL
Objective (maximize): 0.876997
Method:              FBA
Model:               iML1515
Fluxes:              432 non-zero / 2712 total
Constraints:         2
  Environment:       2
Shadow prices:       1877 available
============================================================
```

Flux Balance Analysis (FBA) runs when no method is given, or with `method='FBA'`. Other phenotype simulation methods are selected with one of the identifiers:

- Flux Balance Analysis: `method = 'FBA'`
- Parsimonious FBA:`method = 'pFBA'`
- Minimization of Metabolic Adjustment:`method = 'MOMA'`
- Linear MOMA: `method = 'lMOMA'`
- Regulatory on/off minimization of metabolic flux: `method = 'ROOM'`

MOMA, lMOMA and ROOM compare a mutant's fluxes with a reference, by default the simulator's pFBA fluxes; `simulate(..., reference=fluxes)` sets another one. A reference may hold only some reactions, e.g. the fluxes a ¹³C study measured: then only those reactions count, on COBRApy and on REFRAMED. MOMA is a quadratic problem: on COBRApy it needs a solver for quadratic problems, such as CPLEX or Gurobi, and REFRAMED versions without a quadratic MOMA run lMOMA instead, with a warning. GERM models do not support these methods.

```python
# pFBA
result = simul.simulate(method = 'pFBA')
result
```

```
============================================================
Simulation Result
============================================================
Status:              OPTIMAL
Objective (minimize): 768.978
Method:              pFBA
Model:               iML1515
Fluxes:              436 non-zero / 2712 total
Constraints:         2
  Environment:       2
Shadow prices:       1877 available
============================================================
```



## Reaction fluxes

Besides the objective value and the solver status, the result holds the reaction fluxes as a dictionary:

```python
result.fluxes 
```

```
OrderedDict([('CYTDK2', 0.0),
             ('XPPT', 0.0),
             ('HXPRT', 0.0),
             ('NDPK5', 0.0),
             ('SHK3Dr', 0.3339060610638284),
             ('NDPK6', 0.0),
             ('NDPK8', 0.0),
             ...])
```

The non-zero fluxes are also available as a data frame:

```python
result.dataframe
```

```
                                   Flux rate
Reaction ID                                 
SHK3Dr                          3.339061e-01
DHORTS                         -2.898162e-01
OMPDC                           2.898162e-01
G5SD                            1.936708e-01
CS                              6.886526e+00
...                                      ...
PAI2T                           3.925019e-04
PETNT181pp                      7.566488e-20
BIOMASS_Ec_iML1515_core_75p37M  8.761202e-01
ENLIPAabctex                    1.704001e-20
DM_mththf_c                    -7.457454e-21

[547 rows x 1 columns]
```

`find` lists the non-zero fluxes of the reactions that match a regular expression, sorted by value with `sort=True`:

```python
result.find('glc|BIOMASS', sort=True)
```

```
                                Flux rate
Reaction ID                              
EX_glc__D_e                     -10.00000
BIOMASS_Ec_iML1515_core_75p37M    0.87612
```

Individual reaction fluxes are read from the dictionary. For example, the flux of the *Prephenate dehydratase* reaction in the previous pFBA simulation is:

```python
result.fluxes['PPNDH']
```

```
0.16231441204188227
```



## Retrieving and setting the model objective

The simulation objective, when running FBA or pFBA phenotype simulations, is, by default, the model objective which can be seen using the simulator.

```python
simul.objective
```

```
{'BIOMASS_Ec_iML1515_core_75p37M': 1.0}
```

The simulator may also be used to change the model objective, for example, to optimize the ATP maintenance requirement (ATPM):

```python
simul.objective = 'ATPM'
# or
simul.objective = {'ATPM':1}
```

The dictionary form defines the objective as a linear expression whose coefficients are the dictionary values.



## Adding additional constraints to phenotype simulations

Simulations may include additional constraints on reaction fluxes. The previous pFBA simulation shows that the organism does not produce L-tyrosine:

```python
simul.objective = 'BIOMASS_Ec_iML1515_core_75p37M'
result.fluxes['EX_tyr__L_e']
```

```
0.0
```

Additional constraints can make the organism produce this aromatic amino acid. We may alter, for example, the bounds of the *3-dehydroquinate dehydratase* reaction, among others, starting from its initial bounds:

```python
# initial bounds
simul.get_reaction_bounds('DHQTi')
```

```
(0.0, 1000.0)
```

```python
# additional constraint
constraints = {'DHQTi' : (1, 10000),
               'PPC'   : (16.791619483918264, 10000),
               'FADRx2': (0, 0.0),
               'PDH'   : (0, 0)}

# run a pFBA simulation accounting with the new constraint
result = simul.simulate(method='pFBA',constraints=constraints)

result.fluxes['EX_tyr__L_e']
```

```
0.7365574300745229
```

We also need to verify that the organism continues to grow:

```python
result.fluxes['BIOMASS_Ec_iML1515_core_75p37M']
```

```
0.6912344174010663
```

It is also possible to plot the production envelope:

```python
from mewpy.visualization.envelope import plot_flux_envelope
plot_flux_envelope(simul,'BIOMASS_Ec_iML1515_core_75p37M','EX_tyr__L_e',constraints = constraints)
```

![](envelope.png)

The `simulate` method includes additional parameters, such as the optimization direction. For a full description please refer to the module documentation. 



## Flux Variability Analysis

The simulator also runs Flux Variability Analysis (FVA), here for L-tyrosine, returning a dictionary:

```python
simul.FVA(reactions=['EX_tyr__L_e'])
```

```
{'EX_tyr__L_e': [0.0, 0.572122699386546]}
```

or a data frame:

```python
simul.FVA(reactions=['EX_tyr__L_e'], format='df')
```

```
   Reaction ID  Minimum   Maximum
0  EX_tyr__L_e      0.0  0.572123
```



By default, the objective is kept at 90% or more of its optimum (`obj_frac=0.9`). For example, one might want to consider a fraction of 10% from optimal growth:

```python
simul.FVA(reactions=['EX_tyr__L_e'], obj_frac=0.1)
```

```
{'EX_tyr__L_e': [0.0, 5.086103322262068]}
```

FVA considers the simulator's environmental conditions. Additional constraints may be added, or changed, such as the ones used above:

```python
simul.FVA(reactions=['EX_tyr__L_e'], constraints=constraints)
```

```
{'EX_tyr__L_e': [0.0, 1.1883110320941834]}
```

## Genes and reactions essentiality

Gene and reaction essentiality tests list, respectively, the genes and reactions whose deletion leaves the objective value (growth) below `min_growth`, 0.01 by default.

```python
simul.essential_reactions()[:12]
```

```
['SHK3Dr',
 'DHORTS',
 'OMPDC',
 'G5SD',
 'CS',
 'ICDHyr',
 'PPA',
 'APRAUR',
 'DB4PS',
 'ALAR',
 'RBFK',
 'ASPTA']
```

```python
simul.essential_genes()[:12]
```

```
['b3368',
 'b1260',
 'b0109',
 'b1094',
 'b2472',
 'b1261',
 'b3770',
 'b0243',
 'b0414',
 'b3041',
 'b0025',
 'b1662']
```



## Dynamic FBA

`DynamicFBA` simulates a batch culture with the static optimisation approach (Varma & Palsson 1994; Mahadevan et al. 2002). At each time step, uptake kinetics bound the exchange reactions of the tracked compounds, a pFBA (`method=` sets another method) gives the growth rate and the exchange fluxes, and biomass and concentrations are integrated over the step. Uptake never exceeds what is left in the medium. A tracked compound without kinetics, e.g. a product, keeps the uptake bound of the model. The simulator's medium and constraints apply at every step, here an oxygen uptake of at most 15 mmol/gDW/h. Leave the tracked compounds out of the simulator's medium: a simulation ignores constraints on medium reactions, so their uptake would not follow the kinetics. The same code runs on COBRApy, REFRAMED and GERM models, and reaction ids can be given with or without the SBML prefix. The example uses the *E. coli* core model, `tests/data/e_coli_core.xml` in the repository.

```python
from cobra.io import read_sbml_model
from mewpy.simulation import get_simulator, DynamicFBA, MichaelisMenten

model = read_sbml_model('tests/data/e_coli_core.xml')
simul = get_simulator(model, envcond={'EX_o2_e': (-15, 1000)})

dfba = DynamicFBA(simul,
                  initial_biomass=0.01,                          # gDW/L
                  initial_concentrations={'EX_glc__D_e': 10,     # mM
                                          'EX_ac_e': 0},
                  uptake={'EX_glc__D_e': MichaelisMenten(vmax=10, km=0.015),   # mmol/gDW/h, mM
                          'EX_ac_e': MichaelisMenten(vmax=3, km=0.5)})
result = dfba.simulate(t_end=12, dt=0.1)
result
```

```
DFBAResult(status='stopped at t = 9 h: INFEASIBLE', time 0-9 h, steps=90, biomass 0.01 -> 0.759329 gDW/L)
```

With oxygen limiting, the cells secrete acetate while glucose lasts. Glucose runs out between 6.0 and 6.1 h, and the cells then grow on the acetate. At 9 h the remaining acetate can no longer meet the maintenance demand (ATPM), the step has no feasible solution, and the simulation stops; `result.status` reports it.

```python
result.dataframe.iloc[[0, 30, 60, 75, 90]]
```

```
    time   biomass  EX_glc__D_e   EX_ac_e
0    0.0  0.010000    10.000000  0.000000
30   3.0  0.085782     8.943846  0.720470
60   6.0  0.730284     0.000434  6.679637
75   7.5  0.748493     0.000000  3.648327
90   9.0  0.759329     0.000000  0.910551
```

`result.plot()` draws biomass and concentrations over time. An uptake rate can be any callable that takes the compound's concentration (mM) and returns its maximal uptake rate (mmol/gDW/h), e.g. to add inhibition.

### Community dynamic FBA

`CommunityDynamicFBA` simulates members that share the medium, each with its own biomass, model, uptake kinetics and constraints. In each step every member is simulated on its own. When the members together would take more of a compound than what is left, each member gets a share proportional to its demand. A member without a feasible solution in a step, e.g. one that cannot cover its maintenance, does not grow in that step, and the simulation stops when no member has one. Members interact only through the medium: what one member secretes is available to the others from the next step. The members can be a `CommunityModel` or a dict, and may use different backends; compounds are keyed by the exchange id without the SBML prefix.

Here a producer takes up glucose and secretes acetate, and a specialist without glucose transport grows on that acetate:

```python
from mewpy.simulation import CommunityDynamicFBA

o2 = {'EX_o2_e': (-15, 1000)}
producer = get_simulator(read_sbml_model('tests/data/e_coli_core.xml'), envcond=o2)
specialist = get_simulator(read_sbml_model('tests/data/e_coli_core.xml'), envcond=o2)

community = CommunityDynamicFBA({'producer': producer, 'specialist': specialist},
                                initial_biomass={'producer': 0.01, 'specialist': 0.01},
                                initial_concentrations={'EX_glc__D_e': 10, 'EX_ac_e': 0},
                                uptake={'producer': {'EX_glc__D_e': MichaelisMenten(vmax=10, km=0.015)},
                                        'specialist': {'EX_ac_e': MichaelisMenten(vmax=3, km=0.5)}},
                                constraints={'producer': {'EX_ac_e': (0, 1000)},    # no acetate uptake
                                             'specialist': {'GLCpts': 0}})          # no glucose transport
result = community.simulate(t_end=24, dt=0.1)
result
```

```
CommunityDFBAResult(status='completed', time 0-24 h, steps=240, biomass producer 0.01 -> 0.730284, specialist 0.01 -> 0.0140559 gDW/L)
```

```python
result.dataframe.iloc[[0, 34, 60, 120, 240]]
```

```
     time  producer  specialist  EX_glc__D_e   EX_ac_e
0     0.0  0.010000    0.010000    10.000000  0.000000
34    3.4  0.114245    0.010000     8.547183  0.990994
60    6.0  0.730284    0.010278     0.000434  6.614942
120  12.0  0.730284    0.011417     0.000434  6.433791
240  24.0  0.730284    0.014056     0.000434  6.010896
```

The specialist cannot cover its maintenance, and so does not grow, until acetate reaches about 1 mM at 3.4 h. The producer stops growing at 6 h: the glucose left is too little to cover its maintenance. `result.growth_rate` and `result.exchange_fluxes` hold each member's growth rate and exchange fluxes per step.



## Kinetic simulation

`KineticSimulation` integrates the ordinary differential equations of a kinetic model, an `ODEModel` read from SBML with `load_ODEModel`. The example uses the Chassagnole et al. (2002) model of the *E. coli* central carbon metabolism, `tests/data/chassagnole2002.xml` in the repository (also in `examples/models/kinetic/`):

```python
from mewpy.io.sbml import load_ODEModel
from mewpy.simulation.kinetic import KineticSimulation

kmodel = load_ODEModel('tests/data/chassagnole2002.xml')
ksim = KineticSimulation(kmodel, timeout=60)
result = ksim.simulate()
print(result.status, round(result.fluxes['vPTS'], 4), round(result.concentrations['cglcex'], 4))
```

```
ODEStatus.OPTIMAL 0.2004 0.0557
```

By default the integration runs from 0 to 1e9, which approximates the steady state; `result.fluxes` holds the reaction rates and `result.concentrations` the metabolite concentrations at the end. With a timeout, in seconds (6000 by default), the integration runs in a separate process and is stopped at the timeout; `timeout=0` runs it in the calling process.

`simulate` also takes `parameters`, parameter values that replace the model's, `factors`, which multiply parameter values (a factor of 0 on a maximal rate removes the enzyme), `initcon`, initial concentrations, and `t_points`, the integration time points:

```python
result = ksim.simulate(factors={'rmaxPGI': 0})     # no phosphoglucose isomerase
print(result.status, round(result.fluxes['vPGI'], 4), round(result.fluxes['vG6PDH'], 4))
```

```
ODEStatus.OPTIMAL 0 0.1975
```

Without phosphoglucose isomerase, the glucose 6-phosphate dehydrogenase rate rises from 0.1395 to 0.1975 mM/s.

```python
result = ksim.simulate(t_points=[0, 1, 2, 5, 10, 20])
result.get_y('cpep').round(3)       # the PEP concentration at each time point
```

```
array([2.67 , 1.813, 1.528, 0.653, 0.854, 0.864])
```

`result.plot(['cpep', 'cg6p'])` plots concentrations over time.



## Hybrid kinetic/constraint-based simulation

`HybridSimulation` constrains a constraint-based model with the rates of a kinetic model. `mapping` pairs a kinetic reaction with a reaction of the constraint-based model and a sense: 1, or -1 when their directions are opposite. A rate is converted from mM/s to mmol/gDW/h as rate × 3600 / `gDW`, the cell density in gDW per litre of cell volume (564 by default, for *E. coli*). Here the PTS rate of the Chassagnole model sets the glucose uptake of the *E. coli* core model:

```python
from cobra.io import read_sbml_model
from mewpy.simulation import get_simulator
from mewpy.simulation.hybrid import HybridSimulation

simul = get_simulator(read_sbml_model('tests/data/e_coli_core.xml'))
hybrid = HybridSimulation(kmodel, simul, mapping={'vPTS': ('EX_glc__D_e', -1)})
result = hybrid.simulate(method='pFBA')
print(round(result.fluxes['EX_glc__D_e'], 4), round(result.fluxes['Biomass_Ecoli_core'], 4))
```

```
-1.2792 0.0743
```

`simulate` runs the kinetic model, finds the fluxes that deviate least from the converted rates and from the maximal growth (a linear MOMA over the mapped reactions and the biomass reaction), fixes the mapped reactions to those fluxes, and runs `method`, pFBA by default. With `amplitude=a`, the mapped reactions are bounded to the converted rate ± a/2 instead. `nsamples` samples kinetic rates with maximal rates drawn from a log-normal distribution, and `simulate_distribution` bounds the mapped reactions between two percentiles of such samples. `HybridGeckoSimulation` does the same for GECKO models, bounding the enzyme usage reactions with the kinetic rates, maximal rates and kcat values.
