# Optimization Tasks

MEWpy uses Evolutionary Computation to solve Metabolic Engineering and Systems Biology problems. For example, it finds sets of metabolic modifications that favor a phenotypic objective, such as the production of a target metabolite, and it optimizes the composition of microbial communities, the medium composition, or enzyme kinetic rates. An optimization problem combines a model with one or more optimization objectives, and takes additional arguments that depend on the problem:

```python
problem = Problem(model, [f_obj_1, f_obj_2, f_obj_3], **kwargs)
```

where `Problem` is one of the classes in `mewpy.problems` and each objective is an evaluation function from `mewpy.optimization.evaluation`.

## Optimization Engines

MEWpy runs the Evolutionary Algorithms (EAs) of two packages, jMetalPy and inspyred. When both are installed, it uses jMetalPy. Single objective problems are solved with a Genetic Algorithm (GA) or Simulated Annealing (SA). Problems with two or more objectives are solved with the Non-dominated Sorting Genetic Algorithm (NSGA-II) and, with jMetalPy, also with the Strength Pareto Evolutionary Algorithm (SPEA2) or NSGA-III, which is intended for problems with three or more objectives.

The following example runs NSGA-III on a defined problem:

```python
from mewpy.optimization import EA
ea = EA(problem, max_generations=100, algorithm='NSGAIII')
final_pop = ea.run()
```

`run` returns the final population, a list of solutions, and `ea.dataframe()` returns it as a pandas DataFrame. The [Options](options.md) page describes how to choose the engine and the algorithm.

## Strain Optimization Strategies and Methods

MEWpy supports different models, some integrating omic data, such as enzymatic expression (GECKO and sMOMENT models), others integrating regulatory networks. These models let MEWpy suggest modifications in distinct layers: reaction fluxes, gene expression, enzyme availability, regulatory genes and transcription factors. Modifications follow two strategies: deletions and over/under expression of targets.

### Reactions

**Deletion** of a reaction sets its lower and upper bounds to 0.

**Under expression** of a reaction by a level between 0 and 1 limits its flux to that fraction of a reference flux value `v`. For `v > 0` the bounds become `(0, level * v)`, and for `v < 0` they become `(level * v, 0)`.

**Over expression** by a level greater than 1 forces the flux to at least that multiple of the reference flux. For `v > 0` the bounds become `(level * v, UB)`, and for `v < 0` they become `(-UB, level * v)`, where `UB` is `ModelConstants.REACTION_UPPER_BOUND` (10000).

Both keep the direction of the reference flux. A level of 1, or a reaction with no flux in the reference, adds no constraint, and a level of 0 is a deletion. When a model splits a reversible reaction into two irreversible reactions, the one with less reference flux is deleted.

**Reference flux values** are, by default, the wild type flux distribution: a pFBA simulation of the model under the problem's environmental conditions (`envcond`) and constraints. A problem may instead use a reference given by the user, or a two-step process: the deletions in a solution are applied first, and a pFBA simulation of the resulting strain gives the reference for the over/under expressions in that solution.

The examples below use the E. coli core model, `tests/data/e_coli_core.xml` in the MEWpy repository:

```python
from reframed.io.sbml import load_cbmodel
from mewpy.problems import ROUProblem
from mewpy.optimization.evaluation import BPCY, TargetFlux

model = load_cbmodel('tests/data/e_coli_core.xml')
objectives = [BPCY('R_Biomass_Ecoli_core', 'R_EX_succ_e'), TargetFlux('R_EX_succ_e')]
medium = {'R_EX_glc__D_e': (-10, 1000), 'R_EX_o2_e': (-5, 1000)}
```

To use a customized reference, pass a dictionary of reaction: flux value as the `reference` parameter, here the fluxes of an FBA simulation:

```python
from mewpy.simulation import get_simulator

reference = get_simulator(model, envcond=medium).simulate(method='FBA').fluxes
problem = ROUProblem(model, objectives, envcond=medium, reference=reference)
```

To use the two-step process, set `twostep=True`:

```python
problem = ROUProblem(model, objectives, envcond=medium, twostep=True)
```

If that pFBA simulation is neither optimal nor suboptimal, the solution uses the wild type or user reference.

### Genes

A solution that over or under expresses a set of genes is propagated to the reactions they catalyze, converting gene expression levels into constraints over the reaction fluxes (see Fig. 2). The Boolean operators AND and OR in Gene-Protein-Reaction (GPR) rules are replaced by functions, by default MIN and MAX, which give each reaction an expression level. That level is applied to the reaction's reference flux as described for reactions. `GOUProblem` takes other functions with `operators=(and_function, or_function)`.

![](mewpy-2.png)

Fig.2: *Propagation of gene expression values to constraints over reaction fluxes*.

As an example, let us suppose a reaction R has the GPR rule '(G1 or G2) and G3'.

```python
gpr = '(G1 or G2) and G3'
```

MEWpy converts the GPR rule into a parsing tree (using a Boolean syntax) that is later used for evaluation.

```python
from mewpy.util.parsing import build_tree, Boolean

tree = build_tree(gpr, Boolean)
```

Printing the tree shows the rule with explicit parentheses:

```python
print(tree)
```

```
((G1 or G2) and G3)
```

`print_node` shows the tree structure. It writes to the `mewpy.util.parsing` logger at DEBUG level:

```python
import logging
logging.basicConfig(format="%(message)s")
logging.getLogger("mewpy.util.parsing").setLevel(logging.DEBUG)

tree.print_node()
```

```
|____&
	|____|
		|____G1
		|____G2
	|____G3
```

To evaluate the **deletion** of a gene, for example G1, define the list of *ON* genes, G2 and G3, and create a `BooleanEvaluator`:

```python
from mewpy.util.parsing import BooleanEvaluator
be = BooleanEvaluator(['G2', 'G3'])
print(tree.evaluate(be.f_operand, be.f_operator))
```

```
True
```

That is, deleting G1 has no impact on the reaction.

Now let us suppose that the only *OFF* gene is G3:

```python
be = BooleanEvaluator(['G1', 'G2'])
print(tree.evaluate(be.f_operand, be.f_operator))
```

```
False
```

As expected, deleting G3 knocks out reaction R.

**Over or under expressing genes** follows a similar process, but the Boolean operators are replaced with functions. A `GeneEvaluator` takes the gene expression levels and the functions that replace the AND and OR operators. Genes that are not in the dictionary have a level of 1.

```python
from mewpy.util.parsing import GeneEvaluator
# G1 and G3 are under expressed to half and a quarter.
# G2 is over expressed to twice its expression value.
genes = {'G1': 0.5, 'G2': 2, 'G3': 0.25}

# min replaces AND
# max replaces OR
evaluator = GeneEvaluator(genes, min, max)
```

Now evaluate the GPR rule of reaction R with this *evaluator*:

```python
print(tree.evaluate(evaluator.f_operand, evaluator.f_operator))
```

```
0.25
```

That is, the reaction flux is under expressed to 0.25 of its reference value.

### Enzymes

Genome-Scale Metabolic Models, although efficient in modeling an organism's gene-protein-reaction (GPR) interactions, are oblivious to other equally important factors (e.g. enzyme kinetics and abundance, Transcription Factors, signaling, etc.) which affect cells' metabolism and transcriptional regulation. Incorporating such elements as additional constraints leads to more accurate phenotype predictions. GECKO (GSMM with enzymatic constraints using kinetic and omics data) and sMOMENT (short MetabOlic Modelling with ENzyme kineTics) models enhance GSMMs by imposing soft constraints as upper bounds on fluxes representing the overall enzyme usage. GEMs stoichiometric matrices are extended by adding new rows that represent the enzymes and new columns that represent each enzyme's usage. Kinetic information is included, in each added row, as stoichiometric coefficients in the form of enzymes' turnover (kcat) values. Also, an identity sub-matrix over the added rows and columns models each enzyme usage upper bound, respecting therefore the constraints on each flux.

The strain optimization method modifies the organism's metabolism by over and under expressing proteins (see Fig. 3), which is reflected in the bounds of the catalyzed reactions. `GeckoKOProblem` and `GeckoOUProblem`, defined as in the previous sections, deliver solutions that delete, or over or under express, proteins whose usage is limited by enzyme upper bound constraints and by the protein pool inside the cell.

![](mewpy-3.png)

Fig.3 :*Protein over- and under expression.*

Kinetic models of metabolism may also be explored to engineer mutant strains. Changes to maximum velocity parameter values can be explored as proxies for enzyme concentrations. `KineticKOProblem` and `KineticOUProblem` delete, or over or under express, the maximum velocity parameters of a kinetic model.

### Regulatory Networks

A metabolic model alone has a significant limitation in revealing condition-specific metabolic activity because gene regulation plays an important role in constraining the particular metabolism available under any given condition. Also, the complex cross-talking mechanisms between gene regulation and metabolism are not captured by a metabolic model alone. As such, methods that systematically integrate a transcriptional regulatory network and a metabolic network have been put forward, including Regulatory Flux Balance Analysis (RFBA), the Steady-state Regulatory FBA (SRFBA), Probabilistic Regulation of Metabolism (PROM), Integrated Deduced REgulation And Metabolism (IDREAM), and CoRegFlux.

Computational strain optimization has benefited from the integration of regulatory and metabolic networks, proposing more effective metabolic engineering strategies. OptORF is one of the earliest strategies which aims to identify the optimal metabolic and regulatory gene deletions as well as gene over-expression that maximize biochemical production at the maximum cellular growth under transcriptional regulatory constraints. OptORF is a bi-level optimization problem. While the inner problem maximizes growth under the given gene deletions, the regulatory states are determined by the constraints of the outer problem. One limitation of OptORF is that it is based on a mechanical (Boolean) representation of the regulatory network and therefore ignores the range of possible regulatory intensities between regulatory factors and the target genes. Also, Boolean networks can only suggest the manipulation of transcription factors by knockouts (ON to OFF) and cannot guide more quantitative adjustment of transcriptional regulation. MEWpy implements OptORF as `OptORFProblem`, which deletes genes and regulators of an integrated model and evaluates the result with RFBA. Regulatory rules are often complex and may include conditions. The next example shows how MEWpy evaluates such a rule:

```python
# boolean example with conditions.
expression = "( (Lrp AND NOT (leu_L_e_>0)) OR NOT(((GlnG AND GlnB AND GlnD) AND RpoN) \
             AND ((glu_L_e_>0) OR (arg_L_e_>0) OR (asp_L_e_>0) OR (his_L_e_>0) OR (pro_L_e_>0) )))"
t = build_tree(expression, Boolean)

true_list = ['GlnG']
# dict of variables values
values = {'leu_L_e_': 1, 'glu_L_e_': 7, "arg_L_e_": 0.5, "asp_L_e_": 2, "his_L_e_": 0, "pro_L_e_": 0}

# propositions in the list are evaluated as True and the remaining as False.
evaluator = BooleanEvaluator(true_list, values)

res = t.evaluate(evaluator.f_operand, evaluator.f_operator)
print(res)
```

```
True
```

OptRAM (Optimization of Regulatory And Metabolic Networks) is a strain design algorithm that can identify combinatorial optimization strategies including over-expression, knockdown or knockout of both metabolic genes and transcription factors. OptRAM is based on the IDREAM integrated network framework, which makes it able to deduce a regulatory network from data. The original OptRAM uses a single objective Simulated Annealing algorithm to explore the solution search space. MEWpy implements it as `OptRamProblem`, which can also be solved with the multi-objective EAs.
