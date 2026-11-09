
MEWpy examples
============

This folder contains examples on how MEWpy can be used to perform phenotype simulations and Computational Strain Optimization.

## Jupyter Notebooks:

- [Simulation](01-simulation.ipynb): load models with REFRAMED, COBRApy or as GERM models, and run FBA, pFBA, lMOMA, ROOM, FVA and essentiality analysis, with the same code on each backend.
- [Optimization](02-optimization.ipynb): gene over/under-expression strain design for succinate in *E. coli*, plus inspecting and filtering the solutions.
- [Kinetic](03-kinetic.ipynb): load and simulate a kinetic model, scan parameters, and optimize enzyme levels for serine production.
- [ROUProblem](04-ROUproblem.ipynb): reaction-bound modifications that favour L-tyrosine production in yeast.
- [GOUProblem](05-GOUproblem.ipynb): gene expression modifications for succinate production in *E. coli*, compared with gene knockouts.
- [GeckoKOProblem](06-GeckoKOProblem.ipynb): enzyme knockouts and over/under-expression on an enzyme-constrained yeast model (GECKO) for succinate.
- [CoFactorSwap](07-CoFactorSwap.ipynb): swapping the cofactor specificity of reactions to increase threonine production.
- [Community](08-community.ipynb): build community models, simulate them with SteadyCom and SteadyComVA, compute SMETANA scores, and mix COBRApy and REFRAMED members.
- [CrossFeeding](09-crossfeeding.ipynb): design strains of a community that depend on each other through cross-feeding.
- [E. coli & yeast, tyrosine](10-Ecoli&Yeast-Tyrosine%20dependent.ipynb): simulate an *E. coli*–yeast co-culture on xylose and search for *E. coli* knockouts that make it supply L-tyrosine to the yeast.
- [Omics](11-omics.ipynb): integrate yeast transcriptomics with the iMM904 model using eFlux, GIMME and iMAT, build context-specific models (saved as SBML), and run the same analysis on COBRApy and REFRAMED models.
- [GERM models](12-GERM_Models.ipynb): build and edit integrated metabolic-regulatory (GERM) models.
- [GERM analysis](13-GERM_Models_analysis.ipynb): RFBA, SRFBA, PROM and CoRegFlux analyses of GERM models.

## Python scripts

- [load.py](scripts/): Examples on how to load models and get a model specific phenotype simulator.
- [geneopt.py](scripts/): Examples on how to run CSO thar modify genes expression.
- [gpr_eval.py](scripts/): Examples on how to evaluate GPRs.
- [geckosimulation.py](scripts/): Examples on how to run phenotype simulations on GECKO models.
- [cobra_geckoopt.py](scripts/): Examples on how to run CSOs using the original GECKO model.
- [geckoopt.py](scripts/): Examples on how to run CSOs using the original the REFRAME based GECKO model.
- [eacomparison.py](scripts/): Compares the performance of distincts MOEAs on solving a GECKO OU problem.
- [smoment.py](scripts/): Examples on how to run CSOs using sMOMENT and AUTOPACMEN GECKO like models for E.coli.
- [kinetic.py](scripts/): Example on how to use MEWpy for strain design using a kinetic model.
- [optorf.py](scripts/): OptORF examples.
- [optram.py](scripts/): OptRAM example.
- [germ_models_analysis.py](scripts/): Examples on how to perform GERM Model analysis in MEWpy.

More information can be found in the MEWpy documentation https://mewpy.readthedocs.io.