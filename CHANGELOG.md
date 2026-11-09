# Changelog

## 1.1.0 (unreleased)

### Results that change: re-check work done with earlier versions

"Versions" is the earliest release found to contain the defect. It was checked by looking for the defective code in the source at
the 0.1.36 and 1.0.0 release dates. "Not checked" means the code was different, so the earlier behaviour was not verified.

| Method | What was wrong | Versions | Now |
|---|---|---|---|
| GERM `SRFBA`, `slim_srfba`, and `ifva`/`isingle_*` (default `method='srfba'`) | The GPR and regulatory rows never reached the solver, so the results were plain FBA | 1.0.0 | Regulation is enforced; e.g. knocking out the only GAPD gene in *E. coli* core gives 0 growth |
| GERM `fva`, `ifva`, `slim_*(objective=...)` | The requested objective was ignored, so FVA returned the optimum as both bounds of every reaction | 1.0.0 | Ranges match COBRApy |
| GERM `RFBA` default regulatory state | Every regulator defaulted to 1 (metabolites "present" without uptake), and regulator targets weren't updated first. Growth on *E. coli* core was 0 | 1.0.0 | 0.8514 on *E. coli* core, as in the documentation |
| GERM `PROM`, `CoRegFlux`, `slim_coregflux` | Knockouts and gene states changed no bounds (id mismatches were silently skipped). `slim_coregflux` always returned 0. Infeasible steps were reported as OPTIMAL | 1.0.0 | Applied; real status reported |
| GERM legacy models: `gene.ko()` | It didn't close the reactions whose GPR becomes false | 1.0.0 | Closed for FBA, pFBA, SRFBA, PROM and CoRegFlux |
| GERM `reaction.lower_bound = 0` / `upper_bound = 0` | 0 was treated as "unset" and replaced by ±10000 | 0.1.x and 1.0.0 | 0 is set |
| `iMAT` | Inverted big-M constraints made results infeasible or forced flux patterns. Very large bounds (e.g. ±999999) leaked flux through solver tolerances | 1.0.0 (earlier formulation not checked) | Correct big-M with M capped at 10000; the objective reports the number of agreements |
| `GIMME(build_model=True)` | Raised `TypeError` | 1.0.0 | Builds a growing context-specific model |
| `eFlux(build_model=True)` | The returned model didn't carry the expression-scaled bounds | 1.0.0 (not checked earlier) | It does, and reproduces the eFlux objective |
| `mewpy.cobra.util.add_enzyme_constraints` | Every model it built was infeasible: the enzyme-usage step also turned the protein draw reactions into consumers | not checked | Feasible; same result on COBRApy and REFRAMED |
| `SteadyCom`, `SteadyComVA` | Ran the merged-biomass formulation (attribute typo), so abundances summed to 2 | 0.1.x and 1.0.0 | Abundances sum to 1; growth rates differ |
| SMETANA `mip_score`, `mp_score`, `mro_score` | MIP compared the community with itself; medium compounds were never excluded; D/L isomers were merged; environments were not applied | 1.0.0 (not checked earlier) | Corrected |
| COBRApy `lMOMA`, `MOMA`, `ROOM` | They compared against the mutant's own optimum, not the wild type or the given `reference` | 0.1.x and 1.0.0 | They use the reference |
| REFRAMED `simulate(scalefactor=...)` | Results were divided by the factor although the LP was not scaled | 1.0.0 (not checked earlier) | Correct |
| EA with a **minimised** objective | Failed or lethal candidates scored 0.0, the best value | 0.1.x and 1.0.0 | They score `EAConstants.MAX_THRESHOLD` |
| `BPCY_FVA`, `TargetFlux` | `BPCY_FVA` added a constraint to the candidate seen by later objectives. A misspelled `TargetFlux` reaction silently returned the solver objective | not checked | Fixed; the unknown reaction raises |
| GERM SBML writer | The upper-bound parameter got the lower bound's id, so re-reading a written model set ub = lb | 1.0.0 | Round trip keeps bounds |

### Behaviour changes

- **Ids:** simulators accept reaction, metabolite and gene ids with or without the SBML prefix (REFRAMED `R_`/`M_`/`G_`, COBRApy none). Results stay keyed by the backend's native ids. Regulatory networks loaded into a `RegulatoryExtension` are renamed to the simulator's ids (`model.id_map` lists the renames). Expression data, community members, model-similarity functions, problem definitions (targets, `non_target`, constraints, `envcond`, `reference`, seeds) and evaluation functions (biomass, product and target ids) are matched the same way, so the same data and code run on both backends. Problems report native ids: a REFRAMED problem defined with `"PFL"` reports `"R_PFL"`. Community models can mix COBRApy and REFRAMED members.
- **Imports:** `import mewpy` no longer configures logging (call `mewpy.configure_logging("INFO")`). It also no longer sets the multiprocessing start method or the recursion limit. The kinetic and EA worker pools use a local `fork` context where available.
- **Logic:** Boolean rules give AND precedence over OR, and And/Or/Not return 0/1 for any numeric inputs.
- **iMAT:** `iMAT(build_model=True)` requires 90% of wild-type growth by default (`growth_frac=0` restores the old behaviour).
- **Omics inputs:** `eFlux` raises on negative expression values (use linear scale). Genes with NaN expression are treated as unmeasured.
- **TargetFlux, OptORF:** `TargetFlux` raises for an unknown reaction. OptORF solutions now list the knocked-out regulators and genes.
- **Candidate size:** problems cap the candidate size at the number of targets (this used to crash or hang). `max_candidate_size`/`min_candidate_size` are accepted as deprecated aliases of `candidate_max_size`/`candidate_min_size` (they used to be ignored).
- **Targets:** target lists are sorted, so index-based candidates from older runs map to different targets.
- **Communities:** `SteadyCom` switches the community to `merge_biomasses=False` (this rebuilds the model). With `merge_biomasses=True`, the merged-biomass ids of REFRAMED members now use the community's id prefix.
- **GERM:** `SimulatorBasedMetabolicModel.add`/`remove` and its container setters raise `NotImplementedError` (they were silent no-ops). Copies and unpickled models own their variables, and undo/redo restore the previous state.
- **Signatures:** COBRApy `Simulator.FVA()` has `processes=1` (no process pool by default). `HybridGeckoSimulation(protein_prefix=None)` derives the prefix from the backend. `Environment.from_compounds(prefix=None)` asks the simulator when `model=` is given.

### Fixed

- **Installation:** `mewpy.model`, `mewpy.io.sbml`, `mewpy.problems` and kinetic simulation failed to import on fresh installs, because the undeclared `six` package is no longer installed transitively.
- **jMetal:** the default jMetal engine failed with jmetalpy ≥ 1.10 (`create_solution(rng)`), and multi-objective runs with `algorithm="GA"`/`"SA"` crashed.
- **COBRApy simulator:** it failed whenever MEWpy's default solver was SCIP, and errors while creating a simulator were reported as "no defined simulator".
- **SCIP:** changes after a solve raised, re-solving then segfaulted, empty constraints failed, and solution pools were ignored.
- **Other solvers:** optlang/GLPK aborted on infinite bounds; optlang ignored a 0 bound and the time limit; CPLEX MIP gap parameters were swapped (untested: CPLEX not available).
- **Kinetic models:** scientific notation and unary minus in rate laws; function definitions keyed by id with argument substitution; assignment rules; `initialAmount`; reactions without a kinetic law; `deriv(t, y)` ignoring `y`.
- **I/O:** the regulatory SBML round trip lost constant rules; the COBRApy/REFRAMED model engines stored string coefficients and charges and dropped compartment names; repeated species stoichiometry; JSON, copy and pickle lost compartment names; self-regulated TRN targets became two objects in the CSV readers.
- **GECKO / cobra utilities:** sMOMENT models accept the enzyme-delivery prefix with or without `R_`; `GeckoModel` accepts a DataFrame of protein properties and adds a missing protein pool; `convert_gpr_to_dnf` now converts; `convert_to_irreversible` handles forced reverse flux; `add_enzyme_constraints` uses the given kcats; pFBA crashed on REFRAMED models.
- **GERM models:** `copy()` shared containers with the original; `deepcopy()` and `reset()` emptied models; `attach=True` now resynchronises analyses after model changes; save/load round trip.
- **Optimization:** OptRAM and OptORF simplification; jMetal OU crossover duplicated targets; HybridGecko applied only the first enzyme constraint; `KcatOptProblem` modified the user's model; `nodaemon` evaluator; pools leaked on errors; seeds are honoured via `EA(..., seed=...)`.
- **Community:** `mass_flow`; exchanges and cross-feeding for COBRApy members; `balance_exchanges=False`; community models can mix COBRApy and REFRAMED members; similarity functions work across backends.
- **Utilities and plotting:** undo/redo history; `elements()` with fractional formulas; undirected metabolic graphs; envelope `flip_x`; matplotlib 3.11 compatibility in `visualization.plot`.

### Added

- `Simulator.prefix()`, `resolve()`, `resolve_many()`, `resolve_keys()`, `to_native()` and `to_canonical()` for backend-independent ids.
- `EA(..., seed=)`, `iMAT(..., biomass=, growth_frac=)`, `RegulatoryExtension.get()`, `.id_map`, `.attach()`/`.detach()`, and `decode_solver_solution(..., minimize=)`.
- An optional dependency group `mewpy[examples]` (geckopy, escher), and the notebook `examples/11-omics.ipynb` (eFlux, GIMME and iMAT on yeast transcriptomics; context-specific models; COBRApy and REFRAMED).

### Packaging

- `requires-python >= 3.9`; CI tests 3.9–3.12.
- `cplex` was removed from the `test`/`dev` extras: the PyPI build is the size-limited Community Edition, which fails on genome-scale models.
- `ray` was added to `requirements.txt`, and the bumpversion configuration was fixed.

### Removed

- Internal development notes in `docs/` and `examples/` (they were shipped in the sdist), the duplicate `README.rst`, the CPLEX/SCIP pickling probes `examples/scripts/test_mp*.py`, `tests/benchmark_community_building.py`, and the tracked `PKG-INFO`.
