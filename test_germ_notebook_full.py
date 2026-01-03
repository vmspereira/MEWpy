"""Comprehensive test of GERM_Models_analysis notebook"""

import os
import warnings
from pathlib import Path

# Suppress warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from mewpy.io import Engines, Reader, read_model
from mewpy.germ.analysis import *
from mewpy.omics import ExpressionSet
from mewpy.solvers import set_default_solver

# Set SCIP as default solver
set_default_solver("scip")
print("=" * 80)
print("COMPREHENSIVE GERM NOTEBOOK TEST")
print("=" * 80)
print()

# Setup paths
path = Path("examples/models/germ")

# Define readers
core_gem_reader = Reader(Engines.MetabolicSBML, path.joinpath("e_coli_core.xml"))
core_trn_reader = Reader(
    Engines.BooleanRegulatoryCSV,
    path.joinpath("e_coli_core_trn.csv"),
    sep=",",
    id_col=0,
    rule_col=2,
    aliases_cols=[1],
    header=0,
)

imc1010_gem_reader = Reader(Engines.MetabolicSBML, path.joinpath("iJR904.xml"))
imc1010_trn_reader = Reader(
    Engines.BooleanRegulatoryCSV,
    path.joinpath("iMC1010.csv"),
    sep=",",
    id_col=0,
    rule_col=4,
    aliases_cols=[1, 2, 3],
    header=0,
)

inj661_gem_reader = Reader(Engines.MetabolicSBML, path.joinpath("iNJ661.xml"))
inj661_trn_reader = Reader(
    Engines.TargetRegulatorRegulatoryCSV,
    path.joinpath("iNJ661_trn.csv"),
    sep=";",
    target_col=0,
    regulator_col=1,
    header=None,
)
inj661_gene_expression_path = path.joinpath("iNJ661_gene_expression.csv")

imm904_gem_reader = Reader(Engines.MetabolicSBML, path.joinpath("iMM904.xml"))
imm904_trn_reader = Reader(
    Engines.CoExpressionRegulatoryCSV,
    path.joinpath("iMM904_trn.csv"),
    sep=",",
    target_col=2,
    co_activating_col=3,
    co_repressing_col=4,
    header=0,
)

# =============================================================================
# Test 1: Basic SRFBA
# =============================================================================
print("Test 1: Basic SRFBA with E. coli core")
print("-" * 80)
try:
    model = read_model(core_gem_reader, core_trn_reader)
    srfba = SRFBA(model).build()
    solution = srfba.optimize()
    print(f"✓ SRFBA initialized: {srfba}")
    print(f"  Objective value: {solution.objective_value}")
    print(f"  Status: {solution.status}")
except Exception as e:
    print(f"✗ FAILED: {e}")
print()

# =============================================================================
# Test 2: FBA and pFBA
# =============================================================================
print("Test 2: FBA and pFBA with metabolic model")
print("-" * 80)
try:
    met_model = read_model(core_gem_reader)

    # FBA
    fba_result = FBA(met_model).build().optimize()
    print(f"✓ FBA objective: {fba_result.objective_value}")

    # slim FBA
    slim_result = slim_fba(met_model)
    print(f"✓ slim_fba result: {slim_result}")

    # Simulator
    from mewpy.simulation import get_simulator, SimulationMethod
    simulator = get_simulator(met_model)
    sim_result = simulator.simulate()
    print(f"✓ Simulator FBA: {sim_result.objective_value}")

    # pFBA
    pfba_result = pFBA(met_model).build().optimize()
    print(f"✓ pFBA objective: {pfba_result.objective_value}")

    # Simulator pFBA
    sim_pfba = simulator.simulate(method=SimulationMethod.pFBA)
    print(f"✓ Simulator pFBA: {sim_pfba.objective_value}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
print()

# =============================================================================
# Test 3: FVA
# =============================================================================
print("Test 3: FVA analysis")
print("-" * 80)
try:
    met_model = read_model(core_gem_reader)
    fva_result = fva(met_model)
    print(f"✓ FVA computed: {fva_result.shape[0]} reactions")
    print(f"  Columns: {list(fva_result.columns)}")
except Exception as e:
    print(f"✗ FAILED: {e}")
print()

# =============================================================================
# Test 4: Single deletions
# =============================================================================
print("Test 4: Single reaction and gene deletions")
print("-" * 80)
try:
    met_model = read_model(core_gem_reader)

    rxn_del = single_reaction_deletion(met_model)
    print(f"✓ Reaction deletions: {rxn_del.shape[0]} reactions analyzed")

    gene_del = single_gene_deletion(met_model)
    print(f"✓ Gene deletions: {gene_del.shape[0]} genes analyzed")
except Exception as e:
    print(f"✗ FAILED: {e}")
print()

# =============================================================================
# Test 5: Regulatory truth table
# =============================================================================
print("Test 5: Regulatory truth table")
print("-" * 80)
try:
    reg_model = read_model(core_trn_reader)
    truth_table = regulatory_truth_table(reg_model)
    print(f"✓ Truth table computed: {truth_table.shape}")
except Exception as e:
    print(f"✗ FAILED: {e}")
print()

# =============================================================================
# Test 6: RFBA with iMC1010
# =============================================================================
print("Test 6: RFBA with iMC1010 model")
print("-" * 80)
try:
    model = read_model(imc1010_gem_reader, imc1010_trn_reader)
    model.objective = {"BiomassEcoli": 1}

    # Check if model is feasible with FBA first
    fba_check = FBA(model).build().optimize()
    print(f"  FBA feasibility check: {fba_check.objective_value}")

    if fba_check.objective_value > 0:
        # Try to find conflicts
        try:
            repressed_genes, repressed_reactions = find_conflicts(model)
            print(f"✓ find_conflicts completed")
            print(f"  Repressed genes: {len(repressed_genes)}")
            print(f"  Repressed reactions: {len(repressed_reactions)}")
        except Exception as e:
            print(f"  find_conflicts failed: {e}")

        # Try RFBA with empty initial state
        initial_state = {}
        rfba = RFBA(model).build()
        solution = rfba.optimize(initial_state=initial_state)
        print(f"✓ RFBA objective: {solution.objective_value}")
    else:
        print(f"  Skipping RFBA tests - model not feasible")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
print()

# =============================================================================
# Test 7: SRFBA with iMC1010
# =============================================================================
print("Test 7: SRFBA with iMC1010 model")
print("-" * 80)
try:
    model = read_model(imc1010_gem_reader, imc1010_trn_reader)
    model.objective = {"BiomassEcoli": 1}

    srfba = SRFBA(model).build()
    solution = srfba.optimize()
    print(f"✓ SRFBA objective: {solution.objective_value}")

    # slim version
    slim_result = slim_srfba(model)
    print(f"✓ slim_srfba result: {slim_result}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
print()

# =============================================================================
# Test 8: iFVA with iMC1010
# =============================================================================
print("Test 8: iFVA with iMC1010 model")
print("-" * 80)
try:
    model = read_model(imc1010_gem_reader, imc1010_trn_reader)
    model.objective = {"BiomassEcoli": 1}

    reactions_ids = list(model.reactions)[:5]
    ifva_result = ifva(
        model, fraction=0.9, reactions=reactions_ids, method="srfba"
    )
    print(f"✓ iFVA computed: {ifva_result.shape}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
print()

# =============================================================================
# Test 9: PROM with iNJ661
# =============================================================================
print("Test 9: PROM with iNJ661 model")
print("-" * 80)
try:
    model = read_model(inj661_gem_reader, inj661_trn_reader)
    model.objective = {"biomass_Mtb_9_60atp_test_NOF": 1}

    # Load expression data
    expression = ExpressionSet.from_csv(
        file_path=inj661_gene_expression_path, sep=";", index_col=0, header=None
    )
    quantile_expression, binary_expression = expression.quantile_pipeline()
    print(f"  Expression data loaded: {quantile_expression.shape}")

    # Compute probabilities
    initial_state, missed = target_regulator_interaction_probability(
        model, expression=quantile_expression, binary_expression=binary_expression
    )
    print(f"✓ PROM probabilities computed: {len(initial_state)} interactions")
    print(f"  Missed interactions: {len(missed)}")

    # Run PROM
    prom = PROM(model).build()
    solution = prom.optimize(initial_state=initial_state)
    print(f"✓ PROM optimization completed")
    print(f"  Number of KO solutions: {len(solution.solutions)}")

    # Test slim version
    slim_result = slim_prom(model, initial_state=initial_state, regulator="Rv0001")
    print(f"✓ slim_prom result: {slim_result}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
print()

# =============================================================================
# Test 10: CoRegFlux with iMM904
# =============================================================================
print("Test 10: CoRegFlux with iMM904 model")
print("-" * 80)
try:
    model = read_model(imm904_gem_reader, imm904_trn_reader)
    model.objective = {"BIOMASS_SC5_notrace": 1}

    # Load data
    expression = ExpressionSet.from_csv(
        path.joinpath("iMM904_gene_expression.csv"), sep=";", index_col=0, header=0
    ).dataframe
    influence = ExpressionSet.from_csv(
        path.joinpath("iMM904_influence.csv"), sep=";", index_col=0, header=0
    ).dataframe
    experiments = ExpressionSet.from_csv(
        path.joinpath("iMM904_experiments.csv"), sep=";", index_col=0, header=0
    ).dataframe
    print(f"  Expression: {expression.shape}")
    print(f"  Influence: {influence.shape}")
    print(f"  Experiments: {experiments.shape}")

    # Predict gene expression
    gene_expression_prediction = predict_gene_expression(
        model=model, influence=influence, expression=expression, experiments=experiments
    )
    print(f"✓ Gene expression prediction: {gene_expression_prediction.shape}")

    # Run CoRegFlux
    initial_state = list(gene_expression_prediction.to_dict().values())
    co_reg_flux = CoRegFlux(model).build()
    solution = co_reg_flux.optimize(initial_state=initial_state[0])
    print(f"✓ CoRegFlux steady-state: {solution.objective_value}")

    # Test slim version
    slim_result = slim_coregflux(model, initial_state=initial_state[0])
    print(f"✓ slim_coregflux result: {slim_result}")

    # Test dynamic simulation
    metabolites = {"glc__D_e": 16.6, "etoh_e": 0}
    biomass = 0.45
    time_steps = list(range(1, 5))  # Fewer time steps for faster test

    solution = co_reg_flux.optimize(
        initial_state=initial_state,
        metabolites=metabolites,
        biomass=biomass,
        time_steps=time_steps,
    )
    print(f"✓ CoRegFlux dynamic: {len(solution.solutions)} time steps")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
print()

# =============================================================================
# Test 11: Solution conversion methods
# =============================================================================
print("Test 11: Solution conversion methods")
print("-" * 80)
try:
    model = read_model(core_gem_reader, core_trn_reader)
    srfba = SRFBA(model).build()
    solution = srfba.optimize()

    # to_frame
    df = solution.to_frame()
    print(f"✓ to_frame(): {df.shape}")

    # to_summary
    summary = solution.to_summary()
    print(f"✓ to_summary():")
    print(f"  metabolic: {summary.metabolic.shape}")
    print(f"  regulatory: {summary.regulatory.shape}")
    print(f"  inputs: {summary.inputs.shape}")
    print(f"  outputs: {summary.outputs.shape}")
    print(f"  objective: {summary.objective}")
except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
print()

# =============================================================================
# Test 12: Attached workflow
# =============================================================================
print("Test 12: Attached workflow with model modifications")
print("-" * 80)
try:
    model = read_model(core_gem_reader, core_trn_reader)
    srfba = SRFBA(model, attach=True).build()

    wt_solution = srfba.optimize()
    print(f"  Wild-type objective: {wt_solution.objective_value}")

    # Apply knockout
    model.regulators["b3261"].ko()
    ko_solution = srfba.optimize()
    print(f"  KO objective: {ko_solution.objective_value}")

    # Restore
    model.undo()
    restored_solution = srfba.optimize()
    print(f"✓ Restored objective: {restored_solution.objective_value}")
except Exception as e:
    print(f"✗ FAILED: {e}")
print()

print("=" * 80)
print("COMPREHENSIVE TEST COMPLETE")
print("=" * 80)
