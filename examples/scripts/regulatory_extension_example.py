"""
Example: Using RegulatoryExtension with the New Architecture

This example demonstrates the new RegulatoryExtension architecture that
eliminates internal metabolic storage and makes regulatory networks extend
external metabolic models (COBRApy/reframed).

Key benefits:
- No metabolic data duplication
- Works with any COBRApy or reframed model
- Clean separation: metabolic (external) vs regulatory (GERM)
- All metabolic operations delegated to simulator
- Backwards compatible with legacy GERM models
"""
import os
from pathlib import Path

from mewpy.io import read_model, Engines, Reader
from mewpy.simulation import get_simulator
from mewpy.germ.models import RegulatoryExtension, from_cobra_model_with_regulation
from mewpy.germ.analysis import RFBA, SRFBA


def example_1_regulatory_extension_from_simulator():
    """
    Example 1: Create RegulatoryExtension from a simulator (no regulatory network).

    This demonstrates the basic delegation pattern where all metabolic
    operations are delegated to the external simulator.
    """
    print("\n" + "=" * 80)
    print("Example 1: RegulatoryExtension from Simulator (No Regulatory Network)")
    print("=" * 80)

    # Step 1: Load metabolic model using COBRApy
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    model_path = path.joinpath('models', 'germ', 'e_coli_core.xml')

    # Can use cobra directly
    import cobra
    cobra_model = cobra.io.read_sbml_model(str(model_path))
    print(f"\n✓ Loaded COBRApy model: {cobra_model.id}")
    print(f"  - Reactions: {len(cobra_model.reactions)}")
    print(f"  - Genes: {len(cobra_model.genes)}")

    # Step 2: Create simulator
    simulator = get_simulator(cobra_model)
    print(f"\n✓ Created simulator: {type(simulator).__name__}")

    # Step 3: Create RegulatoryExtension (without regulatory network)
    extension = RegulatoryExtension(simulator)
    print(f"\n✓ Created RegulatoryExtension: {extension}")
    print(f"  - Reactions (delegated): {len(extension.reactions)}")
    print(f"  - Genes (delegated): {len(extension.genes)}")
    print(f"  - Metabolites (delegated): {len(extension.metabolites)}")
    print(f"  - Has regulatory network: {extension.has_regulatory_network()}")

    # Step 4: Access metabolic data (delegated to simulator)
    rxn_data = extension.get_reaction('ACALD')
    print("\n✓ Access metabolic data (delegated):")
    print(f"  - Reaction: {rxn_data.get('id')}")
    print(f"  - Name: {rxn_data.get('name')}")
    print(f"  - GPR: {rxn_data.get('gpr', 'None')}")

    # Step 5: Run FBA via simulator (no regulatory network yet)
    print("\n✓ Running FBA via simulator:")
    solution = extension.simulator.simulate()
    print(f"  - Status: {solution.status}")
    if solution.objective_value is not None:
        print(f"  - Objective: {solution.objective_value:.6f}")
    else:
        print(f"  - Objective: N/A (status: {solution.status})")

    return extension


def example_2_regulatory_extension_with_regulatory_network():
    """
    Example 2: Create RegulatoryExtension with a regulatory network.

    This demonstrates the full integrated GERM model using the new architecture.
    """
    print("\n" + "=" * 80)
    print("Example 2: RegulatoryExtension with Regulatory Network")
    print("=" * 80)

    # Step 1: Load both metabolic and regulatory models
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'germ')

    metabolic_path = str(reg_path.joinpath('e_coli_core.xml'))
    regulatory_path = str(reg_path.joinpath('e_coli_core_trn.csv'))

    # Load metabolic model
    import cobra
    cobra_model = cobra.io.read_sbml_model(metabolic_path)
    simulator = get_simulator(cobra_model)

    # Load regulatory model
    from mewpy.germ.models import RegulatoryModel
    regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                               regulatory_path,
                               sep=',',
                               id_col=0,
                               rule_col=2,
                               aliases_cols=[1],
                               header=0)
    regulatory_model = read_model(regulatory_reader)

    print("\n✓ Loaded models:")
    print(f"  - Metabolic: {cobra_model.id} ({len(cobra_model.reactions)} reactions)")
    print(f"  - Regulatory: {len(regulatory_model.interactions)} interactions")

    # Step 2: Create integrated model using RegulatoryExtension
    integrated = RegulatoryExtension(simulator, regulatory_network=regulatory_model)

    print(f"\n✓ Created integrated model: {integrated}")
    print(f"  - Reactions (delegated): {len(integrated.reactions)}")
    print(f"  - Genes (delegated): {len(integrated.genes)}")
    print(f"  - Regulators (stored): {len(integrated.regulators)}")
    print(f"  - Targets (stored): {len(integrated.targets)}")
    print(f"  - Interactions (stored): {len(integrated.interactions)}")
    print(f"  - Has regulatory network: {integrated.has_regulatory_network()}")

    # Step 3: Access regulatory data
    print("\n✓ Regulatory network iteration:")
    for i, (int_id, interaction) in enumerate(integrated.yield_interactions()):
        if i < 3:  # Show first 3
            print(f"  - {interaction.target.id}: {len(interaction.regulators)} regulators")
        elif i == 3:
            print(f"  - ... and {len(integrated.interactions) - 3} more")
            break

    # Step 4: Run RFBA with regulatory constraints
    print("\n✓ Running RFBA with regulatory network:")
    rfba = RFBA(integrated)
    rfba.build()

    # Steady-state RFBA
    solution = rfba.optimize()
    print("  - Steady-state:")
    print(f"    - Status: {solution.status}")
    print(f"    - Objective: {solution.objective_value:.6f}")

    # Dynamic RFBA
    solution_dynamic = rfba.optimize(dynamic=True)
    print("  - Dynamic:")
    if hasattr(solution_dynamic, 'solutions'):
        print(f"    - Iterations: {len(solution_dynamic.solutions)}")

    # Step 5: Run SRFBA (MILP-based steady-state)
    print("\n✓ Running SRFBA with regulatory network:")
    srfba = SRFBA(integrated)
    srfba.build()
    solution = srfba.optimize()
    print(f"  - Status: {solution.status}")
    print(f"  - Objective: {solution.objective_value:.6f}")

    return integrated


def example_3_factory_functions():
    """
    Example 3: Using factory functions for convenience.

    This demonstrates the convenience factory functions for creating
    RegulatoryExtension instances.
    """
    print("\n" + "=" * 80)
    print("Example 3: Using Factory Functions")
    print("=" * 80)

    # Method 1: from_cobra_model_with_regulation
    print("\n✓ Method 1: from_cobra_model_with_regulation()")
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    model_path = path.joinpath('models', 'germ', 'e_coli_core.xml')

    import cobra
    cobra_model = cobra.io.read_sbml_model(str(model_path))

    extension = from_cobra_model_with_regulation(cobra_model)
    print(f"  - Created: {extension}")
    print(f"  - Reactions: {len(extension.reactions)}")

    # Method 2: load_integrated_model
    print("\n✓ Method 2: load_integrated_model() [Note: Requires implementation]")
    print("  - This would load both metabolic and regulatory from files")
    print("  - Example: load_integrated_model('model.xml', 'regulatory.csv', backend='cobra')")

    return extension


def example_4_delegation_vs_legacy():
    """
    Example 4: Comparing RegulatoryExtension delegation with legacy models.

    This shows the architectural difference between the new and old approaches.
    """
    print("\n" + "=" * 80)
    print("Example 4: New Architecture vs Legacy (Comparison)")
    print("=" * 80)

    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    reg_path = path.joinpath('models', 'germ')

    # New Architecture: RegulatoryExtension
    print("\n✓ New Architecture (RegulatoryExtension):")
    import cobra
    cobra_model = cobra.io.read_sbml_model(str(reg_path.joinpath('e_coli_core.xml')))
    simulator = get_simulator(cobra_model)
    extension = RegulatoryExtension(simulator)

    print(f"  - Type: {type(extension).__name__}")
    print(f"  - Metabolic data: Delegated to {type(simulator).__name__}")
    print("  - Regulatory data: Stored in RegulatoryExtension")
    print("  - Memory: No duplication (single source of truth)")

    # Legacy Architecture: read_model
    print("\n✓ Legacy Architecture (read_model):")
    metabolic_reader = Reader(Engines.MetabolicSBML, str(reg_path.joinpath('e_coli_core.xml')))
    regulatory_reader = Reader(Engines.BooleanRegulatoryCSV,
                               str(reg_path.joinpath('e_coli_core_trn.csv')),
                               sep=',',
                               id_col=0,
                               rule_col=2,
                               aliases_cols=[1],
                               header=0)
    legacy_model = read_model(metabolic_reader, regulatory_reader)

    print(f"  - Type: {type(legacy_model).__name__}")
    print("  - Metabolic data: Stored internally as GERM variables")
    print("  - Regulatory data: Stored internally")
    print("  - Memory: Some duplication")

    print("\n✓ Backwards compatibility:")
    print("  - Both work with RFBA, SRFBA, PROM, CoRegFlux")
    print("  - Legacy models still supported")
    print("  - No breaking changes")


def example_5_prom_analysis():
    """Example 5: PROM analysis with RegulatoryExtension."""
    print("\n" + "=" * 80)
    print("Example 5: PROM Analysis")
    print("=" * 80)

    from mewpy.germ.analysis import PROM

    # Load model with regulatory network
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    model_path = path.joinpath('models', 'germ', 'e_coli_core.xml')
    reg_path = path.joinpath('models', 'germ', 'e_coli_core_trn.csv')

    model = RegulatoryExtension.from_sbml(
        str(model_path),
        str(reg_path),
        regulatory_format='csv',
        sep=',',
        flavor='reframed'
    )

    print("\n✓ PROM (Probabilistic Regulation of Metabolism):")
    print(f"  - Regulators: {len(model.regulators)}")
    print(f"  - Targets: {len(model.targets)}")

    # Create PROM instance
    prom = PROM(model).build()
    print(f"  - PROM method: {prom.method}")
    print(f"  - Synchronized: {prom.synchronized}")

    # Test single regulator knockout
    regulators = list(model.regulators.keys())[:2]
    if regulators:
        print(f"\n✓ Testing knockout of {len(regulators)} regulators...")

        # Run with default probabilities
        result = prom.optimize(regulators=regulators)

        for regulator in regulators:
            sol = result.solutions[f"ko_{regulator}"]
            if sol.objective_value is not None:
                print(f"  - {regulator}: objective = {sol.objective_value:.4f}")
            else:
                print(f"  - {regulator}: status = {sol.status}")

    print("\n✓ PROM Features:")
    print("  - Probabilistic regulatory constraints")
    print("  - Predicts transcriptional perturbation effects")
    print("  - FVA-based maximum rate computation")
    print("  - Works with RegulatoryExtension API")


def example_6_coregflux_analysis():
    """Example 6: CoRegFlux analysis with RegulatoryExtension."""
    print("\n" + "=" * 80)
    print("Example 6: CoRegFlux Analysis")
    print("=" * 80)

    from mewpy.germ.analysis import CoRegFlux

    # Load model with regulatory network
    path = Path(os.path.dirname(os.path.realpath(__file__))).parent
    model_path = path.joinpath('models', 'germ', 'e_coli_core.xml')
    reg_path = path.joinpath('models', 'germ', 'e_coli_core_trn.csv')

    model = RegulatoryExtension.from_sbml(
        str(model_path),
        str(reg_path),
        regulatory_format='csv',
        sep=',',
        flavor='reframed'
    )

    print("\n✓ CoRegFlux (Co-Regulation and Flux):")
    print(f"  - Regulators: {len(model.regulators)}")
    print(f"  - Targets: {len(model.targets)}")

    # Create CoRegFlux instance
    coregflux = CoRegFlux(model).build()
    print(f"  - Synchronized: {coregflux.synchronized}")

    # Create gene state (all genes active)
    genes = list(model.targets.keys())[:10]
    initial_state = {gene: 1.0 for gene in genes}

    print(f"\n✓ Running steady-state simulation with {len(initial_state)} genes...")
    result = coregflux.optimize(initial_state=initial_state)

    print(f"  - Status: {result.status}")
    if result.objective_value is not None:
        print(f"  - Objective: {result.objective_value:.4f}")
    else:
        print(f"  - Objective: N/A (status: {result.status})")

    # Test with reduced expression
    initial_state_reduced = {gene: 0.5 for gene in genes}
    result_reduced = coregflux.optimize(initial_state=initial_state_reduced)

    print("\n✓ With 50% reduced expression:")
    if result_reduced.objective_value is not None and result.objective_value is not None:
        print(f"  - Objective: {result_reduced.objective_value:.4f}")
        if result.objective_value > 0:
            print(f"  - Growth reduction: {(1 - result_reduced.objective_value/result.objective_value)*100:.1f}%")
    else:
        print(f"  - Status: {result_reduced.status}")

    # Dynamic simulation
    print("\n✓ Dynamic simulation with 3 time steps...")
    initial_states = [
        {gene: 1.0 for gene in genes},
        {gene: 0.8 for gene in genes},
        {gene: 0.6 for gene in genes}
    ]
    time_steps = [0.1, 0.2, 0.3]

    dynamic_result = coregflux.optimize(
        initial_state=initial_states,
        time_steps=time_steps
    )

    print(f"  - Time points simulated: {len(dynamic_result.solutions)}")
    for time_key, sol in dynamic_result.solutions.items():
        if sol.objective_value is not None:
            print(f"  - {time_key}: objective = {sol.objective_value:.4f}")
        else:
            print(f"  - {time_key}: status = {sol.status}")

    print("\n✓ CoRegFlux Features:")
    print("  - Linear regression-based gene expression prediction")
    print("  - Continuous gene state constraints")
    print("  - Dynamic simulation support")
    print("  - Metabolite and biomass tracking")
    print("  - Works with RegulatoryExtension API")


def main():
    """Run all examples."""
    print("\n" + "=" * 80)
    print("REGULATORY EXTENSION - NEW ARCHITECTURE EXAMPLES")
    print("=" * 80)
    print("\nDemonstrating the new GERM architecture that eliminates")
    print("metabolic data duplication using the decorator pattern.")

    # Run examples
    example_1_regulatory_extension_from_simulator()
    example_2_regulatory_extension_with_regulatory_network()
    example_3_factory_functions()
    example_4_delegation_vs_legacy()
    example_5_prom_analysis()
    example_6_coregflux_analysis()

    print("\n" + "=" * 80)
    print("All examples completed successfully!")
    print("=" * 80)
    print("\nKey Takeaways:")
    print("✓ RegulatoryExtension delegates all metabolic operations")
    print("✓ No metabolic data duplication")
    print("✓ Works with both COBRApy and reframed")
    print("✓ Clean separation: metabolic (external) vs regulatory (GERM)")
    print("✓ Backwards compatible with legacy models")
    print("✓ PROM and CoRegFlux fully functional with RegulatoryExtension")
    print("=" * 80 + "\n")


if __name__ == '__main__':
    main()
