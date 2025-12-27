#!/usr/bin/env python
"""
Comprehensive multiprocessing test across different combinations:
- EA engines: jmetal vs inspyred
- Simulators: reframed vs cobrapy
- Solvers: cplex vs scip
"""
import sys
import random
sys.path.insert(0, '/Users/vpereira01/Mine/MEWpy/src')

import cobra
from reframed.io.sbml import load_cbmodel
from mewpy.simulation import set_default_solver, get_simulator
from mewpy.problems import ROUProblem
from mewpy.optimization.evaluation import BPCY
from mewpy.optimization import EA, set_default_engine
from mewpy.simulation import SimulationMethod


def test_combination(ea_engine, simulator_type, solver_name):
    """Test a specific combination of EA engine, simulator, and solver"""
    print(f'\n{"="*80}')
    print(f'Testing: EA={ea_engine.upper()}, Simulator={simulator_type.upper()}, Solver={solver_name.upper()}')
    print("="*80)

    try:
        # Set engine and solver
        set_default_engine(ea_engine)
        set_default_solver(solver_name)
        print(f'✓ Set engine to {ea_engine}, solver to {solver_name}')

        # Load model based on simulator type
        if simulator_type == 'cobra':
            model = cobra.io.load_model('textbook')
            BIOMASS_ID = 'Biomass_Ecoli_core'
            PRODUCT_ID = 'EX_ac_e'
            GLC = 'EX_glc__D_e'
            O2 = 'EX_o2_e'
            print(f'✓ Loaded COBRA model: {len(model.reactions)} reactions')
        else:  # reframed
            model = load_cbmodel('../models/ec/iJO1366SL.xml', flavor='cobra')
            BIOMASS_ID = 'R_Ec_biomass_iJO1366_core_53p95M'
            PRODUCT_ID = 'R_EX_ac_LPAREN_e_RPAREN_'
            GLC = 'R_EX_glc_LPAREN_e_RPAREN_'
            O2 = 'R_EX_o2_LPAREN_e_RPAREN_'
            print(f'✓ Loaded REFRAMED model: {len(model.reactions)} reactions')

        # Create problem
        envcond = {GLC: (-10.0, 100000.0), O2: (-10.0, 100000.0)}
        evaluator = BPCY(BIOMASS_ID, PRODUCT_ID, method=SimulationMethod.lMOMA)

        problem = ROUProblem(
            model,
            fevaluation=[evaluator],
            envcond=envcond,
            candidate_max_size=3
        )
        print(f'✓ Created problem')

        # Test with multiprocessing
        print(f'Testing with mp=True (2 workers, 2 generations)...')
        ea = EA(problem, max_generations=2, mp=True, visualizer=False)

        # Try to run
        final_pop = ea.run()
        print(f'✓ SUCCESS: Completed with {len(final_pop)} solutions')
        print(f'  Sample fitness: {final_pop[0].fitness if final_pop else "N/A"}')
        return True, None

    except Exception as e:
        error_type = type(e).__name__
        error_msg = str(e)[:200]
        print(f'✗ FAILED: {error_type}')
        print(f'  Error: {error_msg}')

        # Classify error
        if 'pickle' in error_msg.lower() or 'swig' in error_msg.lower():
            print(f'  >>> PICKLING ERROR')
            return False, 'PICKLING'
        elif 'size limit' in error_msg.lower() or '1016' in error_msg:
            print(f'  >>> CPLEX COMMUNITY EDITION LIMIT')
            return False, 'SIZE_LIMIT'
        else:
            print(f'  >>> OTHER ERROR')
            return False, 'OTHER'


def run_all_tests():
    """Run tests for all combinations"""

    # Define test matrix
    ea_engines = ['jmetal', 'inspyred']
    simulator_types = ['cobra', 'reframed']
    solvers = ['cplex', 'scip']

    results = []

    for ea in ea_engines:
        for sim in simulator_types:
            for solver in solvers:
                success, error_type = test_combination(ea, sim, solver)
                results.append({
                    'ea': ea,
                    'simulator': sim,
                    'solver': solver,
                    'success': success,
                    'error': error_type
                })

    # Print summary
    print('\n\n')
    print('='*80)
    print('SUMMARY')
    print('='*80)
    print(f'{"EA":<10} {"Simulator":<12} {"Solver":<8} {"Result":<20}')
    print('-'*80)

    for r in results:
        status = '✓ SUCCESS' if r['success'] else f'✗ FAIL ({r["error"]})'
        print(f'{r["ea"]:<10} {r["simulator"]:<12} {r["solver"]:<8} {status:<20}')

    # Analysis
    print('\n' + '='*80)
    print('ANALYSIS')
    print('='*80)

    # Group by error type
    pickling_errors = [r for r in results if r['error'] == 'PICKLING']
    size_limit_errors = [r for r in results if r['error'] == 'SIZE_LIMIT']
    successes = [r for r in results if r['success']]

    if pickling_errors:
        print(f'\nPickling errors ({len(pickling_errors)}/{len(results)}):')
        for r in pickling_errors:
            print(f'  - {r["ea"]} + {r["simulator"]} + {r["solver"]}')

    if size_limit_errors:
        print(f'\nCPLEX size limit errors ({len(size_limit_errors)}/{len(results)}):')
        for r in size_limit_errors:
            print(f'  - {r["ea"]} + {r["simulator"]} + {r["solver"]}')

    if successes:
        print(f'\nSuccessful combinations ({len(successes)}/{len(results)}):')
        for r in successes:
            print(f'  - {r["ea"]} + {r["simulator"]} + {r["solver"]}')

    # Patterns
    print('\nPatterns:')

    # Check if solver matters
    cplex_fails = [r for r in results if r['solver'] == 'cplex' and not r['success']]
    scip_fails = [r for r in results if r['solver'] == 'scip' and not r['success']]
    print(f'  - CPLEX failures: {len(cplex_fails)}/{len([r for r in results if r["solver"] == "cplex"])}')
    print(f'  - SCIP failures: {len(scip_fails)}/{len([r for r in results if r["solver"] == "scip"])}')

    # Check if EA matters
    jmetal_fails = [r for r in results if r['ea'] == 'jmetal' and not r['success']]
    inspyred_fails = [r for r in results if r['ea'] == 'inspyred' and not r['success']]
    print(f'  - JMetal failures: {len(jmetal_fails)}/{len([r for r in results if r["ea"] == "jmetal"])}')
    print(f'  - Inspyred failures: {len(inspyred_fails)}/{len([r for r in results if r["ea"] == "inspyred"])}')

    # Check if simulator matters
    cobra_fails = [r for r in results if r['simulator'] == 'cobra' and not r['success']]
    reframed_fails = [r for r in results if r['simulator'] == 'reframed' and not r['success']]
    print(f'  - COBRA failures: {len(cobra_fails)}/{len([r for r in results if r["simulator"] == "cobra"])}')
    print(f'  - REFRAMED failures: {len(reframed_fails)}/{len([r for r in results if r["simulator"] == "reframed"])}')


if __name__ == '__main__':
    run_all_tests()
