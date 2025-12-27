#!/usr/bin/env python
"""Test multiprocessing with CPLEX vs SCIP to investigate pickling issues"""
import sys
import pickle
sys.path.insert(0, '/Users/vpereira01/Mine/MEWpy/src')

import cobra
from mewpy.simulation import set_default_solver, get_simulator


def test_solver_pickling(solver_name):
    """Test if a solver object can be pickled"""
    print(f'\n=== Testing {solver_name.upper()} Solver Pickling ===')

    # Set solver first
    set_default_solver(solver_name)

    # Load model using cobra (works for both solvers)
    model = cobra.io.load_model('textbook')
    envcond = {'EX_glc__D_e': (-10, 100000), 'EX_o2_e': (-10, 100000)}

    try:
        simul = get_simulator(model, envcond=envcond)
    except Exception as e:
        print(f'✗ Failed to create simulator: {e}')
        return False

    print(f'Created simulator with {solver_name} solver')
    print(f'Simulator type: {type(simul).__name__}')
    print(f'Solver: {simul.solver}')

    # Try to pickle the simulator
    try:
        pickled = pickle.dumps(simul)
        print(f'✓ Simulator can be pickled ({len(pickled)} bytes)')
        unpickled = pickle.loads(pickled)
        print(f'✓ Simulator can be unpickled')
        return True
    except Exception as e:
        print(f'✗ Simulator CANNOT be pickled: {type(e).__name__}')
        print(f'  Error: {str(e)[:300]}')
        if 'swig' in str(e).lower():
            print('  >>> SWIG object detected - this is the CPLEX pickling issue <<<')
        return False


if __name__ == '__main__':
    # Test both solvers
    cplex_pickles = test_solver_pickling('cplex')
    scip_pickles = test_solver_pickling('scip')

    print('\n' + '='*70)
    print('SUMMARY')
    print('='*70)
    print(f'CPLEX simulator pickling: {"✓ WORKS" if cplex_pickles else "✗ FAILS"}')
    print(f'SCIP simulator pickling:  {"✓ WORKS" if scip_pickles else "✗ FAILS"}')

    if not cplex_pickles:
        print('\n' + '-'*70)
        print('CONCLUSION: CPLEX solver objects cannot be pickled')
        print('-'*70)
        print('The CPLEX solver uses SWIG-generated Python bindings that create')
        print('SwigPyObject instances which are not picklable. This prevents')
        print('multiprocessing from working with standard Pool.map() approaches.')
        print()
        print('SOLUTIONS:')
        print('  1. Use SCIP solver: set_default_solver("scip")')
        print('  2. Use Ray evaluator: EA(problem, mp=True) with ray installed')
        print('  3. Disable multiprocessing: EA(problem, mp=False)')
        print('-'*70)
