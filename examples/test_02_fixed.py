"""
Test notebook 02 with infeasibility handling fix.
This should now continue instead of crashing when solutions are infeasible.
"""
import warnings
warnings.filterwarnings('ignore')

from cobra.io import read_sbml_model
from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems import GOUProblem
from mewpy.simulation import get_simulator
from mewpy.util.constants import EAConstants

# Enable debug to see infeasibility warnings
EAConstants.DEBUG = True

print("=" * 80)
print("TESTING INFEASIBILITY HANDLING IN NOTEBOOK 02")
print("=" * 80)
print()

# Load E. coli core model
print("Loading E. coli core model...")
model = read_sbml_model('models/ec/e_coli_core.xml.gz')
simul = get_simulator(model)

# Get biomass and product
BIOMASS = simul.biomass_reaction
PRODUCT = 'EX_ac_e'

print(f"Biomass reaction: {BIOMASS}")
print(f"Product reaction: {PRODUCT}")
print()

# Define anaerobic conditions (as in notebook 02)
anaerobic = {'EX_o2_e': (0, 0)}
print(f"Environment: Anaerobic {anaerobic}")
print()

# Define evaluation functions
f1 = BPCY(BIOMASS, PRODUCT, method='pFBA')
f2 = WYIELD(BIOMASS, PRODUCT)

print("Evaluation functions:")
print(f"  f1: BPCY (Biomass-Product Coupled Yield)")
print(f"  f2: WYIELD (Product yield)")
print()

# Create problem
print("Creating GOUProblem with small candidate size...")
problem = GOUProblem(
    simul,
    [f1, f2],
    envcond=anaerobic,
    candidate_max_size=10  # Small to reduce infeasibility
)

print(f"Problem has {len(problem.target_list)} candidate genes")
print()

# Run EA with multiprocessing
print("Running EA with multiprocessing for 5 generations...")
print("(This should now CONTINUE even with infeasible solutions)")
print()

ea = EA(
    problem,
    max_generations=5,
    mp=True
)

try:
    solutions = ea.run(simplify=False)

    print()
    print("=" * 80)
    print("✅ SUCCESS! EA COMPLETED WITHOUT CRASHING")
    print("=" * 80)
    print()
    print(f"Found {len(solutions)} solutions")

    if len(solutions) > 0:
        print()
        print("Best solution:")
        sol = solutions[0]
        print(f"  Genes modified: {len(sol.values)}")
        print(f"  Fitness: {sol.fitness}")

        # Try to simulate it
        result = problem.simulate(solution=sol.values)
        print(f"  Biomass flux: {result.objective_value:.6f}")
    else:
        print()
        print("⚠️  No feasible solutions found (all candidates were infeasible)")
        print("This is expected for this configuration (anaerobic + pFBA on E. coli core)")

except Exception as e:
    print()
    print("=" * 80)
    print("❌ FAILED - EA crashed with exception:")
    print("=" * 80)
    print(f"{type(e).__name__}: {e}")
    print()
    import traceback
    traceback.print_exc()

print()
print("Test complete.")
