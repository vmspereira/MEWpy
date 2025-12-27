"""
Test notebook 02 with harder configuration (more likely to generate infeasible solutions).
"""
import warnings
warnings.filterwarnings('ignore')

from cobra.io import read_sbml_model
from mewpy.optimization import EA
from mewpy.optimization.evaluation import BPCY, WYIELD
from mewpy.problems import GOUProblem
from mewpy.simulation import get_simulator
from mewpy.util.constants import EAConstants

# Enable debug to see warnings
EAConstants.DEBUG = True

print("=" * 80)
print("TESTING WITH HARDER CONFIGURATION (More infeasibility expected)")
print("=" * 80)
print()

# Load model
model = read_sbml_model('models/ec/e_coli_core.xml.gz')
simul = get_simulator(model)

BIOMASS = simul.biomass_reaction
PRODUCT = 'EX_ac_e'

# Anaerobic conditions
anaerobic = {'EX_o2_e': (0, 0)}

# Evaluation functions
f1 = BPCY(BIOMASS, PRODUCT, method='pFBA')
f2 = WYIELD(BIOMASS, PRODUCT)

# LARGER candidate size (more likely to hit infeasibility)
print("Creating GOUProblem with LARGE candidate size (30 genes)...")
problem = GOUProblem(
    simul,
    [f1, f2],
    envcond=anaerobic,
    candidate_max_size=30  # Large - many modifications likely infeasible
)

print(f"Problem has {len(problem.target_list)} candidate genes")
print(f"Max modifications per solution: 30")
print()

print("Running EA for 10 generations...")
print("Expecting many infeasible solutions during search...")
print()

ea = EA(
    problem,
    max_generations=10,
    mp=True
)

solutions = ea.run(simplify=False)

print()
print("=" * 80)
print("✅ EA COMPLETED SUCCESSFULLY")
print("=" * 80)
print()
print(f"Found {len(solutions)} feasible solutions")

if len(solutions) > 0:
    # Show top 5 solutions
    print()
    print("Top 5 solutions:")
    for i, sol in enumerate(solutions[:5]):
        print(f"  {i+1}. {len(sol.values)} modifications, fitness: {sol.fitness}")

    # Simulate best
    best = solutions[0]
    result = problem.simulate(solution=best.values)
    print()
    print(f"Best solution biomass: {result.objective_value:.6f}")
else:
    print()
    print("No feasible solutions found.")
    print("This might happen if all random candidates are infeasible.")
    print("The EA continued anyway instead of crashing!")

print()
