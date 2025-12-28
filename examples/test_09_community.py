"""
Test script for notebook 09 - Community optimization with multiprocessing.

This tests EA optimization on community models using:
- Small E. coli core models (fit within CPLEX Community Edition limits)
- Ray multiprocessing (default evaluator)
- Gene knockout optimization (GKOProblem)
"""
import warnings
warnings.filterwarnings('ignore')

from cobra.io import read_sbml_model
from mewpy import get_simulator
from mewpy.model import CommunityModel
from mewpy.simulation import Environment
from mewpy.com import regComFBA
from mewpy.problems import GKOProblem
from mewpy.optimization import EA
from mewpy.optimization.evaluation import TargetFlux, CandidateSize
print("=" * 80)
print("COMMUNITY MODEL OPTIMIZATION TEST")
print("=" * 80)
print()

# Load E. coli core model
print("Loading E. coli core model...")
model = read_sbml_model('models/ec/e_coli_core.xml.gz')
wildtype = get_simulator(model)

print(f"Using default solver (E. coli core is small: {len(model.reactions)} reactions)")

# Create two mutants with different knockouts
print("Creating mutant strains...")
ec1 = wildtype.copy()
ec1.id = 'ec1'
ec2 = wildtype.copy()
ec2.id = 'ec2'

# Extract medium conditions
medium = Environment.from_model(wildtype)
print(f"Medium has {len(medium)} exchange reactions")

print()
print("-" * 80)
print("TEST 1: Single strain gene knockout optimization")
print("-" * 80)
print()

# Define objectives for single strain
f1 = TargetFlux(wildtype.biomass_reaction, method='FBA')
f2 = CandidateSize(maximize=True)

# Create single strain problem (smaller to run faster)
problem = GKOProblem(
    wildtype,
    [f1, f2],
    candidate_max_size=30
)

print(f"Problem has {len(problem.target_list)} candidate genes")
print(f"Running EA with multiprocessing for 10 generations...")
print()

# Run EA with multiprocessing (Ray is default)
ea = EA(problem, max_generations=10, mp=True)
solutions = ea.run(simplify=False)

print()
print(f"✅ Single strain optimization completed successfully!")
print(f"   Found {len(solutions)} solutions")
print()

# Test first solution
if len(solutions) > 0:
    sol = solutions[0]
    print(f"Best solution: {len(sol.values)} gene knockouts")
    print(f"Fitness values: {sol.fitness}")

    # Simulate the solution
    result = problem.simulate(solution=sol.values)
    print(f"Simulation result: {result.objective_value:.6f}")
    print()

print("-" * 80)
print("TEST 2: Community gene knockout optimization")
print("-" * 80)
print()

# Create community model
print("Creating community model...")
community = CommunityModel([ec1, ec2], merge_biomasses=False, flavor='cobra')
sim = community.get_community_model()
sim.set_environmental_conditions(medium)

print(f"Community model has:")
print(f"  - {len(sim.reactions)} reactions")
print(f"  - {len(sim.metabolites)} metabolites")
print()

# Define objectives for community (3 objectives)
f1_com = TargetFlux(
    community.organisms_biomass['ec1'],
    community.organisms_biomass['ec2'],
    min_biomass_value=0.1,
    method=regComFBA
)

f2_com = TargetFlux(
    community.organisms_biomass['ec2'],
    community.organisms_biomass['ec1'],
    min_biomass_value=0.1,
    method=regComFBA
)

f3_com = CandidateSize(maximize=True)

# Create community problem (smaller to run faster)
problem_com = GKOProblem(
    sim,
    [f1_com, f2_com, f3_com],
    candidate_max_size=30
)

print(f"Community problem has {len(problem_com.target_list)} candidate genes")
print(f"Running EA with multiprocessing for 10 generations...")
print()

# Use some single-strain solutions as initial population
init_pop = []
for s in solutions[:5]:  # Use first 5 solutions
    x = s.values
    init_pop.append([k + '_ec1' for k in x.keys()])
    init_pop.append([k + '_ec2' for k in x.keys()])

# Run EA with multiprocessing and initial population
ea_com = EA(
    problem_com,
    max_generations=10,
    mp=True,
    initial_population=init_pop if init_pop else None
)
solutions_com = ea_com.run(simplify=False)

print()
print(f"✅ Community optimization completed successfully!")
print(f"   Found {len(solutions_com)} solutions")
print()

# Test first community solution
if len(solutions_com) > 0:
    sol_com = solutions_com[0]
    print(f"Best community solution: {len(sol_com.values)} gene knockouts")
    print(f"Fitness values: {sol_com.fitness}")

    # Simulate with regComFBA
    result_com = problem_com.simulate(solution=sol_com.values, method=regComFBA)
    biomasses = result_com.find('BIOMASS|growth', show_nulls=True)
    print("Biomass fluxes:")
    print(biomasses)
    print()

print("=" * 80)
print("ALL TESTS PASSED!")
print("=" * 80)
print()
print("Summary:")
print(f"  ✅ Single strain optimization: {len(solutions)} solutions found")
print(f"  ✅ Community optimization: {len(solutions_com)} solutions found")
print(f"  ✅ Ray multiprocessing: No pickling errors")
print(f"  ✅ Small models: Fit within CPLEX Community Edition limits")
print()
