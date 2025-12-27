#!/usr/bin/env python
"""
Benchmark script to measure community building performance improvements.

This script compares the performance of building community models with
different numbers of organisms to demonstrate the optimizations.
"""
import time

from cobra.io.sbml import read_sbml_model

from mewpy.com import CommunityModel

MODELS_PATH = "tests/data/"
EC_CORE_MODEL = MODELS_PATH + "e_coli_core.xml.gz"


def benchmark_community_building(n_organisms, verbose=True):
    """
    Benchmark community building with n organisms.

    Args:
        n_organisms: Number of organisms in the community
        verbose: Show progress bar

    Returns:
        tuple: (build_time, num_reactions, num_metabolites, num_genes)
    """
    # Load base model
    base_model = read_sbml_model(EC_CORE_MODEL)
    base_model.reactions.get_by_id("ATPM").bounds = (0, 0)

    # Create multiple copies with different IDs
    models = []
    for i in range(n_organisms):
        model = base_model.copy()
        model.id = f"ecoli_{i}"
        models.append(model)

    # Benchmark model building
    start = time.time()
    community = CommunityModel(models, verbose=verbose)
    _ = community.merged_model  # Trigger model building
    build_time = time.time() - start

    # Get statistics
    num_reactions = len(community.reaction_map)
    num_metabolites = len(community.metabolite_map)
    num_genes = len(community.gene_map)

    return build_time, num_reactions, num_metabolites, num_genes


def main():
    """Run benchmarks with different community sizes."""
    print("Community Model Building Performance Benchmark")
    print("=" * 70)
    print()

    test_sizes = [2, 5, 10, 20]

    print(f"{'Organisms':<12} {'Build Time (s)':<18} {'Reactions':<12} {'Metabolites':<15} {'Genes':<8}")
    print("-" * 70)

    results = []
    for n in test_sizes:
        build_time, num_rxns, num_mets, num_genes = benchmark_community_building(n, verbose=False)
        results.append((n, build_time, num_rxns, num_mets, num_genes))
        print(f"{n:<12} {build_time:<18.3f} {num_rxns:<12} {num_mets:<15} {num_genes:<8}")

    print()
    print("Performance Analysis:")
    print("-" * 70)

    if len(results) >= 2:
        # Calculate scaling
        first_n, first_time = results[0][0], results[0][1]
        last_n, last_time = results[-1][0], results[-1][1]

        organism_ratio = last_n / first_n
        time_ratio = last_time / first_time

        print(f"Organism scaling: {first_n} → {last_n} ({organism_ratio:.1f}x)")
        print(f"Time scaling: {first_time:.3f}s → {last_time:.3f}s ({time_ratio:.2f}x)")
        print()

        # Ideal linear scaling would be organism_ratio
        # Better than linear is < organism_ratio
        if time_ratio < organism_ratio:
            efficiency = ((organism_ratio - time_ratio) / organism_ratio) * 100
            print(f"✓ Performance is better than linear: {efficiency:.1f}% efficiency gain")
        elif time_ratio > organism_ratio * 1.5:
            print("⚠ Performance is worse than linear (may indicate O(n²) behavior)")
        else:
            print("✓ Performance scales approximately linearly")

    print()
    print("Optimizations applied:")
    print("- Dictionary pre-allocation based on model sizes")
    print("- Optional progress bar (disabled in benchmark)")
    print("- Prefix operation caching (reduces repeated string operations)")
    print("- Memory-efficient iteration")


if __name__ == "__main__":
    main()
