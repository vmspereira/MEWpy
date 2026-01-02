from collections import defaultdict
from typing import TYPE_CHECKING, Dict, Optional, Sequence, Tuple, Union

import pandas as pd

from mewpy.solvers import solver_prefers_fresh_instance
from mewpy.util.constants import ModelConstants

from .analysis_utils import run_method_and_decode
from .fba import _FBA
from .pfba import _pFBA

if TYPE_CHECKING:
    from mewpy.germ.models import MetabolicModel, Model, RegulatoryModel


def fva(
    model: Union["Model", "MetabolicModel", "RegulatoryModel"],
    fraction: float = 1.0,
    reactions: Sequence[str] = None,
    objective: Union[str, Dict[str, float]] = None,
    constraints: Dict[str, Tuple[float, float]] = None,
) -> pd.DataFrame:
    """
    Flux Variability Analysis (FVA) of a metabolic model.
    FVA is a method to determine the minimum and maximum fluxes for each reaction in a metabolic model.
    It can be used to identify the reactions that are limiting the growth of a cell.
    In MEWpy, FVA is performed by solving a linear problem for each reaction in the model.

    :param model: a metabolic model to be simulated
    :param fraction: the fraction of the optimal solution to be used as the upper bound for the objective function
    (default: 1.0)
    :param reactions: the reactions to be simulated (default: all reactions in the model)
    :param objective: the objective function to be used for the simulation (default: the default objective)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :return: a pandas DataFrame with the minimum and maximum fluxes for each reaction

    Performance Note:
        This function is optimized for the current solver. With SCIP, it creates fresh
        FBA instances per reaction to avoid state management overhead. With CPLEX/Gurobi,
        it reuses a single FBA instance for better performance.
    """
    if not reactions:
        reactions = model.reactions.keys()

    if not constraints:
        constraints = {}

    if objective:
        if hasattr(objective, "keys"):
            obj = next(iter(objective.keys()))
        else:
            obj = str(objective)

    else:
        obj = next(iter(model.objective)).id

    # Get optimal objective value
    _fba = _FBA(model).build()
    objective_value, _ = run_method_and_decode(method=_fba, objective=objective, constraints=constraints)
    constraints[obj] = (fraction * objective_value, ModelConstants.REACTION_UPPER_BOUND)

    # Check if we should use fresh instances (SCIP) or reuse (CPLEX/Gurobi)
    use_fresh_instance = solver_prefers_fresh_instance()

    if not use_fresh_instance:
        # CPLEX/Gurobi: Create once and reuse
        fba = _FBA(model).build()

    result = defaultdict(list)
    for rxn in reactions:
        if use_fresh_instance:
            # SCIP: Create fresh FBA instances for each min/max
            fba_min = _FBA(model).build()
            min_val, _ = run_method_and_decode(
                method=fba_min, objective={rxn: 1.0}, constraints=constraints, minimize=True
            )

            fba_max = _FBA(model).build()
            max_val, _ = run_method_and_decode(
                method=fba_max, objective={rxn: 1.0}, constraints=constraints, minimize=False
            )
        else:
            # CPLEX/Gurobi: Reuse FBA instance
            min_val, _ = run_method_and_decode(method=fba, objective={rxn: 1.0}, constraints=constraints, minimize=True)
            max_val, _ = run_method_and_decode(
                method=fba, objective={rxn: 1.0}, constraints=constraints, minimize=False
            )

        result[rxn].append(min_val)
        result[rxn].append(max_val)

    return pd.DataFrame.from_dict(data=result, orient="index", columns=["minimum", "maximum"])


def single_gene_deletion(
    model: Union["Model", "MetabolicModel", "RegulatoryModel"],
    genes: Sequence[str] = None,
    constraints: Dict[str, Tuple[float, float]] = None,
) -> pd.DataFrame:
    """
    Single gene deletion analysis of a metabolic model.
    Single gene deletion analysis is a method to determine the effect of deleting each gene in a metabolic model.
    It can be used to identify the genes that are essential for the growth of a cell.
    In MEWpy, single gene deletion analysis is performed by solving a linear problem for each gene in the model.
    A gene knockout can switch off reactions associated with the gene, only if the gene is essential for the reaction.

    :param model: a metabolic model to be simulated
    :param genes: the genes to be simulated (default: all genes in the model)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :return: a pandas DataFrame with the fluxes for each gene

    Performance Note:
        This function is optimized for the current solver. With SCIP, it creates fresh
        FBA instances per deletion to avoid state management overhead. With CPLEX/Gurobi,
        it reuses a single FBA instance for better performance.
    """
    if not constraints:
        constraints = {}

    if not genes:
        genes = model.yield_genes()
    else:
        genes = [model.genes[gene] for gene in genes if gene in model.genes]

    # Check if we should use fresh instances (SCIP) or reuse (CPLEX/Gurobi)
    use_fresh_instance = solver_prefers_fresh_instance()

    # Get wild-type result
    if use_fresh_instance:
        wt_fba = _FBA(model).build()
        wt_objective_value, wt_status = run_method_and_decode(method=wt_fba, constraints=constraints)
    else:
        # Reuse FBA instance for all deletions (CPLEX/Gurobi)
        fba = _FBA(model).build()
        wt_objective_value, wt_status = run_method_and_decode(method=fba, constraints=constraints)

    state = {gene.id: max(gene.coefficients) for gene in model.yield_genes()}

    result = {}
    for gene in genes:

        gene_coefficient = state.pop(gene.id, 0.0)
        state[gene.id] = 0.0

        gene_constraints = {}
        for reaction in gene.yield_reactions():

            if reaction.gpr.is_none:
                continue

            gpr_eval = reaction.gpr.evaluate(values=state)

            if gpr_eval:
                continue

            gene_constraints[reaction.id] = (0.0, 0.0)

        if gene_constraints:
            if use_fresh_instance:
                # SCIP: Create fresh FBA instance for each deletion
                # This avoids freeTransform() overhead and is more stable
                gene_fba = _FBA(model).build()
                solution, status = run_method_and_decode(
                    method=gene_fba, constraints={**constraints, **gene_constraints}
                )
            else:
                # CPLEX/Gurobi: Reuse FBA instance (they handle modifications efficiently)
                solution, status = run_method_and_decode(method=fba, constraints={**constraints, **gene_constraints})

            result[gene.id] = [solution, status]

        else:
            result[gene.id] = [float(wt_objective_value), str(wt_status)]

        state[gene.id] = gene_coefficient

    return pd.DataFrame.from_dict(data=result, orient="index", columns=["growth", "status"])


def single_reaction_deletion(
    model: Union["Model", "MetabolicModel", "RegulatoryModel"],
    reactions: Sequence[str] = None,
    constraints: Dict[str, Tuple[float, float]] = None,
) -> pd.DataFrame:
    """
    Single reaction deletion analysis of a metabolic model.
    Single reaction deletion analysis is a method to determine the effect of deleting each reaction
    in a metabolic model.
    It can be used to identify the reactions that are essential for the growth of a cell.
    In MEWpy, single reaction deletion analysis is performed by solving a linear problem for each reaction in the model.

    :param model: a metabolic model to be simulated
    :param reactions: the reactions to be simulated (default: all reactions in the model)
    :param constraints: additional constraints to be used for the simulation (default: None)
    :return: a pandas DataFrame with the fluxes for each reaction

    Performance Note:
        This function is optimized for the current solver. With SCIP, it creates fresh
        FBA instances per deletion to avoid state management overhead. With CPLEX/Gurobi,
        it reuses a single FBA instance for better performance.
    """
    if not reactions:
        reactions = model.reactions.keys()

    if not constraints:
        constraints = {}

    # Check if we should use fresh instances (SCIP) or reuse (CPLEX/Gurobi)
    use_fresh_instance = solver_prefers_fresh_instance()

    if not use_fresh_instance:
        # CPLEX/Gurobi: Create once and reuse
        fba = _FBA(model).build()

    result = {}
    for reaction in reactions:
        reaction_constraints = {reaction: (0.0, 0.0)}

        if use_fresh_instance:
            # SCIP: Create fresh FBA instance for each deletion
            reaction_fba = _FBA(model).build()
            solution, status = run_method_and_decode(
                method=reaction_fba, constraints={**constraints, **reaction_constraints}
            )
        else:
            # CPLEX/Gurobi: Reuse FBA instance
            solution, status = run_method_and_decode(method=fba, constraints={**constraints, **reaction_constraints})

        result[reaction] = [solution, status]

    return pd.DataFrame.from_dict(data=result, orient="index", columns=["growth", "status"])
