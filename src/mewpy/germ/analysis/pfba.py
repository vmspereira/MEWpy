from typing import Dict, Union

from mewpy.germ.analysis.fba import _FBA
from mewpy.germ.models import MetabolicModel, Model, RegulatoryModel
from mewpy.solvers import solver_instance
from mewpy.solvers.solution import Solution, Status
from mewpy.solvers.solver import Solver


class _pFBA(_FBA):
    """
    Parsimonious Flux Balance Analysis (pFBA) using pure simulator-based approach.

    This implementation uses simulators as the foundation and minimizes total flux
    while maintaining optimal objective value.
    """

    def __init__(
        self,
        model: Union[Model, MetabolicModel, RegulatoryModel],
        solver: Union[str, Solver, None] = None,
        build: bool = False,
        attach: bool = False,
    ):
        """
        Parsimonious Flux Balance Analysis (pFBA) of a metabolic model.
        Pure simulator-based implementation.

        This pFBA implementation was heavily inspired by pFBA implementation of reframed python package. Take a look at
        the source: https://github.com/cdanielmachado/reframed and https://reframed.readthedocs.io/en/latest/

        For more details consult: https://doi.org/10.1038/msb.2010.47

        :param model: a MetabolicModel, RegulatoryModel or GERM model. The model is used to retrieve
        the simulator for optimization
        :param solver: A Solver, CplexSolver, GurobiSolver or OptLangSolver instance.
        Alternatively, the name of the solver is also accepted.
        The solver interface will be used to load and solve a linear problem in a given solver.
        If none, a new solver is instantiated.
        :param build: Whether to build the linear problem upon instantiation. Default: False
        :param attach: Whether to attach the linear problem to the model upon instantiation. Default: False
        """
        super().__init__(model=model, solver=solver, build=build, attach=attach)
        self.method = "pFBA"

    def build(self, fraction: float = None, constraints: Dict = None):
        """
        Build the pFBA problem using pure simulator approach.

        :param fraction: Fraction of optimal objective value to maintain. Default: None (exact optimal)
        :param constraints: Optional constraints to apply when finding optimal objective value
        """
        # Get simulator - support both RegulatoryExtension and legacy models
        if hasattr(self.model, "simulator"):
            simulator = self.model.simulator
        else:
            from mewpy.simulation import get_simulator

            try:
                simulator = get_simulator(self.model)
            except:
                simulator = self.model

        # Step 1: Create a temporary solver to find optimal objective value
        temp_solver = solver_instance(simulator)

        # Set up the biomass objective
        biomass_objective = {var.id: value for var, value in self.model.objective.items()}

        # Solve FBA to get optimal objective value (with constraints if provided)
        fba_solution = temp_solver.solve(linear=biomass_objective, minimize=False, constraints=constraints)

        if fba_solution.status != Status.OPTIMAL:
            raise RuntimeError(f"FBA failed with status: {fba_solution.status}")

        # Step 2: Create a fresh solver for pFBA optimization
        # (SCIP doesn't allow adding constraints after solving)
        self._solver = solver_instance(simulator)

        # Add constraint to maintain objective at optimal level (or fraction thereof)
        if fraction is None:
            constraint_value = fba_solution.fobj
        else:
            constraint_value = fba_solution.fobj * fraction

        # Add biomass constraint to maintain optimal growth
        self._solver.add_constraint("pfba_biomass_constraint", biomass_objective, "=", constraint_value)
        self._solver.update()

        # Step 3: Set up minimization objective (sum of absolute fluxes)
        minimize_objective = {}

        # Get all reactions from simulator
        reactions = simulator.reactions

        for r_id in reactions:
            lb, ub = simulator.get_reaction_bounds(r_id)
            if lb < 0:  # Reversible reaction - split into positive and negative parts
                pos_var = f"{r_id}_pos"
                neg_var = f"{r_id}_neg"

                # Add auxiliary variables for absolute value
                self._solver.add_variable(pos_var, 0, float("inf"), update=False)
                self._solver.add_variable(neg_var, 0, float("inf"), update=False)

                # Add constraint: r_id = pos_var - neg_var
                self._solver.add_constraint(f"split_{r_id}", {r_id: 1, pos_var: -1, neg_var: 1}, "=", 0, update=False)

                # Add to minimization objective
                minimize_objective[pos_var] = 1
                minimize_objective[neg_var] = 1
            else:  # Irreversible reaction
                minimize_objective[r_id] = 1

        self._solver.update()

        # Store the minimization objective
        self._linear_objective = minimize_objective
        self._minimize = True

        # Track the constraints used during build so we know when to rebuild
        self._build_constraints = constraints

        # Mark as synchronized
        self._synchronized = True

        # Return self for chaining
        return self

    def optimize(self, fraction: float = None, solver_kwargs: Dict = None, **kwargs) -> Solution:
        """
        Optimize the pFBA problem using pure simulator approach.

        :param fraction: Fraction of optimal objective value to maintain. Default: None (exact optimal)
        :param solver_kwargs: A dictionary of keyword arguments to be passed to the solver.
        :return: A Solution instance.
        """
        if not solver_kwargs:
            solver_kwargs = {}

        # Check if constraints are provided
        current_constraints = solver_kwargs.get("constraints") if "constraints" in solver_kwargs else None

        # Get the constraints used during the last build (if any)
        previous_constraints = getattr(self, "_build_constraints", None)

        # Need to rebuild if:
        # 1. fraction is provided
        # 2. not synchronized
        # 3. constraints changed (either added, removed, or modified)
        constraints_changed = current_constraints != previous_constraints

        if fraction is not None or not self.synchronized or constraints_changed:
            # Rebuild pFBA with the current constraints
            self.build(fraction=fraction, constraints=current_constraints)

        # Make a copy to avoid modifying the original
        solver_kwargs_copy = solver_kwargs.copy()

        # Remove conflicting arguments that we set explicitly
        solver_kwargs_copy.pop("linear", None)
        solver_kwargs_copy.pop("minimize", None)

        # Solve the parsimonious problem
        solution = self.solver.solve(linear=self._linear_objective, minimize=self._minimize, **solver_kwargs_copy)

        # Filter out auxiliary variables from solution if present
        if hasattr(solution, "values") and solution.values:
            # Keep only original reaction variables (not _pos/_neg auxiliary ones)
            filtered_values = {k: v for k, v in solution.values.items() if not ("_pos" in k or "_neg" in k)}
            # Create a new solution with filtered values
            solution.values = filtered_values

        return solution
