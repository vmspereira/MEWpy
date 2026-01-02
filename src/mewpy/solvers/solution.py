# Copyright (C) 2019- Centre of Biological Engineering,
#     University of Minho, Portugal

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
##############################################################################
LP Solution interface

Adapted by Vitor Pereira from Daniel Machado's REFRAMED
https://github.com/cdanielmachado/reframed
##############################################################################
"""
import re
from enum import Enum

from mewpy.simulation import SimulationResult, Simulator, SStatus, get_simulator


class Status(Enum):
    """Enumeration of possible solution status."""

    OPTIMAL = "Optimal"
    UNKNOWN = "Unknown"
    SUBOPTIMAL = "Suboptimal"
    UNBOUNDED = "Unbounded"
    INFEASIBLE = "Infeasible"
    INF_OR_UNB = "Infeasible or Unbounded"


status_mapping = {
    Status.OPTIMAL: SStatus.OPTIMAL,
    Status.UNKNOWN: SStatus.UNKNOWN,
    Status.SUBOPTIMAL: SStatus.SUBOPTIMAL,
    Status.UNBOUNDED: SStatus.UNBOUNDED,
    Status.INFEASIBLE: SStatus.INFEASIBLE,
    Status.INF_OR_UNB: SStatus.INF_OR_UNB,
}


class Solution(object):
    """Stores the results of an optimization.

    Instantiate without arguments to create an empty Solution representing a failed optimization.
    """

    def __init__(
        self,
        status=Status.UNKNOWN,
        message=None,
        fobj=None,
        values=None,
        shadow_prices=None,
        reduced_costs=None,
        method=None,
        model=None,
        simulator=None,
        objective_direction="maximize",
        objective_value=None,
    ):
        # Handle backward compatibility: objective_value is an alias for fobj
        if objective_value is not None and fobj is None:
            fobj = objective_value

        self.status = status
        self.message = message
        self.fobj = fobj
        self.values = values or {}
        self.shadow_prices = shadow_prices or {}
        self.reduced_costs = reduced_costs or {}
        # Additional attributes for ModelSolution compatibility
        self._method = method
        self._model = model
        self._simulator = simulator
        self._objective_direction = objective_direction

    @property
    def objective_value(self):
        """Backward compatibility alias for fobj attribute (ModelSolution API)."""
        return self.fobj

    @property
    def x(self):
        """Backward compatibility alias for values attribute (ModelSolution API)."""
        return self.values

    @property
    def method(self):
        """The analysis method used to obtain the solution."""
        return self._method

    @property
    def model(self):
        """The model used to obtain the solution."""
        return self._model

    @property
    def simulator(self):
        """The simulator used to obtain the solution."""
        return self._simulator

    @property
    def objective_direction(self):
        """The direction of the objective function."""
        return self._objective_direction

    def __str__(self):
        """Simple string representation for backward compatibility."""
        return f"Objective: {self.fobj}\nStatus: {self.status.value}\n"

    def __repr__(self):
        """Rich representation matching notebook output format."""
        # Get method name
        method_name = self._method if self._method else "Solution"

        # Get status string
        status_str = self.status.value if self.status else "unknown"

        # Format objective value
        if self.fobj is not None:
            obj_str = f"{self.fobj:.10g}"  # Use general format to avoid scientific notation for small values
        else:
            obj_str = "None"

        # Build output lines
        lines = []
        lines.append(f"{method_name} Solution")
        lines.append(f" Objective value: {obj_str}")
        lines.append(f" Status: {status_str}")

        return "\n".join(lines)

    def _repr_html_(self):
        """HTML representation for Jupyter notebooks."""
        from mewpy.util.html_repr import render_html_table

        rows = []

        # Method
        if self._method:
            rows.append(("Method", self._method))

        # Status
        if self.status:
            status_str = self.status.value
            # Add color based on status
            if self.status == Status.OPTIMAL:
                status_str = f'<span style="color: green; font-weight: bold;">{status_str}</span>'
            elif self.status == Status.INFEASIBLE:
                status_str = f'<span style="color: red; font-weight: bold;">{status_str}</span>'
            elif self.status == Status.SUBOPTIMAL:
                status_str = f'<span style="color: orange; font-weight: bold;">{status_str}</span>'
            rows.append(("Status", status_str))

        # Objective value
        if self.fobj is not None:
            rows.append(("Objective", f"{self.fobj:.6g}"))

        # Objective direction
        if self._objective_direction:
            rows.append(("Direction", self._objective_direction.capitalize()))

        # Number of variables
        if self.values:
            rows.append(("Variables", str(len(self.values))))

        # Shadow prices count (if available)
        if self.shadow_prices:
            rows.append(("Shadow prices", str(len(self.shadow_prices))))

        # Reduced costs count (if available)
        if self.reduced_costs:
            rows.append(("Reduced costs", str(len(self.reduced_costs))))

        # Message (if available)
        if self.message:
            # Truncate long messages
            msg = self.message if len(self.message) <= 100 else self.message[:97] + "..."
            rows.append(("Message", msg))

        title = f"{self._method} Solution" if self._method else "Solution"
        return render_html_table(title, rows)

    def to_dataframe(self):
        """Convert solution to *pandas.DataFrame*

        Creates a DataFrame with values (fluxes) and optionally includes
        shadow prices and reduced costs if available.

        Returns:
            pandas.DataFrame: DataFrame with 'value' column and optional
                'shadow_price' and 'reduced_cost' columns
        """
        try:
            import pandas as pd
        except ImportError:
            raise RuntimeError("Pandas is not installed.")

        # Start with values
        df = pd.DataFrame(self.values.values(), columns=["value"], index=self.values.keys())

        # Add shadow prices if available
        if self.shadow_prices:
            shadow_series = pd.Series(self.shadow_prices, name="shadow_price")
            df = df.join(shadow_series, how="left")

        # Add reduced costs if available
        if self.reduced_costs:
            reduced_series = pd.Series(self.reduced_costs, name="reduced_cost")
            df = df.join(reduced_series, how="left")

        return df

    def to_series(self):
        """Convert solution values to pandas Series (ModelSolution compatibility)."""
        try:
            import pandas as pd
        except ImportError:
            raise RuntimeError("Pandas is not installed.")
        return pd.Series(self.values)

    def to_frame(self, dimensions=None):
        """Basic to_frame method for ModelSolution compatibility."""
        # For basic compatibility, just return the dataframe version
        return self.to_dataframe()

    def to_summary(self, dimensions=None):
        """
        Convert solution to a Summary object.

        This provides basic compatibility with the old ModelSolution API.
        Returns a simplified Summary with basic solution information.

        :param dimensions: Ignored for compatibility (kept for API consistency)
        :return: Summary object with dataframes of solution values
        """
        try:
            import pandas as pd
        except ImportError:
            raise RuntimeError("Pandas is not installed.")

        # Import Summary class
        from mewpy.germ.solution.summary import Summary

        # Create basic dataframe with all values
        df = self.to_dataframe()

        # Try to separate into inputs/outputs if model is available
        inputs = pd.DataFrame()
        outputs = pd.DataFrame()
        objective_df = pd.DataFrame()
        metabolic = pd.DataFrame()
        regulatory = pd.DataFrame()

        # Basic objective information
        if self.fobj is not None:
            objective_df = pd.DataFrame(
                [[self.fobj, self._objective_direction]],
                columns=["value", "direction"],
                index=[self._method if self._method else "objective"],
            )

        # If model is available, try to categorize variables
        if self._model is not None and hasattr(self._model, "reactions"):
            # Try to separate metabolic and regulatory variables
            metabolic_ids = set()
            regulatory_ids = set()
            exchange_ids = set()

            # Get reaction IDs (metabolic) and identify exchange reactions
            if hasattr(self._model, "reactions"):
                try:
                    reactions_dict = (
                        self._model.reactions
                        if hasattr(self._model.reactions, "keys")
                        else {r: r for r in self._model.reactions}
                    )
                    metabolic_ids = set(reactions_dict.keys())

                    # Identify exchange reactions (boundary reactions)
                    for rxn_id, rxn in reactions_dict.items():
                        try:
                            if hasattr(rxn, "boundary") and rxn.boundary:
                                exchange_ids.add(rxn_id)
                        except:
                            pass
                except:
                    pass

            # Get regulator/target IDs (regulatory)
            if hasattr(self._model, "regulators"):
                try:
                    regulatory_ids.update(
                        self._model.regulators.keys()
                        if hasattr(self._model.regulators, "keys")
                        else self._model.regulators
                    )
                except:
                    pass
            if hasattr(self._model, "targets"):
                try:
                    regulatory_ids.update(
                        self._model.targets.keys() if hasattr(self._model.targets, "keys") else self._model.targets
                    )
                except:
                    pass

            # Separate values into metabolic and regulatory (excluding exchange reactions from metabolic)
            metabolic_values = {k: v for k, v in self.values.items() if k in metabolic_ids and k not in exchange_ids}
            regulatory_values = {k: v for k, v in self.values.items() if k in regulatory_ids}

            # Separate exchange reactions into inputs (negative flux) and outputs (positive flux)
            input_values = {k: v for k, v in self.values.items() if k in exchange_ids and v < 0}
            output_values = {k: v for k, v in self.values.items() if k in exchange_ids and v > 0}

            if metabolic_values:
                metabolic = pd.DataFrame(metabolic_values.values(), columns=["value"], index=metabolic_values.keys())
            if regulatory_values:
                regulatory = pd.DataFrame(regulatory_values.values(), columns=["value"], index=regulatory_values.keys())
            if input_values:
                inputs = pd.DataFrame(input_values.values(), columns=["value"], index=input_values.keys())
            if output_values:
                outputs = pd.DataFrame(output_values.values(), columns=["value"], index=output_values.keys())

        # If we couldn't separate, put everything in metabolic
        if metabolic.empty and regulatory.empty and not df.empty:
            metabolic = df

        return Summary(
            inputs=inputs, outputs=outputs, objective=objective_df, df=df, metabolic=metabolic, regulatory=regulatory
        )

    @classmethod
    def from_solver(cls, method, solution, **kwargs):
        """Create a Solution from another solution object (ModelSolution compatibility)."""
        minimize = kwargs.pop("minimize", False)
        objective_direction = "minimize" if minimize else "maximize"

        return cls(
            status=getattr(solution, "status", Status.UNKNOWN),
            message=getattr(solution, "message", None),
            fobj=getattr(solution, "fobj", getattr(solution, "objective_value", 0)),
            values=getattr(solution, "values", {}),
            shadow_prices=getattr(solution, "shadow_prices", {}),
            reduced_costs=getattr(solution, "reduced_costs", {}),
            method=method,
            objective_direction=objective_direction,
            **kwargs,
        )


def to_simulation_result(model, objective_value, constraints, sim, solution, method=None):
    res = SimulationResult(
        model.model if isinstance(model, Simulator) else model,
        objective_value,
        status=status_mapping[solution.status],
        fluxes=solution.values,
        envcond=sim.environmental_conditions,
        model_constraints=sim._constraints.copy(),
        simul_constraints=constraints,
        method=method,
    )
    return res


def print_values(value_dict, pattern=None, sort=False, abstol=1e-9):

    values = [(key, value) for key, value in value_dict.items() if abs(value) > abstol]

    if pattern:
        re_expr = re.compile(pattern)
        values = [x for x in values if re_expr.search(x[0]) is not None]

    if sort:
        values.sort(key=lambda x: x[1])

    entries = (f"{r_id:<12} {val: .6g}" for (r_id, val) in values)

    print("\n".join(entries))


def print_balance(values, m_id, model, sort=False, percentage=False, abstol=1e-9):

    sim = get_simulator(model)
    inputs = sim.get_metabolite_producers(m_id)
    outputs = sim.get_metabolite_consumers(m_id)

    fwd_in = [
        (r_id, sim.get_reaction(r_id).stoichiometry[m_id] * values[r_id], "--> o")
        for r_id in inputs
        if values[r_id] > 0
    ]
    rev_in = [
        (r_id, sim.get_reaction(r_id).stoichiometry[m_id] * values[r_id], "o <--")
        for r_id in outputs
        if values[r_id] < 0
    ]
    fwd_out = [
        (r_id, sim.get_reaction(r_id).stoichiometry[m_id] * values[r_id], "o -->")
        for r_id in outputs
        if values[r_id] > 0
    ]
    rev_out = [
        (r_id, sim.get_reaction(r_id).stoichiometry[m_id] * values[r_id], "<-- o")
        for r_id in inputs
        if values[r_id] < 0
    ]

    flux_in = [x for x in fwd_in + rev_in if x[1] > abstol]
    flux_out = [x for x in fwd_out + rev_out if -x[1] > abstol]

    if sort:
        flux_in.sort(key=lambda x: x[1], reverse=True)
        flux_out.sort(key=lambda x: x[1], reverse=False)

    if percentage:
        turnover = sum([x[1] for x in flux_in])
        flux_in = [(x[0], x[1] / turnover, x[2]) for x in flux_in]
        flux_out = [(x[0], x[1] / turnover, x[2]) for x in flux_out]
        print_format = "[ {} ] {:<12} {:< 10.2%}"
    else:
        print_format = "[ {} ] {:<12} {:< 10.6g}"

    lines = (print_format.format(x[2], x[0], x[1]) for x in flux_in + flux_out)

    print("\n".join(lines))
