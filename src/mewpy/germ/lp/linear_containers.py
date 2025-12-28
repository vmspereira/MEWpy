from typing import Iterable, List, Union

from mewpy.solvers.solver import VarType

from .linear_utils import Node


class VariableContainer:

    def __init__(
        self,
        name: str,
        sub_variables: List[str],
        lbs: List[Union[int, float]],
        ubs: List[Union[int, float]],
        variables_type: List[VarType],
    ):
        """

        Internal use only

        A container for variables, since multiple linear variables might have to be created out of a single variable
        during problem building. It is the main object for variable management in the linear problem object

        :param name: the name of variable. It is used as key/identifier of the variable in the linear problem
        :param sub_variables: the name of all sub-variables.
        If there is a single variable that does not have sub-variables, this attribute should be filled
        with the name of the variable only
        :param lbs: a list of the lower bounds of all sub-variables
        :param ubs: a list of the upper bounds of all sub-variables
        :param variables_type: a list of the variable types of all sub-variables
        """

        if not name:
            name = None

        if not sub_variables:
            sub_variables = []

        if not lbs:
            lbs = []

        if not ubs:
            ubs = []

        if not variables_type:
            variables_type = []

        self.name = name
        self.sub_variables = sub_variables
        self.lbs = lbs
        self.ubs = ubs
        self.variables_type = variables_type

    def __len__(self):
        return len(self.sub_variables)

    def __str__(self):
        return f"Variable {self.name}"

    def __repr__(self):
        """Rich representation showing variable container details."""
        lines = []
        lines.append("=" * 60)
        lines.append(f"VariableContainer: {self.name}")
        lines.append("=" * 60)

        # Sub-variables count
        try:
            sub_var_count = len(self.sub_variables)
            lines.append(f"{'Sub-variables:':<20} {sub_var_count}")

            # Show first few sub-variables
            if sub_var_count > 0 and sub_var_count <= 3:
                for sv in self.sub_variables:
                    lines.append(f"{'  -':<20} {sv}")
            elif sub_var_count > 3:
                for sv in self.sub_variables[:3]:
                    lines.append(f"{'  -':<20} {sv}")
                lines.append(f"{'  ...':<20} and {sub_var_count - 3} more")
        except:
            pass

        # Bounds summary
        try:
            if len(self.lbs) > 0:
                lb_min = min(self.lbs)
                lb_max = max(self.lbs)
                ub_min = min(self.ubs)
                ub_max = max(self.ubs)

                if lb_min == lb_max and ub_min == ub_max:
                    lines.append(f"{'Bounds:':<20} ({lb_min:.4g}, {ub_max:.4g})")
                else:
                    lines.append(f"{'Lower bounds:':<20} [{lb_min:.4g}, {lb_max:.4g}]")
                    lines.append(f"{'Upper bounds:':<20} [{ub_min:.4g}, {ub_max:.4g}]")
        except:
            pass

        # Variable types
        try:
            if len(self.variables_type) > 0:
                type_counts = {}
                for vt in self.variables_type:
                    type_name = vt.name if hasattr(vt, "name") else str(vt)
                    type_counts[type_name] = type_counts.get(type_name, 0) + 1

                if len(type_counts) == 1:
                    type_name = list(type_counts.keys())[0]
                    lines.append(f"{'Type:':<20} {type_name}")
                else:
                    types_str = ", ".join(f"{k}: {v}" for k, v in type_counts.items())
                    lines.append(f"{'Types:':<20} {types_str}")
        except:
            pass

        lines.append("=" * 60)
        return "\n".join(lines)

    def __eq__(self, other: "VariableContainer"):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __iter__(self):

        self.__i = -1
        self.__n = len(self.sub_variables)

        return self

    def __next__(self):

        self.__i += 1

        if self.__i < self.__n:
            return self.sub_variables[self.__i]

        raise StopIteration

    def keys(self):

        return (var for var in self.sub_variables)

    def values(self):

        return ((lb, ub, var_type) for lb, ub, var_type in zip(self.lbs, self.ubs, self.variables_type))

    def items(self):

        return (
            (var, (lb, ub, var_type))
            for var, lb, ub, var_type in zip(self.sub_variables, self.lbs, self.ubs, self.variables_type)
        )

    def to_node(self):

        return Node(value=self.name, length=len(self.sub_variables))


class ConstraintContainer:

    def __init__(self, name, coefs, lbs, ubs):
        """

        Internal use only

        A container for constraints, since multiple constraints might have to be created out of a single constraint
        during problem building. It is the main object for constraint management in the linear problem object.

        For instance, the linearization of a single gpr can yield multiple rows/constraints to be added
        to the linear problem.
        Nevertheless, if one wants to replace/remove this gpr,
        a full mapping of the rows/constraints that must be replaced/removed can be found here
        and in the rows linked listed engine of the linear problem.

        :param name: the name of variable. It is used as key/identifier of the constraint in the linear problem
        :param coefs: a list of dictionaries having the coefficients for the constraint.
        That is, each dictionary in the list stands for a row in the linear problem. Dictionaries of coefficients
        must contain variable identifier value/coef pairs
        :param lbs: a list of the lower bounds of all coefficients
        :param ubs: a list of the upper bounds of all coefficients
        """

        if not name:
            name = None

        if not coefs:
            coefs = []

        if not lbs:
            lbs = []

        if not ubs:
            ubs = []

        self.name = name
        self.coefs = coefs
        self.lbs = lbs
        self.ubs = ubs

    def __len__(self):
        return len(self.coefs)

    def __str__(self):
        return f"Constraint {self.name}"

    def __repr__(self):
        """Rich representation showing constraint container details."""
        lines = []
        lines.append("=" * 60)
        lines.append(f"ConstraintContainer: {self.name}")
        lines.append("=" * 60)

        # Number of constraints/rows
        try:
            constraint_count = len(self.coefs)
            lines.append(f"{'Constraints:':<20} {constraint_count}")
        except:
            pass

        # Bounds summary
        try:
            if len(self.lbs) > 0:
                lb_min = min(self.lbs)
                lb_max = max(self.lbs)
                ub_min = min(self.ubs)
                ub_max = max(self.ubs)

                if lb_min == lb_max and ub_min == ub_max:
                    lines.append(f"{'Bounds:':<20} ({lb_min:.4g}, {ub_max:.4g})")
                else:
                    lines.append(f"{'Lower bounds:':<20} [{lb_min:.4g}, {lb_max:.4g}]")
                    lines.append(f"{'Upper bounds:':<20} [{ub_min:.4g}, {ub_max:.4g}]")
        except:
            pass

        # Coefficient statistics
        try:
            if len(self.coefs) > 0:
                # Count total variables involved
                all_vars = set()
                for coef_dict in self.coefs:
                    all_vars.update(coef_dict.keys())

                lines.append(f"{'Variables:':<20} {len(all_vars)}")

                # Average non-zero coefficients per constraint
                avg_nonzero = sum(len(c) for c in self.coefs) / len(self.coefs)
                lines.append(f"{'Avg non-zero/row:':<20} {avg_nonzero:.1f}")
        except:
            pass

        lines.append("=" * 60)
        return "\n".join(lines)

    def __eq__(self, other: "ConstraintContainer"):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __iter__(self):

        self.__i = -1
        self.__n = len(self.coefs)

        return self

    def __next__(self):

        self.__i += 1

        if self.__i < self.__n:
            return self.coefs[self.__i]

        raise StopIteration

    def keys(self):

        return (i for i in range(len(self.coefs)))

    def values(self):

        return ((coef, lb, ub) for coef, lb, ub in zip(self.coefs, self.lbs, self.ubs))

    def items(self):

        return ((i, (coef, lb, ub)) for i, (coef, lb, ub) in enumerate(zip(self.coefs, self.lbs, self.ubs)))

    def to_node(self):

        return Node(value=self.name, length=len(self.coefs))


def concat_constraints(constraints: Iterable[ConstraintContainer], name: str = None):
    """
    Internal use only.

    Concatenates a list of constraints into a single constraint container.
    :param constraints: a list of constraint containers
    :param name: the name of the new constraint container
    :return: a new constraint container
    """
    coefs = []
    lbs = []
    ubs = []

    for cnt in constraints:
        coefs.extend(cnt.coefs)
        lbs.extend(cnt.lbs)
        ubs.extend(cnt.ubs)

    return ConstraintContainer(name=name, coefs=coefs, lbs=lbs, ubs=ubs)
