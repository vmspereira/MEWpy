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
Kinetic Modeling Module

Authors: Vitor Pereira
##############################################################################
"""
import re
import warnings
from collections import OrderedDict
from typing import Any, Dict, List

import numexpr as ne
import numpy as np

from mewpy.util.parsing import Arithmetic, Latex, build_tree
from mewpy.util.utilities import AttrDict


class Compartment(object):
    """Class for modeling compartments."""

    def __init__(self, comp_id: str, name: str = None, external: bool = False, size: float = 1.0):
        """Initialize a Compartment.

        Args:
            comp_id (str): A valid unique identifier
            name (str): Compartment name (optional)
            external (bool): Is external (default: False)
            size (float): Compartment size (default: 1.0)
        """
        self.id = comp_id
        self.name = name if name is not None else comp_id
        self.size = size
        self.external = external
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


class Metabolite(object):
    """Class for modeling metabolites."""

    def __init__(self, met_id: str, name: str = None, compartment: str = None):
        """Initialize a Metabolite.

        Args:
            met_id (str): A valid unique identifier
            name (str): Common metabolite name
            compartment (str): Compartment containing the metabolite
        """
        self.id = met_id
        self.name = name if name is not None else met_id
        self.compartment = compartment
        self.metadata = OrderedDict()

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


def calculate_yprime(y, rate: np.array, stoichiometry: Dict[str, float]):
    """Calculate the rate of change for each metabolite.

    Applies stoichiometric coefficients to the reaction rate for each metabolite.
    Negative coefficients indicate substrates, positive indicate products.

    Args:
        y: Dictionary of metabolite concentrations
        rate: The calculated reaction rate
        stoichiometry: Dictionary mapping metabolite IDs to stoichiometric coefficients
                       (negative for substrates, positive for products)

    Returns:
        Dictionary of metabolite rates (y_prime) after applying stoichiometric coefficients
    """
    y_prime = {name: 0 for name in y.keys()}
    for m_id, coeff in stoichiometry.items():
        if m_id in y_prime:
            y_prime[m_id] += coeff * rate

    return y_prime


def check_positive(y_prime: List[float]) -> List[float]:
    """
    Check that substrate values are not negative when they shouldn't be.

    Returns a new list with negative values replaced by zero.
    Does not mutate the input list.

    Args:
        y_prime: List of substrate values

    Returns:
        New list with non-negative values
    """
    return [max(0, val) for val in y_prime]


class Rule(object):
    """Base class for kinetic rules."""

    def __init__(self, r_id: str, law: str, parameters: Dict[str, float] = None):
        """Initialize a Rule.

        Args:
            r_id (str): Reaction/rule identifier
            law (str): The rule string representation
            parameters (Dict[str, float], optional): Parameter values. Defaults to None.
        """
        self.id = r_id
        self.law = law
        self._tree = None
        self.parameters = parameters if parameters is not None else {}

    @property
    def tree(self):
        """Parsing tree of the law.

        Returns:
            Node: Root node of the parsing tree.
        """
        if not self._tree:
            self._tree = build_tree(self.law, Arithmetic)
        return self._tree

    def parse_parameters(self):
        """Returns the list of parameters within the rule.

        Returns:
            list: parameters
        """
        return list(self.tree.get_parameters())

    def rename_parameter(self, old_parameter: str, new_parameter: str):
        """Need to rename in the law, tree and parameters"""

        if old_parameter in self.parse_parameters():
            t = self.tree.replace({old_parameter: new_parameter})
            self.law = t.to_infix()
            self._tree = None
        try:
            self.parameters[new_parameter] = self.parameters[old_parameter]
            del self.parameters[old_parameter]
        except KeyError:
            # Parameter doesn't exist in parameters dict, that's OK
            pass

    def get_parameters(self):
        return self.parameters

    def replace(self, parameters: Dict[str, Any] = None, local: bool = True, infix: bool = True, latex: bool = False):
        """Replaces parameters with values taken from a dictionary.
        If no parameter are given for replacement, returns the string representation of the rule
        built from the parsing tree.

        Args:
            parameters (dict, optional): Replacement dictionary. Defaults to None.
            local (bool, optional): use parameter values defined in the rule
        Returns:
            str: the kinetic rule.
        """
        param = parameters.copy() if parameters else dict()
        if local:
            param.update(self.parameters)
        t = self.tree.replace(param)
        if latex:
            return Latex(t.to_latex()[0])
        elif infix:
            return t.to_infix()
        else:
            return t

    def calculate_rate(self, substrates={}, parameters={}):
        param = {**self.parameters, **substrates, **parameters}

        if len(param.keys()) != len(self.parse_parameters()):
            s = set(self.parse_parameters()) - set(param.keys())
            raise ValueError(f"Values missing for parameters: {s}")
        t = self.replace(param)

        # Convert pow(x, y) to x**y for numexpr compatibility
        t = re.sub(r"pow\s*\(\s*([^,]+)\s*,\s*([^)]+)\s*\)", r"(\1)**(\2)", t)

        # Use numexpr for safe evaluation (prevents code injection)
        try:
            rate = ne.evaluate(t, local_dict={}).item()
        except Exception as e:
            raise ValueError(f"Failed to evaluate rate expression '{t}': {e}")
        return rate

    def __str__(self):
        return self.law.replace(" ", "")

    def __repr__(self):
        return self.replace().replace(" ", "")

    def _repr_latex_(self):
        """Generate LaTeX representation for Jupyter display.

        Returns:
            str: LaTeX formatted string with $$ delimiters, or None if generation fails
        """
        try:
            s, _ = self.tree.to_latex()
            if s:
                # Jupyter needs $$ delimiters to recognize and render LaTeX as math
                return "$$ %s $$" % s
            return None
        except Exception as e:
            # If latex generation fails, return None to fall back to __repr__
            import warnings

            warnings.warn(
                f"Failed to generate LaTeX representation for {self.id}: {type(e).__name__}: {e}", RuntimeWarning
            )
            return None


class KineticReaction(Rule):

    def __init__(
        self,
        r_id: str,
        law: str,
        name: str = None,
        stoichiometry: dict = None,
        parameters: dict = None,
        modifiers: list = None,
        reversible: bool = True,
        functions: dict = None,
    ):
        """Initialize a KineticReaction.

        Args:
            r_id (str): Reaction identifier
            law (str): Kinetic law expression
            name (str, optional): The name of the reaction. Defaults to None.
            stoichiometry (dict, optional): The stoichiometry of the reaction. Defaults to None.
            parameters (dict, optional): Local parameters. Defaults to None.
            modifiers (list, optional): Reaction modifiers. Defaults to None.
            reversible (bool, optional): Reversibility. Defaults to True.
            functions (dict, optional): Functions defined in the model. Defaults to None.
        """
        super(KineticReaction, self).__init__(r_id, law, parameters)
        self.name = name if name else r_id
        self.stoichiometry = stoichiometry if stoichiometry is not None else {}
        self.modifiers = modifiers if modifiers is not None else []
        self.parameter_distributions = {}
        self.reversible = reversible
        self._model = None
        self.functions = {k: v[1] for k, v in functions.items()} if functions else {}

    @property
    def tree(self):
        """Parsing tree of the law.

        Returns:
            Node: Root node of the parsing tree.
        """
        if not self._tree:
            self._tree = build_tree(self.law, Arithmetic)
            self._tree.replace_nodes(self.functions)
        return self._tree

    @property
    def substrates(self):
        return [k for k, v in self.stoichiometry.items() if v < 0]

    @property
    def products(self):
        return [k for k, v in self.stoichiometry.items() if v > 0]

    def rename_parameter(self, old_parameter, new_parameter):
        super().rename_parameter(old_parameter, new_parameter)

        if old_parameter in self.modifiers:
            self.modifiers.append(new_parameter)
            self.modifiers.remove(old_parameter)

        if old_parameter in self.parameter_distributions:
            self.parameter_distributions[new_parameter] = self.parameter_distributions[old_parameter]
            del self.parameter_distributions[old_parameter]

    def set_parameter_distribution(self, param, dist):
        self.parameter_distributions[param] = dist

    def sample_parameter(self, param):
        dist = self.parameter_distributions.get(param, None)
        if not dist:
            raise ValueError(f"The parameter {param} has no associated distribution.")
        return dist.rvs()

    def parse_law(self, map: dict, functions=None, local=True):
        """Auxiliary method invoked by the model to build the ODE system.

        Args:
            map (dict): Dictionary of global paramameters replacements.

        Returns:
            str: kinetic rule.
        """
        m = {p_id: f"p['{p_id}']" for p_id in self.parameters.keys()}
        r_map = map.copy()
        r_map.update(m)

        return self.replace(r_map, local=local)

    def calculate_rate(self, substrates={}, parameters={}):
        # Build parameter dictionary with proper precedence (later values override earlier)
        param = {
            **self._model.get_concentrations(),  # Model defaults
            **self._model.get_parameters(),
            **self.parameters,  # Reaction defaults
            **substrates,  # User defined
            **parameters,
        }
        s = set(self.parse_parameters()) - set(param.keys())
        if s:
            # check for missing parameters distributions
            r = s - set(self.parameter_distributions.keys())
            if r:
                raise ValueError(f"Missing values or distributions for parameters: {r}")
            else:
                for p in s:
                    param[p] = self.parameter_distributions[p].rvs()
        t = self.replace(param)

        # Convert pow(x, y) to x**y for numexpr compatibility
        t = re.sub(r"pow\s*\(\s*([^,]+)\s*,\s*([^)]+)\s*\)", r"(\1)**(\2)", t)

        # Use numexpr for safe evaluation (prevents code injection)
        try:
            rate = ne.evaluate(t, local_dict={}).item()
        except Exception as e:
            raise ValueError(f"Failed to evaluate rate expression '{t}': {e}")
        return rate

    def reaction(self, y, substrates={}, parameters={}):
        """Calculate the reaction's contribution to metabolite rates of change.

        Args:
            y (dict): Dictionary of metabolite concentrations
            substrates (dict, optional): Substrate concentrations. Defaults to {}.
            parameters (dict, optional): Kinetic parameters. Defaults to {}.

        Returns:
            np.array: Array of metabolite rates with stoichiometric coefficients applied
        """
        rate = self.calculate_rate(substrates, parameters)
        y_prime_dic = calculate_yprime(y, rate, self.stoichiometry)
        # y_prime_dic = self.modify_product(y_prime_dic, substrate_names)
        y_prime = np.array(list(y_prime_dic.values()))
        # if not self.reversible:
        #    y_prime = check_positive(y_prime)

        return y_prime

    def set_parameter_defaults_to_mean(self):
        """Sets not defined parameters to the median of a distribution."""
        for name in self.parameter_distributions:
            if name not in self.parameters:
                if isinstance(self.parameter_distributions[name], (list, tuple)):
                    self.parameters[name] = (
                        self.parameter_distributions[name][0] + self.parameter_distributions[name][1]
                    ) / 2
                else:
                    self.parameters[name] = self.parameter_distributions[name].mean()


class ODEModel:
    """ODE-based kinetic model for metabolic systems."""

    def __init__(self, model_id):
        """Initialize an ODEModel.

        Args:
            model_id: Unique identifier for the model
        """
        self.id = model_id
        self.metabolites = OrderedDict()
        self.compartments = OrderedDict()
        # Kinetic rule of each reaction
        self.ratelaws = OrderedDict()
        # Initial concentration of metabolites
        self.concentrations = OrderedDict()
        # Parameter defined as constants
        self.constant_params = OrderedDict()
        # Variable parameters
        self.variable_params = OrderedDict()
        self.assignment_rules = OrderedDict()
        self.function_definition = OrderedDict()

        self._func_str = None
        self._constants = None
        self._m_r_lookup = None

    def __repr__(self):
        """Rich representation showing kinetic model details."""
        lines = []
        lines.append("=" * 60)
        lines.append(f"ODEModel: {self.id}")
        lines.append("=" * 60)

        # Model type
        lines.append(f"{'Type:':<20} ODE-based kinetic model")

        # Metabolites count
        try:
            met_count = len(self.metabolites)
            if met_count > 0:
                lines.append(f"{'Metabolites:':<20} {met_count}")
        except:
            pass

        # Reactions count
        try:
            rxn_count = len(self.ratelaws)
            if rxn_count > 0:
                lines.append(f"{'Reactions:':<20} {rxn_count}")
        except:
            pass

        # Compartments count
        try:
            comp_count = len(self.compartments)
            if comp_count > 0:
                lines.append(f"{'Compartments:':<20} {comp_count}")
        except:
            pass

        # Parameters
        try:
            const_param_count = len(self.constant_params)
            var_param_count = len(self.variable_params)
            total_params = const_param_count + var_param_count

            if total_params > 0:
                lines.append(f"{'Parameters:':<20} {total_params}")
                if const_param_count > 0:
                    lines.append(f"{'  Constant:':<20} {const_param_count}")
                if var_param_count > 0:
                    lines.append(f"{'  Variable:':<20} {var_param_count}")
        except:
            pass

        # Initial concentrations
        try:
            conc_count = len(self.concentrations)
            if conc_count > 0:
                lines.append(f"{'Init. concentrations:':<20} {conc_count}")
        except:
            pass

        # Assignment rules
        try:
            rule_count = len(self.assignment_rules)
            if rule_count > 0:
                lines.append(f"{'Assignment rules:':<20} {rule_count}")
        except:
            pass

        # Function definitions
        try:
            func_count = len(self.function_definition)
            if func_count > 0:
                lines.append(f"{'Functions:':<20} {func_count}")
        except:
            pass

        lines.append("=" * 60)
        return "\n".join(lines)

    def _repr_html_(self):
        """Pandas-like HTML representation for Jupyter notebooks."""
        from mewpy.util.html_repr import render_html_table

        rows = []
        rows.append(("Type", "ODE-based kinetic model"))

        if len(self.metabolites) > 0:
            rows.append(("Metabolites", str(len(self.metabolites))))

        if len(self.ratelaws) > 0:
            rows.append(("Reactions", str(len(self.ratelaws))))

        if len(self.compartments) > 0:
            rows.append(("Compartments", str(len(self.compartments))))

        const_param_count = len(self.constant_params)
        var_param_count = len(self.variable_params)
        total_params = const_param_count + var_param_count

        if total_params > 0:
            rows.append(("Parameters", str(total_params)))
            if const_param_count > 0:
                rows.append(("  Constant", str(const_param_count)))
            if var_param_count > 0:
                rows.append(("  Variable", str(var_param_count)))

        if len(self.concentrations) > 0:
            rows.append(("Init. concentrations", str(len(self.concentrations))))

        if len(self.assignment_rules) > 0:
            rows.append(("Assignment rules", str(len(self.assignment_rules))))

        if len(self.function_definition) > 0:
            rows.append(("Functions", str(len(self.function_definition))))

        return render_html_table(f"ODEModel: {self.id}", rows)

    def _clear_temp(self):
        self._func_str = None

    def add_compartment(self, compartment, replace=True):
        """Add a compartment to the model.
        Arguments:
            compartment (Compartment): compartment to add
            replace (bool): replace previous compartment with same id (default: True)
        """
        if compartment.id in self.compartments and not replace:
            raise RuntimeError(f"Compartment {compartment.id} already exists.")
        self.compartments[compartment.id] = compartment

    def add_metabolite(self, metabolite, replace=True):
        """Add a metabolite to the model.
        Arguments:
            metabolite (Metabolite): metabolite to add
            replace (bool): replace previous metabolite with same id (default: True)
        """

        if metabolite.id in self.metabolites and not replace:
            raise RuntimeError(f"Metabolite {metabolite.id} already exists.")

        if metabolite.compartment not in self.compartments:
            raise RuntimeError(
                f"Metabolite {metabolite.id} \
                has invalid compartment {metabolite.compartment}."
            )

        self.metabolites[metabolite.id] = metabolite

    def set_functions(self, functions):
        self.function_definition = functions

    @property
    def reactions(self):
        return AttrDict(self.ratelaws)

    def get_reaction(self, r_id):
        if r_id not in self.ratelaws:
            raise ValueError(f"Unknown reaction {r_id}")
        r = self.ratelaws[r_id]
        d = {
            "id": r_id,
            "name": r.name,
            "stoichiometry": r.stoichiometry,
            "law": r.law,
            "reversible": r.reversible,
            "parameters": r.parameters,
            "modifiers": r.modifiers,
        }
        return AttrDict(d)

    def get_metabolite(self, m_id):
        if m_id not in self.metabolites:
            raise ValueError(f"Unknown metabolite {m_id}")
        d = {
            "id": m_id,
            "name": self.metabolites[m_id].name,
            "compartment": self.metabolites[m_id].compartment,
            "formula": self.metabolites[m_id].metadata.get("FORMULA", ""),
            "charge": self.metabolites[m_id].metadata.get("CHARGE", ""),
            "y0": self.concentrations[m_id],
        }
        return AttrDict(d)

    def find(self, pattern=None, sort=False):
        """A user friendly method to find reactions in the model.

        :param pattern: The pattern which can be a regular expression,
            defaults to None in which case all entries are listed.
        :type pattern: str, optional
        :param sort: if the search results should be sorted, defaults to False
        :type sort: bool, optional
        :return: the search results
        :rtype: pandas dataframe
        """
        values = list(self.reactions.keys())
        if pattern:
            if isinstance(pattern, list):
                patt = "|".join(pattern)
            else:
                patt = pattern
            re_expr = re.compile(patt)
            values = [x for x in values if re_expr.search(x)]
        if sort:
            values.sort()

        import pandas as pd

        data = [self.get_reaction(x) for x in values]

        if data:
            df = pd.DataFrame(data)
            df = df.set_index(df.columns[0])
        else:
            df = pd.DataFrame()
        return df

    def find_reactions(self, pattern=None, sort=False):
        return self.find(pattern, sort)

    def find_metabolites(self, pattern=None, sort=False):
        """A user friendly method to find metabolites in the model.

        :param pattern: The pattern which can be a regular expression,
            defaults to None in which case all entries are listed.
        :type pattern: str, optional
        :param sort: if the search results should be sorted, defaults to False
        :type sort: bool, optional
        :return: the search results
        :rtype: pandas dataframe
        """
        values = list(self.metabolites.keys())
        if pattern:
            if isinstance(pattern, list):
                patt = "|".join(pattern)
            else:
                patt = pattern
            re_expr = re.compile(patt)
            values = [x for x in values if re_expr.search(x)]
        if sort:
            values.sort()

        import pandas as pd

        data = [self.get_metabolite(x) for x in values]

        if data:
            df = pd.DataFrame(data)
            df = df.set_index(df.columns[0])
        else:
            df = pd.DataFrame()
        return df

    def get_concentrations(self):
        return self.concentrations

    def set_concentration(self, m_id: str, concentration: float):
        """Sets a metabolite initial concentration

        Args:
            m_id (str): Metabolite identifier
            concentration (float): Initial concentration
        """
        if m_id in self.metabolites:
            self.concentrations[m_id] = concentration
        else:
            warnings.warn(f"No such metabolite '{m_id}'", RuntimeWarning)

    def set_ratelaw(self, r_id: str, law: KineticReaction):
        """Define the rate law for a given reaction.

        Args:
            r_id (str): Reaction Identifier
            law (KineticReaction): The reaction rate law.
        """
        law._model = self
        self.ratelaws[r_id] = law

    def get_ratelaw(self, r_id):
        if r_id in self.ratelaws.keys():
            return self.ratelaws[r_id]
        else:
            raise ValueError("Reaction has no rate law.")

    def set_assignment_rule(self, p_id: str, rule: Rule):
        if p_id in self.variable_params or p_id in self.metabolites:
            self.assignment_rules[p_id] = rule
        else:
            warnings.warn(f"No such variable parameter '{p_id}'", RuntimeWarning)

    def set_global_parameter(self, key, value, constant=True):
        if constant:
            self.constant_params[key] = value
        else:
            self.variable_params[key] = value

    def merge_constants(self):
        constants = OrderedDict()

        for c_id, comp in self.compartments.items():
            constants[c_id] = comp.size

        constants.update(self.constant_params)

        for r_id, law in self.ratelaws.items():
            for p_id, value in law.parameters.items():

                full_id = f"{p_id}"
                if full_id in constants:
                    warnings.warn(f"renaming {p_id} to {r_id}_{p_id}")
                    full_id = f"{r_id}_{p_id}"
                    law.rename_parameter(f"{p_id}", f"{r_id}_{p_id}")
                constants[full_id] = value

        self._constants = constants
        return constants

    def metabolite_reaction_lookup(self):
        if not self._m_r_lookup:
            self._m_r_lookup = {m_id: {} for m_id in self.metabolites}

            for r_id, rule in self.ratelaws.items():
                for m_id, coeff in rule.stoichiometry.items():
                    self._m_r_lookup[m_id][r_id] = coeff
        return self._m_r_lookup

    def print_balance(self, m_id, factors=None):
        """Returns a string representation of the mass balance equation
           of a metabolite

        Args:
            m_id (str): The metabolite identifier
            factors (dict, optional): Factors applied to metabolite production (products only).
                                     Used to model reduced enzyme expression or regulatory effects.
                                     Defaults to None.

        Returns:
            str: Mass balance equation
        """
        f = factors.get(m_id, 1) if factors else 1
        c_id = self.metabolites[m_id].compartment
        table = self.metabolite_reaction_lookup()

        terms = []
        for r_id, coeff in table[m_id].items():
            # Apply factor only to products (positive coefficients) to model reduced production
            # while keeping consumption (negative coefficients) unchanged
            v = coeff * f if coeff > 0 else coeff
            terms.append(f"{v:+g} * r['{r_id}']")

        # Check if metabolite is constant and boundary (attributes may not exist)
        is_constant = getattr(self.metabolites[m_id], "constant", False)
        is_boundary = getattr(self.metabolites[m_id], "boundary", False)
        if f == 0 or len(terms) == 0 or (is_constant and is_boundary):
            expr = "0"
        else:
            expr = f"1/p['{c_id}'] * ({' '.join(terms)})"
        return expr

    def get_parameters(self, exclude_compartments=False):
        """Returns a dictionary of the model parameters"""
        if not self._constants:
            self.merge_constants()
        parameters = self._constants.copy()
        if exclude_compartments:
            for c_id in self.compartments:
                del parameters[c_id]
        return parameters

    def find_parameters(self, pattern=None, sort=False):
        """A user friendly method to find parameters in the model.

        :param pattern: The pattern which can be a regular expression,
            defaults to None in which case all entries are listed.
        :type pattern: str, optional
        :param sort: if the search results should be sorted, defaults to False
        :type sort: bool, optional
        :return: the search results
        :rtype: pandas dataframe
        """
        params = self.get_parameters()
        values = list(params.keys())
        if pattern:
            if isinstance(pattern, list):
                patt = "|".join(pattern)
            else:
                patt = pattern
            re_expr = re.compile(patt)
            values = [x for x in values if re_expr.search(x)]
        if sort:
            values.sort()

        import pandas as pd

        data = [(x, params[x]) for x in values]

        if data:
            df = pd.DataFrame(data, columns=["Parameter", "Value"])
            df = df.set_index(df.columns[0])
        else:
            df = pd.DataFrame()
        return df

    def find_functions(self, pattern=None, sort=False):
        """A user friendly method to find functions in the model.

        :param pattern: The pattern which can be a regular expression,
            defaults to None in which case all entries are listed.
        :type pattern: str, optional
        :param sort: if the search results should be sorted, defaults to False
        :type sort: bool, optional
        :return: the search results
        :rtype: pandas dataframe
        """
        params = self.function_definition
        values = list(params.keys())
        if pattern:
            if isinstance(pattern, list):
                patt = "|".join(pattern)
            else:
                patt = pattern
            re_expr = re.compile(patt)
            values = [x for x in values if re_expr.search(x)]
        if sort:
            values.sort()

        import pandas as pd

        data = [(x, ",".join(params[x][0]), str(params[x][1])) for x in values]

        if data:
            df = pd.DataFrame(data, columns=["Name", "Arguments", "Body"])
            df = df.set_index(df.columns[0])
        else:
            df = pd.DataFrame()
        return df

    def deriv(self, t, y):
        """
        Deriv function called by integrate.

        For each step when the model is run, the rate for each reaction is calculated
        and changes in substrates and products calculated, normalized by compartment volume.

        Args:
            t : time, not used in this function but required
            y (list): ordered list of substrate values (the order
            is the same as the metabolites order) at this current timepoint.
               Has the same order as self.run_model_species_names

        Returns:
            y_prime - ordered list the same as y, y_prime is the new set of y's for this timepoint,
                     normalized by compartment volumes.
        """
        p = self.merge_constants()
        m_y = OrderedDict(zip(self.metabolites, y))
        yprime = np.zeros(len(y))
        for _, reaction in self.ratelaws.items():
            yprime += reaction.reaction(m_y, self.get_parameters(), p)

        # Normalize by compartment volume (dC/dt = rate / volume)
        for i, m_id in enumerate(self.metabolites):
            c_id = self.metabolites[m_id].compartment
            volume = p[c_id]
            yprime[i] /= volume

        return yprime.tolist()

    def build_ode(self, factors: dict = None, local: bool = False) -> str:
        """
        Auxiliary function to build the ODE as a string
        to be evaluated by eval, as an alternative to deriv.
        Allows the inclusion of factors to be applied to the parameters

        Args:
            factors (dict): factors to be applied to parameters
            local (bool): enforces the usage of parameter values defined within
            the reactions
        Returns:
            func_str the right-hand side of the system.
        """

        m = {m_id: f"x[{i}]" for i, m_id in enumerate(self.metabolites)}
        c = {c_id: f"p['{c_id}']" for c_id in self.compartments}
        p = {p_id: f"p['{p_id}']" for p_id in self.constant_params}
        v = {p_id: f"v['{p_id}']" for p_id in self.variable_params}
        rmap = OrderedDict({**m, **c, **p, **v})

        parsed_rates = {r_id: ratelaw.parse_law(rmap, local=local) for r_id, ratelaw in self.ratelaws.items()}

        r = {r_id: f"({parsed_rates[r_id]})" for r_id in self.ratelaws.keys()}

        rmap.update(r)

        rate_exprs = [" " * 4 + "r['{}'] = {}".format(r_id, parsed_rates[r_id]) for r_id in self.ratelaws.keys()]

        # Build mass balance equations for each metabolite
        balances = [" " * 8 + self.print_balance(m_id, factors=factors) for m_id in self.metabolites]

        func = "def ode_func(t, x, r, p, v)"
        func_str = (
            func
            + ":\n\n"
            + "\n".join(rate_exprs)
            + "\n\n"
            + "    dxdt = [\n"
            + ",\n".join(balances)
            + "\n"
            + "    ]\n\n"
            + "    return dxdt\n"
        )

        self._func_str = func_str
        return self._func_str

    def get_ode(self, r_dict=None, params=None, factors=None):
        """
        Args:
            r_dict: reaction identifiers to be modified
            params: modified parameters
            factors: factors to be applied to parameters
        """
        p = self.merge_constants()
        if params:
            p.update(params)

        if factors is not None:
            for k, v in factors.items():
                if k in p.keys():
                    p[k] = v * p[k]

        v = self.variable_params.copy()
        r = r_dict if r_dict is not None else dict()

        np.seterr(divide="ignore", invalid="ignore")
        # Use local namespace instead of globals() to prevent pollution and security issues
        local_namespace = {}
        exec(self.build_ode(factors), local_namespace)
        ode_func = local_namespace["ode_func"]

        return lambda t, y: ode_func(t, y, r, p, v)
