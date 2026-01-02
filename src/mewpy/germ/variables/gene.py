from typing import TYPE_CHECKING, Any, Dict, Generator, List, Sequence, Set, Tuple, Union

from mewpy.germ.models.serialization import serialize
from mewpy.util.constants import ModelConstants
from mewpy.util.history import recorder
from mewpy.util.utilities import generator

from .variable import Variable
from .variables_utils import coefficients_setter, initialize_coefficients

if TYPE_CHECKING:
    from .reaction import Reaction


class Gene(Variable, variable_type="gene", register=True, constructor=True, checker=True):

    def __init__(
        self, identifier: Any, coefficients: Sequence[float] = None, reactions: Dict[str, "Reaction"] = None, **kwargs
    ):
        """
        A metabolic gene is regularly associated with metabolic reactions. A metabolic gene is inferred from a metabolic
        model by parsing the Gene-Protein-Reactions associations usually encoded as boolean rules.

        A gene contains multiple coefficients for GPR rule evaluation. It usually takes
        either 0 (inactive) or 1 (active)
        A gene object also holds information regarding the reactions that it is associated with.

        In addition, a metabolic gene can also be represented as a target gene in Metabolic-Regulatory models, as it
        can be inferred from the regulatory rules defining regulatory interactions.

        :param identifier: identifier, e.g. b0001
        :param coefficients: the set of coefficients that this gene can take.
        These coefficients can be expanded later. 0 and 1 are added by default
        :param reactions: the dictionary of reactions to which the gene is associated with
        """
        # Initialize coefficients with defaults (0.0, 1.0) if not provided
        coefficients = initialize_coefficients(coefficients)

        if not reactions:
            reactions = {}

        self._coefficients = coefficients
        self._reactions = reactions

        super().__init__(identifier, **kwargs)

    @property
    def types(self):
        # noinspection PyUnresolvedReferences
        _types = {Gene.variable_type}

        _types.update(super(Gene, self).types)

        return _types

    def __repr__(self):
        """Rich representation showing gene details."""
        lines = []
        lines.append("=" * 60)
        lines.append(f"Gene: {self.id}")
        lines.append("=" * 60)

        # Name if available
        if hasattr(self, "name") and self.name and self.name != self.id:
            lines.append(f"{'Name:':<20} {self.name}")

        # Activity status
        try:
            status = "Active" if self.is_active else "Inactive"
            lines.append(f"{'Status:':<20} {status}")
        except:
            pass

        # Coefficients
        try:
            coef = self.coefficients
            if len(coef) <= 3:
                coef_str = ", ".join(f"{c:.4g}" for c in coef)
            else:
                coef_str = f"{len(coef)} values: [{coef[0]:.4g}, ..., {coef[-1]:.4g}]"
            lines.append(f"{'Coefficients:':<20} {coef_str}")
        except:
            pass

        # Reactions count
        try:
            rxn_count = len(self.reactions)
            lines.append(f"{'Reactions:':<20} {rxn_count}")
        except:
            pass

        lines.append("=" * 60)
        return "\n".join(lines)

    def _repr_html_(self):
        """Pandas-like HTML representation for Jupyter notebooks."""
        from mewpy.util.html_repr import render_html_table

        rows = []

        # Name
        if hasattr(self, "name") and self.name and self.name != self.id:
            rows.append(("Name", self.name))

        # Aliases
        try:
            if hasattr(self, "aliases") and self.aliases:
                aliases_list = sorted(list(self.aliases))
                if len(aliases_list) <= 5:
                    aliases_str = ", ".join(aliases_list)
                else:
                    aliases_str = f"{', '.join(aliases_list[:5])}, ... ({len(aliases_list)} total)"
                rows.append(("Aliases", aliases_str))
        except:
            pass

        # Activity status
        try:
            status = "Active" if self.is_active else "Inactive"
            rows.append(("Status", status))
        except:
            pass

        # Coefficients
        try:
            coef = self.coefficients
            if len(coef) <= 3:
                coef_str = ", ".join(f"{c:.4g}" for c in coef)
            else:
                coef_str = f"{len(coef)} values: [{coef[0]:.4g}, ..., {coef[-1]:.4g}]"
            rows.append(("Coefficients", coef_str))
        except:
            pass

        # Reactions (detailed list)
        try:
            reactions = self.reactions
            if reactions:
                rxn_ids = sorted(list(reactions.keys()))
                if len(rxn_ids) <= 5:
                    rxn_str = ", ".join(rxn_ids)
                else:
                    rxn_str = f"{', '.join(rxn_ids[:5])}, ... ({len(rxn_ids)} total)"
                rows.append(("Reactions", rxn_str))
        except:
            pass

        return render_html_table(f"Gene: {self.id}", rows)

    def _gene_to_html_(self):
        """
        It returns a html dict representation.
        """
        html_dict = {
            "Coefficients": self.coefficients,
            "Active": self.is_active,
            "Reactions": ", ".join(self.reactions),
        }
        return html_dict

    # -----------------------------------------------------------------------------
    # Static attributes
    # -----------------------------------------------------------------------------
    @serialize("coefficients", "coefficients", "_coefficients")
    @property
    def coefficients(self) -> Tuple[float, ...]:
        """
        The gene coefficients
        :return: The gene coefficients
        """
        if hasattr(self, "_bounds"):

            # if it is a reaction, bounds must be returned
            return self._bounds

        # if it is a metabolite, the bounds coefficient of the exchange reaction must be returned
        elif hasattr(self, "exchange_reaction"):

            if hasattr(self.exchange_reaction, "_bounds"):

                # noinspection PyProtectedMember
                return self.exchange_reaction._bounds

        return self._coefficients

    @serialize("reactions", "reactions", "_reactions")
    @property
    def reactions(self) -> Dict[str, "Reaction"]:
        """
        The reactions associated with this gene
        """
        return self._reactions.copy()

    @property
    def is_active(self):
        """
        It checks whether the gene is active or not
        :return: True if the gene is active, False otherwise
        """
        return max(self.coefficients) > ModelConstants.TOLERANCE

    # -----------------------------------------------------------------------------
    # Static attributes setters
    # -----------------------------------------------------------------------------
    @coefficients.setter
    @recorder
    def coefficients(self, value: Union[float, Sequence[float]]):
        """
        The gene coefficients setter
        :param value: The gene coefficients
        :return:
        """
        coefficients_setter(self, value)

    @reactions.setter
    def reactions(self, value: Dict[str, "Reaction"]):
        """
        The reactions associated with this gene
        :param value: The reactions associated with this gene
        :return:
        """
        if not value:
            value = {}

        self._reactions = value

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------
    def yield_reactions(self) -> Generator["Reaction", None, None]:
        """
        It yields all reactions
        """
        return generator(self._reactions)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------
    def ko(self, minimum_coefficient: float = 0.0, history=True):
        """
        It performs a gene knock-out. This is accomplished by setting the coefficient object to a tuple of zeros.

        This operation can be reverted using the model history

        :param minimum_coefficient: The minimum coefficient that represents a gene KO. 0.0 is used as default
        :param history: Whether to register this operation in the model history
        :return: None
        """
        old_coef = tuple(self.coefficients)

        coefficients_setter(self, minimum_coefficient)

        if history:
            self.history.queue_command(
                undo_func=coefficients_setter,
                undo_kwargs={"instance": self, "value": old_coef},
                func=self.ko,
                kwargs={"minimum_coefficient": minimum_coefficient, "history": False},
            )
        return

    def update(
        self,
        coefficients: Union[Set[Union[int, float]], List[Union[int, float]], Tuple[Union[int, float]]] = None,
        reactions: Dict[str, "Reaction"] = None,
        **kwargs,
    ):
        """
        It performs an update operation to this gene.
        The update operation is similar to a dictionary update.

        Note that, some update operations are not registered in history.
        It is strongly advisable to use update outside history context manager

        :param coefficients: The gene coefficients
        :param reactions: The reactions associated with this gene
        :param kwargs: Other arguments for the base variable, such as identifier, name, etc
        :return:
        """
        super(Gene, self).update(**kwargs)

        if coefficients is not None:
            self.coefficients = coefficients

        if reactions is not None:
            self._reactions.update(reactions)
