from io import TextIOWrapper
from typing import TYPE_CHECKING, Dict, Sequence, Tuple, Union

from cobra import Model as Cobra_Model
from reframed import CBModel as Reframed_Model

from mewpy.io import Engines, Reader, read_model

from .problem import AbstractKOProblem

if TYPE_CHECKING:
    from mewpy.germ.models import MetabolicModel, Model, RegulatoryModel
    from mewpy.optimization import EvaluationFunction


def load_optorf(
    regulatory_model: Union[str, TextIOWrapper, Reader],
    metabolic_model: Union[str, TextIOWrapper, Cobra_Model, Reframed_Model, Reader],
    config: dict = None,
    warnings: bool = False,
):
    """
    The standard method to load an OptORF problem.
    A OptORF problem is a KO strain optimization problem
    that can be created from an integrated Metabolic-Regulatory model.

    The integrated model can be assembled from a metabolic model in the following formats:
        - SBML
        - JSON
        - Cobrapy model
        - Reframed model

    and a regulatory model in the following formats:
        - SBML-qual
        - JSON
        - CSV regulatory network

    :param regulatory_model: the regulatory model or a reader for parsing it
    :param metabolic_model: the metabolic model or a reader for parsing it
    :param config: A configuration dictionary with information for parsing the models
    :param warnings: if True, warnings will be printed
    :return: an integrated model for OptORF
    """
    if not config:
        config = {}

    if not isinstance(regulatory_model, Reader):

        if isinstance(regulatory_model, str):

            file_name = regulatory_model

        elif isinstance(regulatory_model, TextIOWrapper):

            file_name = regulatory_model.name

        else:
            raise ImportError("Invalid file type")

        engine = Engines.BooleanRegulatoryCSV

        if file_name.endswith(".xml") or file_name.endswith(".sbml"):
            engine = Engines.RegulatorySBML

        regulatory_model = Reader(engine=engine, io=regulatory_model, **config)

    if not isinstance(metabolic_model, Reader):

        if isinstance(metabolic_model, str):

            engine = Engines.MetabolicSBML

        elif isinstance(metabolic_model, TextIOWrapper):

            engine = Engines.MetabolicSBML

        elif isinstance(metabolic_model, Cobra_Model):

            engine = Engines.CobraModel

        elif isinstance(metabolic_model, Reframed_Model):

            engine = Engines.ReframedModel

        else:
            raise ImportError("Invalid file type")

        metabolic_model = Reader(engine=engine, io=metabolic_model, **config)

    return read_model(regulatory_model, metabolic_model, warnings=warnings)


class OptORFProblem(AbstractKOProblem):

    # regulators that are not genes are left unchanged by the gene resolution
    _target_kind = "gene"

    def __init__(
        self,
        model: Union["Model", "MetabolicModel", "RegulatoryModel"],
        fevaluation: Sequence["EvaluationFunction"],
        initial_state: Dict[str, float] = None,
        **kwargs,
    ):
        """
        OptORF problem using the RFBA implementation analysis.
        The OptORF approach is based on gene and regulator deletion to identify optimization strategies.

        It currently supports KOProblem-based decoding

        For more details consult: https://doi.org/10.1186/1752-0509-4-53
        """
        if isinstance(model, (Cobra_Model, Reframed_Model)):
            raise ValueError(
                f"OptORF is not available for a model of type {type(model)}."
                f"Please use load_optorf() to retrieve an integrated GERM model"
            )

        super(OptORFProblem, self).__init__(model, fevaluation, **kwargs)

        # ignore user defined target list
        self._trg_list = None

        if not initial_state:
            initial_state = {}
        self._initial_state = initial_state.copy()

    def _build_target_list(self):
        """
        Build the target list and initial state for the OptORF problem.
        The target list is the list of regulators in the regulatory layer.
        The initial state is the state of the model before any gene deletion.
        :return:
        """
        # Target list is the combination of genes and regulators available into the mewpy integrated model

        regulators = [
            regulator.id
            for regulator in self.model.yield_regulators()
            if not regulator.is_reaction() and not regulator.is_metabolite()
        ]
        targets = [
            target.id
            for target in self.model.yield_targets()
            if not target.is_reaction() and not target.is_metabolite()
        ]
        genes = [gene.id for gene in self.model.yield_genes() if not gene.is_reaction() and not gene.is_metabolite()]

        # sorted: set order depends on the interpreter hash seed, which breaks seeded runs
        self._trg_list = sorted(set.union(set(regulators), set(targets), set(genes)))

    def _rfba(self):
        """The RFBA of the problem's model, which decodes the candidates (not pickled; rebuilt where needed)"""
        rfba = self.__dict__.get("_rfba_cache")
        if rfba is None:
            from mewpy.germ.analysis import RFBA

            rfba = RFBA(self.model)
            self._rfba_cache = rfba
        return rfba

    def __getstate__(self):
        state = self.__dict__.copy()
        state.pop("_rfba_cache", None)
        return state

    def solution_to_constraints(
        self,
        decoded_candidate: Dict[str, float],
    ) -> Dict[str, Tuple[float, float]]:
        """
        Converts a decoded candidate, the regulators/genes knocked out ({id: 0}), into constraints: the reactions
        that RFBA closes with the problem's initial state and the knockouts (see RFBA.regulatory_constraints). The
        regulators that follow rules settle from that state, and the knockouts win over the decoded states, also for
        regulators and genes that are targets of the regulatory network.

        :param decoded_candidate: the regulators/genes to be knocked out, as returned by decode
        :return: a metabolic state (a set of constraints)
        """
        knockouts = {k: 0.0 for k in decoded_candidate}
        return self._rfba().regulatory_constraints({**self._initial_state, **knockouts})
