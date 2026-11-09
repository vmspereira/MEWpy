"""
RegulatoryExtension: Decorator that adds regulatory network capabilities to Simulators.

This module provides the RegulatoryExtension class, which wraps external simulators
(COBRApy, reframed) and adds regulatory network functionality without duplicating
metabolic data. All metabolic operations are delegated to the wrapped simulator.
"""

import json
import logging
from typing import TYPE_CHECKING, Any, Dict, Generator, List, Optional, Tuple, Union

if TYPE_CHECKING:
    from mewpy.germ.algebra import Expression
    from mewpy.germ.variables import Interaction, Regulator, Target
    from mewpy.simulation.simulation import Simulator

# Runtime imports
from mewpy.germ.algebra import parse_expression
from mewpy.germ.models.model import notify_observer

logger = logging.getLogger(__name__)

# id kinds tried by the TRN reconciliation, in order, and the simulator container of each kind
_ID_KINDS = (("gene", "genes"), ("reaction", "reactions"), ("metabolite", "metabolites"))


class RegulatoryExtension:
    """
    Decorator that extends Simulator with regulatory network capabilities.

    This class wraps a Simulator instance and adds regulatory network functionality
    while delegating ALL metabolic operations to the underlying simulator. This ensures
    no duplication of metabolic data and maintains a clean separation between
    metabolic (external) and regulatory (GERM) concerns.

    Key principles:
    - Stores ONLY regulatory network (regulators, targets, interactions)
    - Delegates ALL metabolic operations to wrapped simulator
    - No GERM metabolic variables created
    - Caches parsed GPR expressions for performance

    Example:
        >>> import cobra
        >>> from mewpy.simulation import get_simulator
        >>> from mewpy.germ.models import RegulatoryExtension, RegulatoryModel
        >>>
        >>> # Load metabolic model
        >>> cobra_model = cobra.io.read_sbml_model('ecoli_core.xml')
        >>> simulator = get_simulator(cobra_model)
        >>>
        >>> # Add regulatory network
        >>> from mewpy.io import Engines, Reader, read_model
        >>> reader = Reader(Engines.BooleanRegulatoryCSV, 'regulatory.csv', sep=',')
        >>> reg_network = read_model(reader)
        >>> integrated = RegulatoryExtension(simulator, reg_network)
        >>>
        >>> # Use in analysis
        >>> from mewpy.germ.analysis import RFBA
        >>> rfba = RFBA(integrated)
        >>> solution = rfba.optimize()
    """

    def __init__(
        self, simulator: "Simulator", regulatory_network: Optional[Any] = None, identifier: Optional[str] = None
    ):
        """
        Initialize RegulatoryExtension with a simulator and optional regulatory network.

        :param simulator: External simulator instance (COBRApy, reframed, etc.)
        :param regulatory_network: RegulatoryModel instance or None
        :param identifier: Optional identifier for the integrated model
        """
        self._simulator = simulator
        self._identifier = identifier or getattr(simulator, "id", "regulatory_extension")

        # Regulatory network storage (ONLY regulatory components)
        self._regulators: Dict[str, "Regulator"] = {}
        self._targets: Dict[str, "Target"] = {}
        self._interactions: Dict[str, "Interaction"] = {}

        # Cache for parsed GPR expressions (performance optimization)
        self._gpr_cache: Dict[str, "Expression"] = {}

        # Observers (analysis methods attached with attach=True) and TRN id -> simulator id renames
        self._simulators: List[Any] = []
        self._id_map: Dict[str, str] = {}
        self._regulatory_network = None

        # Load regulatory network if provided
        if regulatory_network is not None:
            self._load_regulatory_network(regulatory_network)

    # =========================================================================
    # Factory Methods
    # =========================================================================

    @classmethod
    def from_sbml(
        cls,
        metabolic_path: str,
        regulatory_path: str = None,
        regulatory_format: str = "csv",
        flavor: str = "reframed",
        identifier: str = None,
        **regulatory_kwargs,
    ) -> "RegulatoryExtension":
        """
        Create RegulatoryExtension from SBML metabolic model and optional regulatory network.

        This is the most convenient way to create an integrated model from files.

        :param metabolic_path: Path to SBML metabolic model file
        :param regulatory_path: Optional path to regulatory network file
        :param regulatory_format: Format of regulatory file ('csv', 'sbml', 'json'). Default: 'csv'
        :param flavor: Metabolic model library ('reframed' or 'cobra'). Default: 'reframed'
                       reframed is preferred as it is more lightweight than COBRApy.
        :param identifier: Optional identifier for the model
        :param regulatory_kwargs: Additional arguments for regulatory network reader
                                   (e.g., sep=',', id_col=0, rule_col=2 for CSV)
        :return: RegulatoryExtension instance

        Example:
            >>> # Load metabolic model only (uses reframed by default)
            >>> model = RegulatoryExtension.from_sbml('ecoli_core.xml')
            >>>
            >>> # Load with regulatory network from CSV
            >>> model = RegulatoryExtension.from_sbml(
            ...     'ecoli_core.xml',
            ...     'ecoli_core_trn.csv',
            ...     regulatory_format='csv',
            ...     sep=',',
            ...     id_col=0,
            ...     rule_col=2
            ... )
            >>>
            >>> # Use COBRApy instead if needed
            >>> model = RegulatoryExtension.from_sbml('ecoli_core.xml', flavor='cobra')
            >>>
            >>> # Use in analysis
            >>> from mewpy.germ.analysis import RFBA
            >>> rfba = RFBA(model)
            >>> solution = rfba.optimize()
        """
        from mewpy.simulation import get_simulator

        # Load metabolic model (prefer reframed - more lightweight)
        if flavor == "reframed":
            from reframed.io.sbml import load_cbmodel

            metabolic_model = load_cbmodel(metabolic_path)
        elif flavor == "cobra":
            import cobra

            metabolic_model = cobra.io.read_sbml_model(metabolic_path)
        else:
            raise ValueError(f"Unknown flavor: {flavor}. Use 'reframed' (default, lightweight) or 'cobra'")

        # Create simulator
        simulator = get_simulator(metabolic_model)

        # Load regulatory network if provided
        regulatory_network = None
        if regulatory_path:
            regulatory_network = cls._load_regulatory_from_file(regulatory_path, regulatory_format, **regulatory_kwargs)

        return cls(simulator, regulatory_network, identifier)

    @classmethod
    def from_model(
        cls,
        metabolic_model,
        regulatory_path: str = None,
        regulatory_format: str = "csv",
        identifier: str = None,
        **regulatory_kwargs,
    ) -> "RegulatoryExtension":
        """
        Create RegulatoryExtension from COBRApy/reframed model and optional regulatory network.

        :param metabolic_model: COBRApy Model or reframed CBModel instance
        :param regulatory_path: Optional path to regulatory network file
        :param regulatory_format: Format of regulatory file ('csv', 'sbml', 'json'). Default: 'csv'
        :param identifier: Optional identifier for the model
        :param regulatory_kwargs: Additional arguments for regulatory network reader
        :return: RegulatoryExtension instance

        Example:
            >>> import cobra
            >>> cobra_model = cobra.io.read_sbml_model('ecoli_core.xml')
            >>>
            >>> # Create with regulatory network
            >>> model = RegulatoryExtension.from_model(
            ...     cobra_model,
            ...     'ecoli_core_trn.csv',
            ...     sep=','
            ... )
        """
        from mewpy.simulation import get_simulator

        # Create simulator from model
        simulator = get_simulator(metabolic_model)

        # Load regulatory network if provided
        regulatory_network = None
        if regulatory_path:
            regulatory_network = cls._load_regulatory_from_file(regulatory_path, regulatory_format, **regulatory_kwargs)

        return cls(simulator, regulatory_network, identifier)

    @classmethod
    def from_json(cls, json_path: str, identifier: str = None) -> "RegulatoryExtension":
        """
        Create RegulatoryExtension from JSON file containing both metabolic and regulatory data.

        :param json_path: Path to JSON file
        :param identifier: Optional identifier for the model
        :return: RegulatoryExtension instance

        Example:
            >>> model = RegulatoryExtension.from_json('integrated_model.json')
        """
        from mewpy.io import Engines, Reader

        # Use existing JSON reader
        reader = Reader(Engines.JSON, json_path)
        from mewpy.io import read_model

        # Read model using existing infrastructure
        integrated_model = read_model(reader, warnings=False)

        # Extract simulator and regulatory network
        from mewpy.simulation import get_simulator

        simulator = get_simulator(integrated_model)

        # Create RegulatoryExtension with regulatory components
        from mewpy.germ.models.regulatory import RegulatoryModel

        regulatory_network = RegulatoryModel(
            identifier="regulatory",
            interactions=integrated_model.interactions if hasattr(integrated_model, "interactions") else {},
            regulators=integrated_model.regulators if hasattr(integrated_model, "regulators") else {},
            targets=integrated_model.targets if hasattr(integrated_model, "targets") else {},
        )

        return cls(simulator, regulatory_network, identifier)

    @staticmethod
    def _load_regulatory_from_file(file_path: str, file_format: str, **kwargs):
        """
        Load regulatory network from file.

        :param file_path: Path to regulatory network file
        :param file_format: File format ('csv', 'sbml', 'json')
        :param kwargs: Additional arguments for the reader
        :return: RegulatoryModel instance
        """
        from mewpy.io import Engines, Reader, read_model

        # Determine engine based on format
        if file_format.lower() == "csv":
            # Default to BooleanRegulatoryCSV
            engine = Engines.BooleanRegulatoryCSV
        elif file_format.lower() == "sbml":
            engine = Engines.RegulatorySBML
        elif file_format.lower() == "json":
            engine = Engines.JSON
        else:
            raise ValueError(f"Unknown regulatory format: {file_format}. Use 'csv', 'sbml', or 'json'")

        # Create reader and load model
        reader = Reader(engine, file_path, **kwargs)
        regulatory_model = read_model(reader, warnings=False)

        return regulatory_model

    # =========================================================================
    # Properties
    # =========================================================================

    @property
    def id(self) -> str:
        """Model identifier."""
        return self._identifier

    @property
    def simulator(self) -> "Simulator":
        """Access to the underlying simulator."""
        return self._simulator

    # =========================================================================
    # Metabolic Data Access (Delegated to Simulator)
    # =========================================================================

    @property
    def reactions(self):
        """List of reaction IDs from simulator."""
        return self._simulator.reactions

    @property
    def genes(self):
        """List of gene IDs from simulator."""
        return self._simulator.genes

    @property
    def metabolites(self):
        """List of metabolite IDs from simulator."""
        return self._simulator.metabolites

    @property
    def compartments(self):
        """Compartments from simulator."""
        return self._simulator.compartments

    @property
    def medium(self):
        """Medium composition from simulator."""
        return self._simulator.medium

    @property
    def objective(self):
        """Objective function from simulator."""
        return self._simulator.objective

    def get_reaction(self, rxn_id: str) -> Dict[str, Any]:
        """
        Get reaction data from simulator.

        :param rxn_id: Reaction identifier
        :return: Dictionary with reaction data (id, name, lb, ub, stoichiometry, gpr, etc.)
        """
        return self._simulator.get_reaction(rxn_id)

    def get_gene(self, gene_id: str) -> Dict[str, Any]:
        """
        Get gene data from simulator.

        :param gene_id: Gene identifier
        :return: Dictionary with gene data (id, name, reactions)
        """
        return self._simulator.get_gene(gene_id)

    def get_metabolite(self, met_id: str) -> Dict[str, Any]:
        """
        Get metabolite data from simulator.

        :param met_id: Metabolite identifier
        :return: Dictionary with metabolite data (id, name, compartment, formula)
        """
        return self._simulator.get_metabolite(met_id)

    def get_compartment(self, comp_id: str) -> Dict[str, Any]:
        """
        Get compartment data from simulator.

        :param comp_id: Compartment identifier
        :return: Dictionary with compartment data
        """
        return self._simulator.get_compartment(comp_id)

    def get_exchange_reactions(self):
        """Get exchange reactions from simulator."""
        return self._simulator.get_exchange_reactions()

    def get_gene_reactions(self):
        """Get gene-to-reaction mapping from simulator."""
        return self._simulator.get_gene_reactions()

    def get_gpr(self, rxn_id: str) -> str:
        """
        Get GPR rule string from simulator.

        :param rxn_id: Reaction identifier
        :return: GPR rule as string
        """
        return self._simulator.get_gpr(rxn_id)

    # =========================================================================
    # Simulation Methods (Delegated to Simulator)
    # =========================================================================

    def simulate(self, method="FBA", **kwargs):
        """
        Run simulation using the underlying simulator.

        :param method: Simulation method (FBA, pFBA, MOMA, etc.)
        :param kwargs: Additional simulation parameters
        :return: SimulationResult
        """
        return self._simulator.simulate(method=method, **kwargs)

    def FVA(self, **kwargs):
        """
        Flux Variability Analysis using the underlying simulator.

        :param kwargs: FVA parameters
        :return: FVA results
        """
        return self._simulator.FVA(**kwargs)

    def set_objective(self, reaction):
        """
        Set objective function in the underlying simulator.

        :param reaction: Reaction ID or dict of reaction:coefficient
        """
        return self._simulator.set_objective(reaction)

    # =========================================================================
    # Regulatory Network Management
    # =========================================================================

    def _load_regulatory_network(self, regulatory_network):
        """
        Load regulatory network from a RegulatoryModel instance.

        :param regulatory_network: RegulatoryModel instance
        """
        # Import here to avoid circular dependency
        from mewpy.germ.models.regulatory import RegulatoryModel

        if isinstance(regulatory_network, RegulatoryModel):
            # Extract regulators, targets, interactions from RegulatoryModel
            # These are already stored as dictionaries in the regulatory model
            self._regulatory_network = regulatory_network
            self._regulators = regulatory_network.regulators.copy()
            self._targets = regulatory_network.targets.copy()
            self._interactions = regulatory_network.interactions.copy()
            self._reconcile_ids()
            self._warn_if_disconnected()
        else:
            raise TypeError(f"Expected RegulatoryModel, got {type(regulatory_network)}")

    def _reconcile_ids(self, variables=None):
        """
        Rename TRN variables whose ids differ from the simulator ids only by the backend prefix.

        Rule: a TRN id that is not a simulator gene, reaction or metabolite id is renamed to
        ``simulator.resolve(id, kind)`` for the first kind, in the order gene, reaction, metabolite,
        that resolves it to a different id of that kind (the simulator owns its id prefixes, e.g.
        reframed adds or strips G_/R_/M_). Ambiguous ids are not renamed. A rename is skipped when
        the new id is already used by another TRN variable. Symbols inside the regulatory expressions
        are renamed with the variables, so targets match GPR genes and regulators match
        reaction/metabolite ids. The applied renames are available in ``id_map``
        (original TRN id -> simulator id).

        It runs when a network is loaded and whenever regulators, targets or interactions are added.

        :param variables: the TRN variables to reconcile. If None, all variables of the network
        """
        simulator = self._simulator
        resolve = getattr(simulator, "resolve", None)
        if resolve is None:
            # only mewpy Simulators know their id prefixes
            return
        simulator_ids = {kind: set(getattr(simulator, name, None) or ()) for kind, name in _ID_KINDS}
        known = set.union(*simulator_ids.values())

        if variables is None:
            variables = list(self._regulators.values()) + list(self._targets.values())
            for interaction in self._interactions.values():
                variables.append(interaction.target)
                variables.extend(interaction.yield_regulators())
        candidates = {id(variable): variable for variable in variables if variable is not None}

        # distinct objects sharing an id (e.g. a self-regulated target read as a target and a regulator)
        # are renamed together, so they keep a common id
        groups = {}
        for variable in candidates.values():
            groups.setdefault(variable.id, []).append(variable)

        used = set(self._regulators) | set(self._targets) | set(groups)
        renames = {}
        for old_id, members in groups.items():
            if not isinstance(old_id, str) or old_id in known:
                continue
            new_id = None
            for kind, _ in _ID_KINDS:
                try:
                    native = resolve(old_id, kind)
                except ValueError:
                    logger.warning("Regulatory id %s is ambiguous as a simulator %s id; not renamed", old_id, kind)
                    continue
                if native != old_id and native in simulator_ids[kind]:
                    new_id = native
                    break
            # the new id may only be taken by variables renamed earlier from the same TRN id
            if new_id is None or (new_id in used and self._id_map.get(old_id) != new_id):
                continue
            for variable in members:
                renames[id(variable)] = (variable, old_id, new_id)
            used.add(new_id)

        if not renames:
            return

        # expressions holding a renamed variable: regulators know their interactions
        interactions = {}
        for variable, _, _ in renames.values():
            for interaction in getattr(variable, "_interactions", {}).values():
                interactions[id(interaction)] = interaction
        expressions = {}
        for interaction in interactions.values():
            for expression in interaction.yield_expressions():
                expressions[id(expression)] = expression

        for expression in expressions.values():
            expression_variables = expression.variables
            renamed = {name: renames[id(var)][2] for name, var in expression_variables.items() if id(var) in renames}
            if not renamed:
                continue
            for symbol in expression.symbolic.atoms(symbols_only=True):
                if symbol.name in renamed:
                    symbol.value = renamed[symbol.name]
                    symbol.name = symbol.value
            for variable in expression_variables.values():
                if id(variable) in renames:
                    variable._id = renames[id(variable)][2]
            # re-link the expression with the new keys
            expression.variables = {var.id: var for var in expression_variables.values()}

        for variable, _, new_id in renames.values():
            variable._id = new_id

        containers = [self, self._regulatory_network] if self._regulatory_network is not None else [self]
        for owner in containers:
            # keep the containers (and the wrapped RegulatoryModel) keyed by the new identifiers
            owner._regulators = {var.id: var for var in owner._regulators.values()}
            owner._targets = {var.id: var for var in owner._targets.values()}

        renamed_ids = {old_id: new_id for _, old_id, new_id in renames.values()}
        self._id_map.update(renamed_ids)
        logger.info("Renamed %d regulatory ids to match the simulator ids", len(renamed_ids))

    def _warn_if_disconnected(self):
        """Warn when no regulatory target matches a metabolic gene (the TRN cannot affect metabolism)."""
        if not self._targets:
            return
        genes = set(getattr(self._simulator, "genes", None) or ())
        if not genes.intersection(self._targets):
            logger.warning(
                "RegulatoryExtension %s: no regulatory target matches a metabolic gene of the "
                "simulator (%d targets, %d genes); the regulatory network cannot affect metabolism",
                self._identifier,
                len(self._targets),
                len(genes),
            )

    @property
    def id_map(self) -> Dict[str, str]:
        """Renames applied to TRN ids to match the simulator (original TRN id -> simulator id)."""
        return self._id_map.copy()

    # =========================================================================
    # Observer pattern (analysis methods attached to the model)
    # =========================================================================

    @property
    def simulators(self) -> List[Any]:
        """Analysis methods attached to this model."""
        return self._simulators

    def attach(self, simulator: Any):
        """
        Attach an analysis method. It is notified when the regulatory network changes.
        The observer is also attached to the wrapped RegulatoryModel, so changes to regulatory
        variables notify it too.

        :param simulator: the analysis method to be attached
        """
        if simulator not in self._simulators:
            self._simulators.append(simulator)
        network = self._regulatory_network
        if network is not None and simulator not in network.simulators:
            network.attach(simulator)

    def detach(self, simulator: Any):
        """
        Detach an analysis method.

        :param simulator: the analysis method to be detached
        """
        self._simulators.remove(simulator)
        network = self._regulatory_network
        if network is not None and simulator in network.simulators:
            network.detach(simulator)

    def notify(self):
        """Notify all attached analysis methods that the model changed."""
        for simulator in self._simulators:
            notify_observer(simulator)

    def add_regulator(self, regulator: "Regulator"):
        """
        Add a regulator to the regulatory network.

        :param regulator: Regulator instance
        """
        self._regulators[regulator.id] = regulator
        self._reconcile_ids([regulator])
        self.notify()

    def add_target(self, target: "Target"):
        """
        Add a target to the regulatory network.

        :param target: Target instance
        """
        self._targets[target.id] = target
        self._reconcile_ids([target])
        self.notify()

    def add_interaction(self, interaction: "Interaction"):
        """
        Add an interaction to the regulatory network.

        :param interaction: Interaction instance
        """
        self._interactions[interaction.id] = interaction
        self._reconcile_ids([interaction.target, *interaction.yield_regulators()])
        self.notify()

    def remove_regulator(self, regulator_id: str):
        """
        Remove a regulator from the regulatory network.

        :param regulator_id: Regulator identifier
        """
        if regulator_id in self._regulators:
            del self._regulators[regulator_id]
            self.notify()

    def remove_target(self, target_id: str):
        """
        Remove a target from the regulatory network.

        :param target_id: Target identifier
        """
        if target_id in self._targets:
            del self._targets[target_id]
            self.notify()

    def remove_interaction(self, interaction_id: str):
        """
        Remove an interaction from the regulatory network.

        :param interaction_id: Interaction identifier
        """
        if interaction_id in self._interactions:
            del self._interactions[interaction_id]
            self.notify()

    def get_regulator(self, regulator_id: str) -> "Regulator":
        """
        Get a regulator by ID.

        :param regulator_id: Regulator identifier
        :return: Regulator instance
        """
        return self._regulators.get(regulator_id)

    def get_target(self, target_id: str) -> "Target":
        """
        Get a target by ID.

        :param target_id: Target identifier
        :return: Target instance
        """
        return self._targets.get(target_id)

    def get_interaction(self, interaction_id: str) -> "Interaction":
        """
        Get an interaction by ID.

        :param interaction_id: Interaction identifier
        :return: Interaction instance
        """
        return self._interactions.get(interaction_id)

    def get(self, identifier: str, default=None) -> Union["Regulator", "Target", "Interaction", Any]:
        """
        Get a regulator, target or interaction by ID (same lookup interface as GERM models).

        :param identifier: Regulator, target or interaction identifier
        :param default: Value returned when the identifier is not found
        :return: the regulatory variable or the default value
        """
        for container in (self._regulators, self._targets, self._interactions):
            if identifier in container:
                return container[identifier]
        return default

    def yield_regulators(self) -> Generator[Tuple[str, "Regulator"], None, None]:
        """
        Yield all regulators.

        :return: Generator of (regulator_id, Regulator) tuples
        """
        for reg_id, regulator in self._regulators.items():
            yield reg_id, regulator

    def yield_targets(self) -> Generator[Tuple[str, "Target"], None, None]:
        """
        Yield all targets.

        :return: Generator of (target_id, Target) tuples
        """
        for tgt_id, target in self._targets.items():
            yield tgt_id, target

    def yield_interactions(self) -> Generator[Tuple[str, "Interaction"], None, None]:
        """
        Yield all interactions.

        :return: Generator of (interaction_id, Interaction) tuples

        **Note:** Changed in v1.0 to return tuples for API consistency with
        yield_regulators() and yield_targets(). If you need just the interaction
        objects, use: `for _, interaction in model.yield_interactions(): ...`
        """
        for int_id, interaction in self._interactions.items():
            yield int_id, interaction

    def yield_reactions(self) -> Generator[str, None, None]:
        """
        Yield reactions from simulator (for legacy compatibility).

        :return: Generator of reaction IDs
        """
        # For simulator-based models, just yield reaction IDs
        # Legacy GERM models yield reaction objects, but we yield IDs for simplicity
        for rxn_id in self.reactions:
            yield rxn_id

    def yield_metabolites(self) -> Generator[str, None, None]:
        """
        Yield metabolites from simulator (for legacy compatibility).

        :return: Generator of metabolite IDs
        """
        # For simulator-based models, just yield metabolite IDs
        for met_id in self.metabolites:
            yield met_id

    def yield_genes(self) -> Generator[str, None, None]:
        """
        Yield genes from simulator (for legacy compatibility).

        :return: Generator of gene IDs
        """
        # For simulator-based models, just yield gene IDs
        for gene_id in self.genes:
            yield gene_id

    @property
    def regulators(self) -> Dict[str, "Regulator"]:
        """Dictionary of regulators."""
        return self._regulators.copy()

    @property
    def targets(self) -> Dict[str, "Target"]:
        """Dictionary of targets."""
        return self._targets.copy()

    @property
    def interactions(self) -> Dict[str, "Interaction"]:
        """Dictionary of interactions."""
        return self._interactions.copy()

    # =========================================================================
    # GPR Parsing with Caching
    # =========================================================================

    def get_parsed_gpr(self, rxn_id: str) -> "Expression":
        """
        Get parsed GPR expression for a reaction (with caching).

        :param rxn_id: Reaction identifier
        :return: Parsed Expression object
        """
        if rxn_id not in self._gpr_cache:
            gpr_str = self.get_gpr(rxn_id)
            if gpr_str and gpr_str.strip():
                try:
                    self._gpr_cache[rxn_id] = parse_expression(gpr_str)
                except Exception:
                    # If parsing fails, create None expression
                    from mewpy.germ.algebra import Expression

                    self._gpr_cache[rxn_id] = Expression()
            else:
                # Empty GPR
                from mewpy.germ.algebra import Expression

                self._gpr_cache[rxn_id] = Expression()

        return self._gpr_cache[rxn_id]

    def clear_gpr_cache(self):
        """Clear the GPR expression cache."""
        self._gpr_cache.clear()

    # =========================================================================
    # Utility Methods
    # =========================================================================

    def is_metabolic(self) -> bool:
        """Check if model has metabolic component (always True)."""
        return True

    def is_regulatory(self) -> bool:
        """Check if model has regulatory component."""
        return len(self._interactions) > 0

    def has_regulatory_network(self) -> bool:
        """Check if regulatory network is present."""
        return self.is_regulatory()

    # =========================================================================
    # Serialization
    # =========================================================================

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialize to dictionary.

        :return: Dictionary representation
        """
        return {
            "id": self._identifier,
            "simulator_type": type(self._simulator).__name__,
            "simulator_id": getattr(self._simulator, "id", None),
            "regulators": {reg_id: reg.to_dict() for reg_id, reg in self._regulators.items()},
            "targets": {tgt_id: tgt.to_dict() for tgt_id, tgt in self._targets.items()},
            "interactions": {int_id: inter.to_dict() for int_id, inter in self._interactions.items()},
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any], simulator: "Simulator") -> "RegulatoryExtension":
        """
        Deserialize from dictionary.

        :param data: Dictionary representation
        :param simulator: Simulator instance to wrap
        :return: RegulatoryExtension instance
        """
        from mewpy.germ.models.regulatory import RegulatoryModel

        # Rebuild the network as a single model, so interactions, regulators and targets share the same objects
        identifier = data.get("id")
        network = {
            "id": identifier or "regulatory",
            "name": identifier or "",
            "types": ["regulatory"],
            "interactions": data.get("interactions", {}),
            "regulators": data.get("regulators", {}),
            "targets": data.get("targets", {}),
        }
        regulatory_model = RegulatoryModel.from_dict(network, variables=True)

        return cls(simulator, regulatory_model, identifier)

    def save(self, filepath: str):
        """
        Save regulatory network to JSON file.

        Note: This saves only the regulatory network. The metabolic model
        should be saved separately using its native format (SBML, etc.).

        :param filepath: Path to output JSON file
        """
        data = self.to_dict()
        with open(filepath, "w") as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load(cls, filepath: str, simulator: "Simulator") -> "RegulatoryExtension":
        """
        Load regulatory network from JSON file.

        :param filepath: Path to JSON file
        :param simulator: Simulator instance to wrap
        :return: RegulatoryExtension instance
        """
        with open(filepath, "r") as f:
            data = json.load(f)
        return cls.from_dict(data, simulator)

    # =========================================================================
    # String Representation
    # =========================================================================

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"RegulatoryExtension(id='{self._identifier}', "
            f"simulator={type(self._simulator).__name__}, "
            f"regulators={len(self._regulators)}, "
            f"targets={len(self._targets)}, "
            f"interactions={len(self._interactions)})"
        )

    def __str__(self) -> str:
        """Human-readable string."""
        return (
            f"RegulatoryExtension '{self._identifier}'\n"
            f"  Metabolic: {len(self.reactions)} reactions, "
            f"{len(self.genes)} genes, {len(self.metabolites)} metabolites\n"
            f"  Regulatory: {len(self._regulators)} regulators, "
            f"{len(self._targets)} targets, {len(self._interactions)} interactions"
        )
