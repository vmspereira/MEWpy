"""
RegulatoryExtension: Decorator that adds regulatory network capabilities to Simulators.

This module provides the RegulatoryExtension class, which wraps external simulators
(COBRApy, reframed) and adds regulatory network functionality without duplicating
metabolic data. All metabolic operations are delegated to the wrapped simulator.
"""

import json
from typing import TYPE_CHECKING, Any, Dict, Generator, Optional, Tuple, Union

if TYPE_CHECKING:
    from mewpy.germ.algebra import Expression
    from mewpy.germ.variables import Interaction, Regulator, Target
    from mewpy.simulation.simulation import Simulator

# Runtime imports
from mewpy.germ.algebra import parse_expression


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
        >>> reg_network = RegulatoryModel.from_file('regulatory.json')
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
            self._regulators = regulatory_network.regulators.copy()
            self._targets = regulatory_network.targets.copy()
            self._interactions = regulatory_network.interactions.copy()
        else:
            raise TypeError(f"Expected RegulatoryModel, got {type(regulatory_network)}")

    def add_regulator(self, regulator: "Regulator"):
        """
        Add a regulator to the regulatory network.

        :param regulator: Regulator instance
        """
        self._regulators[regulator.id] = regulator

    def add_target(self, target: "Target"):
        """
        Add a target to the regulatory network.

        :param target: Target instance
        """
        self._targets[target.id] = target

    def add_interaction(self, interaction: "Interaction"):
        """
        Add an interaction to the regulatory network.

        :param interaction: Interaction instance
        """
        self._interactions[interaction.id] = interaction

    def remove_regulator(self, regulator_id: str):
        """
        Remove a regulator from the regulatory network.

        :param regulator_id: Regulator identifier
        """
        if regulator_id in self._regulators:
            del self._regulators[regulator_id]

    def remove_target(self, target_id: str):
        """
        Remove a target from the regulatory network.

        :param target_id: Target identifier
        """
        if target_id in self._targets:
            del self._targets[target_id]

    def remove_interaction(self, interaction_id: str):
        """
        Remove an interaction from the regulatory network.

        :param interaction_id: Interaction identifier
        """
        if interaction_id in self._interactions:
            del self._interactions[interaction_id]

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
        from mewpy.germ.variables import Interaction, Regulator, Target

        extension = cls(simulator, identifier=data.get("id"))

        # Load regulators
        for reg_id, reg_data in data.get("regulators", {}).items():
            regulator = Regulator.from_dict(reg_data)
            extension.add_regulator(regulator)

        # Load targets
        for tgt_id, tgt_data in data.get("targets", {}).items():
            target = Target.from_dict(tgt_data)
            extension.add_target(target)

        # Load interactions
        for int_id, int_data in data.get("interactions", {}).items():
            interaction = Interaction.from_dict(int_data)
            extension.add_interaction(interaction)

        return extension

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
