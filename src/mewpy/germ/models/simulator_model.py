"""
SimulatorBasedMetabolicModel: GERM-compatible wrapper for external simulators.

This module provides a MetabolicModel-compatible interface that wraps external
simulators (COBRApy, reframed) to provide the same API as GERM MetabolicModel
but sourcing ALL data directly from the simulator interface (not the underlying model).
"""

import copy
from collections import defaultdict
from typing import TYPE_CHECKING, Any, Dict, Generator, List, Set, Tuple, Union

from mewpy.germ.models.serialization import serialize
from mewpy.util.history import recorder
from mewpy.util.utilities import generator

from .model import Model

if TYPE_CHECKING:
    from mewpy.germ.algebra import Expression
    from mewpy.germ.variables import Gene, Metabolite, Reaction
    from mewpy.simulation.simulation import Simulator

# Import these for runtime use
from mewpy.germ.algebra import Expression, parse_expression


def _rebuild_simulator_model(simulator, identifier, name):
    """Pickle helper: a SimulatorBasedMetabolicModel is fully defined by its simulator."""
    return SimulatorBasedMetabolicModel(simulator, identifier=identifier, name=name)


class SimulatorBasedMetabolicModel(
    Model, model_type="simulator_metabolic", register=True, constructor=False, checker=False
):
    """
    A GERM-compatible MetabolicModel that wraps external simulators.

    This class provides the same interface as MetabolicModel but sources
    ALL data directly from the simulator interface (never accessing the underlying model directly).

    Key principles:
    - All data access goes through simulator methods (get_reaction, get_metabolite, etc.)
    - GERM variables are created on-demand from simulator data
    - Changes are routed through simulator interface when possible
    - Maintains complete compatibility with existing MEWpy ecosystem
    """

    def __init__(self, simulator: "Simulator", identifier: Any = None, **kwargs):
        """
        Initialize a SimulatorBasedMetabolicModel from an external simulator.

        :param simulator: External simulator instance (COBRApy, reframed, etc.)
        :param identifier: Model identifier (defaults to simulator's model ID)
        """
        self._simulator = simulator

        # Use simulator's model ID if no identifier provided
        if identifier is None:
            identifier = getattr(simulator, "id", "simulator_model")

        # Cache for GERM variables created from simulator data
        self._germ_genes_cache = {}
        self._germ_metabolites_cache = {}
        self._germ_reactions_cache = {}
        self._germ_compartments_cache = {}
        # bounds last synchronized with the simulator, per cached GERM reaction
        self._synced_bounds = {}

        # Initialize parent
        super().__init__(identifier, **kwargs)

    # -----------------------------------------------------------------------------
    # Simulator access
    # -----------------------------------------------------------------------------

    @property
    def simulator(self) -> "Simulator":
        """Access to the underlying simulator."""
        return self._simulator

    def _invalidate_caches(self):
        """Invalidate all GERM variable caches."""
        self._germ_genes_cache.clear()
        self._germ_metabolites_cache.clear()
        self._germ_reactions_cache.clear()
        self._germ_compartments_cache.clear()
        self._synced_bounds.clear()

    # -----------------------------------------------------------------------------
    # Copy and pickle: the model is fully defined by its simulator
    # -----------------------------------------------------------------------------
    def __reduce__(self):
        return _rebuild_simulator_model, (self._simulator, self._id, self._name)

    def __copy__(self):
        return SimulatorBasedMetabolicModel(self._simulator, identifier=self._id, name=self._name)

    def __deepcopy__(self, memo):
        simulator = copy.deepcopy(self._simulator, memo)
        return SimulatorBasedMetabolicModel(simulator, identifier=self._id, name=self._name)

    def notify(self):
        """
        It pushes the bounds changed on cached GERM reactions to the simulator and notifies the observers.
        """
        for rxn_id, reaction in self._germ_reactions_cache.items():
            bounds = tuple(reaction.bounds)
            if self._synced_bounds.get(rxn_id) != bounds:
                self._simulator.set_reaction_bounds(rxn_id, bounds[0], bounds[1])
                self._synced_bounds[rxn_id] = bounds

        super().notify()

    # -----------------------------------------------------------------------------
    # GERM Variable Creation from Simulator Data
    # -----------------------------------------------------------------------------

    def _create_germ_metabolite(self, met_id: str) -> "Metabolite":
        """Create a GERM Metabolite from simulator data."""
        if met_id in self._germ_metabolites_cache:
            return self._germ_metabolites_cache[met_id]

        from mewpy.germ.variables import Metabolite

        # Get metabolite data from simulator
        met_data = self._simulator.get_metabolite(met_id)

        # Create GERM metabolite
        germ_met = Metabolite(
            identifier=met_id,
            name=met_data.get("name", met_id),
            compartment=met_data.get("compartment"),
            formula=met_data.get("formula"),
        )

        self._germ_metabolites_cache[met_id] = germ_met
        return germ_met

    def _create_germ_gene(self, gene_id: str) -> "Gene":
        """Create a GERM Gene from simulator data."""
        if gene_id in self._germ_genes_cache:
            return self._germ_genes_cache[gene_id]

        from mewpy.germ.variables import Gene

        # Get gene data from simulator
        gene_data = self._simulator.get_gene(gene_id)

        # Create GERM gene
        germ_gene = Gene(identifier=gene_id, name=gene_data.get("name", gene_id))

        self._germ_genes_cache[gene_id] = germ_gene
        return germ_gene

    def _create_germ_reaction(self, rxn_id: str) -> "Reaction":
        """Create a GERM Reaction from simulator data."""
        if rxn_id in self._germ_reactions_cache:
            return self._germ_reactions_cache[rxn_id]

        from mewpy.germ.algebra import parse_expression
        from mewpy.germ.variables import Reaction

        # Get reaction data from simulator
        rxn_data = self._simulator.get_reaction(rxn_id)

        # Build stoichiometry with GERM metabolites
        stoichiometry = {}
        for met_id, coeff in rxn_data.get("stoichiometry", {}).items():
            germ_met = self._create_germ_metabolite(met_id)
            stoichiometry[germ_met] = coeff

        # Handle GPR
        gpr_str = rxn_data.get("gpr")
        if gpr_str and gpr_str.strip():
            try:
                # Parse GPR expression and create GERM genes as needed
                parsed_gpr = parse_expression(gpr_str)

                # Create gene variables dictionary
                genes = {}
                for gene_id in {symbol.name for symbol in parsed_gpr.atoms(symbols_only=True)}:
                    germ_gene = self._create_germ_gene(gene_id)
                    genes[gene_id] = germ_gene

                # Create Expression with symbolic and variables
                gpr = Expression(symbolic=parsed_gpr, variables=genes)

            except Exception as e:
                # If GPR parsing fails, create empty expression
                print(f"Warning: Failed to parse GPR '{gpr_str}': {e}")
                gpr = Expression()
        else:
            # No GPR - use empty expression
            gpr = Expression()

        # Create GERM reaction
        germ_rxn = Reaction(
            identifier=rxn_id,
            name=rxn_data.get("name", rxn_id),
            stoichiometry=stoichiometry,
            bounds=(rxn_data.get("lb", -1000), rxn_data.get("ub", 1000)),
            gpr=gpr,
        )

        # attached to this model, so bounds changes (bounds, ko, undo) are pushed to the simulator by notify
        germ_rxn.model = self
        self._synced_bounds[rxn_id] = tuple(germ_rxn.bounds)

        self._germ_reactions_cache[rxn_id] = germ_rxn
        return germ_rxn

    def _extract_genes_from_gpr(self, gpr_str: str) -> Set[str]:
        """Extract gene identifiers from GPR string (the symbols of the parsed GPR expression)."""
        return {symbol.name for symbol in parse_expression(gpr_str).atoms(symbols_only=True)}

    # -----------------------------------------------------------------------------
    # MetabolicModel-compatible interface
    # -----------------------------------------------------------------------------

    @serialize("types", None)
    @property
    def types(self):
        """Returns the types of the model."""
        _types = {SimulatorBasedMetabolicModel.model_type}
        _types.update(super().types)
        return _types

    # -----------------------------------------------------------------------------
    # Core properties - route to simulator
    # -----------------------------------------------------------------------------

    @serialize("genes", "genes", None)
    @property
    def genes(self) -> Dict[str, "Gene"]:
        """
        Returns a dictionary with the genes from the simulator.
        Creates GERM Gene objects on-demand.
        """
        genes = {}
        for gene_id in self._simulator.genes:
            genes[gene_id] = self._create_germ_gene(gene_id)
        return genes

    @serialize("metabolites", "metabolites", None)
    @property
    def metabolites(self) -> Dict[str, "Metabolite"]:
        """
        Returns a dictionary with the metabolites from the simulator.
        Creates GERM Metabolite objects on-demand.
        """
        metabolites = {}
        for met_id in self._simulator.metabolites:
            metabolites[met_id] = self._create_germ_metabolite(met_id)
        return metabolites

    @serialize("reactions", "reactions", None)
    @property
    def reactions(self) -> Dict[str, "Reaction"]:
        """
        Returns a dictionary with the reactions from the simulator.
        Creates GERM Reaction objects on-demand.
        """
        reactions = {}
        for rxn_id in self._simulator.reactions:
            reactions[rxn_id] = self._create_germ_reaction(rxn_id)
        return reactions

    @serialize("objective", "objective", None)
    @property
    def objective(self) -> Dict["Reaction", Union[float, int]]:
        """
        Returns the objective function from the simulator.
        """
        objective = {}
        sim_objective = getattr(self._simulator, "objective", {})

        for rxn_id, coeff in sim_objective.items():
            if rxn_id in self._simulator.reactions:
                germ_rxn = self._create_germ_reaction(rxn_id)
                objective[germ_rxn] = coeff

        return objective

    @property
    def compartments(self) -> Dict[str, str]:
        """
        Returns a dictionary with the compartments from the simulator.
        """
        if hasattr(self._simulator, "compartments"):
            compartments = {}
            sim_compartments = self._simulator.compartments

            if isinstance(sim_compartments, dict):
                compartments.update(sim_compartments)
            else:
                # Handle case where compartments is a list
                for comp_id in sim_compartments:
                    comp_data = self._simulator.get_compartment(comp_id)
                    compartments[comp_id] = comp_data.get("name", comp_id)

            return compartments
        else:
            # Infer compartments from metabolites
            compartments = {}
            for met_id in self._simulator.metabolites:
                met_data = self._simulator.get_metabolite(met_id)
                comp_id = met_data.get("compartment")
                if comp_id and comp_id not in compartments:
                    compartments[comp_id] = comp_id
            return compartments

    # -----------------------------------------------------------------------------
    # Setters - modify simulator directly
    # -----------------------------------------------------------------------------

    @compartments.setter
    def compartments(self, value: Dict[str, str]):
        """Set compartments - not supported for simulator models."""
        raise NotImplementedError("Compartments of a SimulatorBasedMetabolicModel cannot be set; edit the simulator")

    @genes.setter
    def genes(self, value: Dict[str, "Gene"]):
        """Set genes - not supported for simulator models."""
        raise NotImplementedError("Genes of a SimulatorBasedMetabolicModel cannot be set; edit the simulator")

    @metabolites.setter
    def metabolites(self, value: Dict[str, "Metabolite"]):
        """Set metabolites - not supported for simulator models."""
        raise NotImplementedError("Metabolites of a SimulatorBasedMetabolicModel cannot be set; edit the simulator")

    @objective.setter
    @recorder
    def objective(self, value: Dict["Reaction", Union[float, int]]):
        """Set objective function on the simulator."""
        if not value:
            value = {}

        if isinstance(value, str):
            value = {self.get(value): 1}
        elif hasattr(value, "types"):
            value = {value: 1}
        elif isinstance(value, dict):
            value = {self.get(var, var): val for var, val in value.items()}
        else:
            raise ValueError(f"{value} is not a valid objective")

        # Convert to simulator format
        linear_obj = {}
        for var, coef in value.items():
            if hasattr(var, "id"):
                linear_obj[var.id] = coef
            else:
                linear_obj[str(var)] = coef

        # Set objective on simulator
        if hasattr(self._simulator, "set_objective"):
            self._simulator.set_objective(linear=linear_obj, minimize=False)

        self._invalidate_caches()

    @reactions.setter
    def reactions(self, value: Dict[str, "Reaction"]):
        """Set reactions - not supported for simulator models."""
        raise NotImplementedError("Reactions of a SimulatorBasedMetabolicModel cannot be set; edit the simulator")

    # -----------------------------------------------------------------------------
    # Dynamic attributes (same logic as MetabolicModel)
    # -----------------------------------------------------------------------------

    @property
    def external_compartment(self) -> Union[str, None]:
        """
        Returns the external compartment of the model.
        """
        if not self.compartments:
            return None

        # Try to use simulator's method if available
        if hasattr(self._simulator, "get_exchange_reactions"):
            exchange_reactions = self._simulator.get_exchange_reactions()

            boundary_compartments = defaultdict(int)

            for rxn_id in exchange_reactions:
                rxn_data = self._simulator.get_reaction(rxn_id)
                stoichiometry = rxn_data.get("stoichiometry", {})

                for met_id in stoichiometry:
                    met_data = self._simulator.get_metabolite(met_id)
                    compartment = met_data.get("compartment")
                    if compartment:
                        boundary_compartments[compartment] += 1

            if boundary_compartments:
                return max(boundary_compartments, key=boundary_compartments.get)

        return None

    def _get_boundaries(self) -> Tuple[Dict[str, "Reaction"], Dict[str, "Reaction"], Dict[str, "Reaction"]]:
        """Returns the boundary reactions of the model."""
        external_compartment = self.external_compartment

        if external_compartment is None:
            return {}, {}, {}

        exchanges = {}
        sinks = {}
        demands = {}

        # Use simulator's exchange reactions if available
        if hasattr(self._simulator, "get_exchange_reactions"):
            exchange_rxn_ids = self._simulator.get_exchange_reactions()

            for rxn_id in exchange_rxn_ids:
                germ_rxn = self._create_germ_reaction(rxn_id)
                exchanges[rxn_id] = germ_rxn
        else:
            # Fallback: identify boundary reactions manually
            for rxn_id in self._simulator.reactions:
                rxn_data = self._simulator.get_reaction(rxn_id)
                stoichiometry = rxn_data.get("stoichiometry", {})

                # Check if it's a boundary reaction (single metabolite)
                if len(stoichiometry) == 1:
                    met_id = list(stoichiometry.keys())[0]
                    met_data = self._simulator.get_metabolite(met_id)
                    compartment = met_data.get("compartment")

                    germ_rxn = self._create_germ_reaction(rxn_id)

                    if compartment == external_compartment:
                        exchanges[rxn_id] = germ_rxn
                    else:
                        # Determine if it's sink or demand based on bounds
                        lb = rxn_data.get("lb", -1000)
                        ub = rxn_data.get("ub", 1000)

                        if lb < 0 and ub > 0:
                            sinks[rxn_id] = germ_rxn
                        else:
                            demands[rxn_id] = germ_rxn

        return exchanges, sinks, demands

    @property
    def demands(self) -> Dict[str, "Reaction"]:
        """Returns the demand reactions of the model."""
        _, _, demands = self._get_boundaries()
        return demands

    @property
    def exchanges(self) -> Dict[str, "Reaction"]:
        """Returns the exchange reactions of the model."""
        exchanges, _, _ = self._get_boundaries()
        return exchanges

    @property
    def sinks(self) -> Dict[str, "Reaction"]:
        """Returns the sink reactions of the model."""
        _, sinks, _ = self._get_boundaries()
        return sinks

    # -----------------------------------------------------------------------------
    # Generators
    # -----------------------------------------------------------------------------

    def yield_compartments(self) -> Generator[str, None, None]:
        """Yields the compartments of the model."""
        return generator(self.compartments)

    def yield_demands(self) -> Generator["Reaction", None, None]:
        """Yields the demand reactions of the model."""
        return generator(self.demands)

    def yield_exchanges(self) -> Generator["Reaction", None, None]:
        """Yields the exchange reactions of the model."""
        return generator(self.exchanges)

    def yield_genes(self) -> Generator["Gene", None, None]:
        """Yields the genes of the model."""
        return generator(self.genes)

    def yield_gprs(self) -> Generator["Expression", None, None]:
        """Yields the GPRs of the model."""
        for rxn in self.yield_reactions():
            if rxn.gpr:
                yield rxn.gpr

    def yield_metabolites(self) -> Generator["Metabolite", None, None]:
        """Yields the metabolites of the model."""
        return generator(self.metabolites)

    def yield_reactions(self) -> Generator["Reaction", None, None]:
        """Yields the reactions of the model."""
        return generator(self.reactions)

    def yield_sinks(self) -> Generator["Reaction", None, None]:
        """Yields the sink reactions of the model."""
        return generator(self.sinks)

    # -----------------------------------------------------------------------------
    # Operations/Manipulations
    # -----------------------------------------------------------------------------

    def get(self, identifier: Any, default=None) -> Union["Gene", "Metabolite", "Reaction"]:
        """
        Returns the object associated with the identifier.
        """
        # Check metabolites first
        if identifier in self._simulator.metabolites:
            return self._create_germ_metabolite(identifier)

        # Check reactions
        if identifier in self._simulator.reactions:
            return self._create_germ_reaction(identifier)

        # Check genes
        if identifier in self._simulator.genes:
            return self._create_germ_gene(identifier)

        # Fall back to parent
        return super().get(identifier=identifier, default=default)

    def add(self, *variables, comprehensive: bool = True, history: bool = True, **kwargs):
        """
        Adding variables is not supported: the metabolic data lives in the simulator.
        Add reactions, metabolites or genes to the simulator (or its underlying model) instead.
        """
        raise NotImplementedError("Variables cannot be added to a SimulatorBasedMetabolicModel; edit the simulator")

    def remove(self, *variables, remove_orphans: bool = False, history: bool = True, **kwargs):
        """
        Removing variables is not supported: the metabolic data lives in the simulator.
        Remove reactions, metabolites or genes from the simulator (or its underlying model) instead.
        """
        raise NotImplementedError("Variables cannot be removed from a SimulatorBasedMetabolicModel; edit the simulator")

    def update(self, compartments=None, objective=None, variables=None, **kwargs):
        """
        Update the model with relevant information.
        """
        if objective is not None:
            self.objective = objective

        # Other updates are typically not supported for simulator models
        self._invalidate_caches()
        super().update(**kwargs)

    # -----------------------------------------------------------------------------
    # Simulation methods - delegate to simulator
    # -----------------------------------------------------------------------------

    def simulate(self, method="FBA", **kwargs):
        """Run simulation using the underlying simulator."""
        return self._simulator.simulate(method=method, **kwargs)

    def FVA(self, **kwargs):
        """Run FVA using the underlying simulator."""
        return self._simulator.FVA(**kwargs)
