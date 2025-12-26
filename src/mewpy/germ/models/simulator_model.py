"""
SimulatorBasedMetabolicModel: GERM-compatible wrapper for external simulators.

This module provides a MetabolicModel-compatible interface that wraps external
simulators (COBRApy, reframed) to provide the same API as GERM MetabolicModel
but sourcing ALL data directly from the simulator interface (not the underlying model).
"""

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
                for gene_id in self._extract_genes_from_gpr(gpr_str):
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

        self._germ_reactions_cache[rxn_id] = germ_rxn
        return germ_rxn

    def _extract_genes_from_gpr(self, gpr_str: str) -> Set[str]:
        """Extract gene identifiers from GPR string."""
        # Simple extraction - in practice might need more sophisticated parsing
        import re

        # Find all identifiers that look like genes (alphanumeric + underscore)
        genes = set(re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", gpr_str))
        # Filter out logical operators
        logical_ops = {"and", "or", "not", "AND", "OR", "NOT", "(", ")"}
        return genes - logical_ops

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

    @serialize("genes", "genes", "_genes")
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

    @serialize("metabolites", "metabolites", "_metabolites")
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

    @serialize("reactions", "reactions", "_reactions")
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

    @serialize("objective", "objective", "_objective")
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
    @recorder
    def compartments(self, value: Dict[str, str]):
        """Set compartments - not directly supported for simulator models."""
        # Simulator compartments are typically read-only
        # Could potentially add compartments to simulator if it supports it
        self._invalidate_caches()

    @genes.setter
    @recorder
    def genes(self, value: Dict[str, "Gene"]):
        """Set genes - not directly supported for simulator models."""
        # Genes are typically derived from reactions in simulators
        self._invalidate_caches()

    @metabolites.setter
    @recorder
    def metabolites(self, value: Dict[str, "Metabolite"]):
        """Set metabolites - not directly supported for simulator models."""
        # Metabolites are typically derived from reactions in simulators
        self._invalidate_caches()

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
    @recorder
    def reactions(self, value: Dict[str, "Reaction"]):
        """Set reactions - not directly supported for simulator models."""
        # Adding/removing reactions from simulators is complex
        self._invalidate_caches()

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

    def add(self, *variables, comprehensive: bool = True, history: bool = True):
        """
        Add variables to the simulator model.
        Note: Adding to simulator models is limited compared to GERM models.
        """
        # For simulator models, adding variables is typically not supported
        # or requires complex simulator-specific operations
        self._invalidate_caches()
        return super().add(*variables, comprehensive=comprehensive, history=history)

    def remove(self, *variables, remove_orphans: bool = False, history: bool = True):
        """
        Remove variables from the simulator model.
        Note: Removing from simulator models is limited compared to GERM models.
        """
        # For simulator models, removing variables is typically not supported
        # or requires complex simulator-specific operations
        self._invalidate_caches()
        return super().remove(*variables, remove_orphans=remove_orphans, history=history)

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
