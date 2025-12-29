# Copyright (C) 2019-2023 Centre of Biological Engineering,
#     University of Minho, Portugal
# Vitor Pereira 2019-

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
Compartmentalized community model.
Can build community models from models loaded from any toolbox
for which there is a Simulator implementation.

Author: Vitor Pereira
##############################################################################
"""
from copy import deepcopy
from typing import TYPE_CHECKING, Dict, List, Union
from warnings import warn

from numpy import inf
from tqdm import tqdm

from mewpy.simulation import get_simulator
from mewpy.util import AttrDict
from mewpy.util.parsing import Boolean, build_tree

if TYPE_CHECKING:
    from cobra.core import Model
    from reframed.core.cbmodel import CBModel

    from mewpy.simulation import Simulator


class CommunityModel:

    EXT_COMP = "e"
    GROWTH_ID = "community_growth"

    def __init__(
        self,
        models: List[Union["Simulator", "Model", "CBModel"]],
        abundances: List[float] = None,
        merge_biomasses: bool = True,
        copy_models: bool = False,
        add_compartments=True,
        balance_exchange=False,
        flavor: str = "reframed",
        verbose: bool = True,
    ):
        """Community Model.

        :param models: A list of metabolic models.
        :param abundances: A list of relative abundances for each model.
            Default None.
        :param merge_biomasses: If a biomass equation is to be build requiring
            each organism to grow in acordance to a relative abundance.
            Default True.
            If no abundance list is provided, all organism will have equal abundance.
        :param add_compartments: If each organism external compartment is to be added
            to the community model. Default True.
        :param balance_exchange: **DEPRECATED - May violate mass conservation.**
            If True, modifies stoichiometric coefficients of exchange metabolites
            based on organism abundances. This approach is mathematically problematic
            as it violates conservation of mass (e.g., 1 mol consumed produces only
            0.3 mol if abundance=0.3). Default False.
            Note: Abundance scaling is already handled through the merged biomass equation.
            This parameter will be removed in a future version.
        :param bool copy_models: if the models are to be copied, default True.
        :param str flavor: use 'cobrapy' or 'reframed. Default 'reframed'.
        :param bool verbose: show progress bar during model merging (default True).
            Set to False for batch processing or when building many communities.
        """
        self.organisms = AttrDict()
        self.model_ids = list({model.id for model in models})
        self._verbose = verbose
        if len(self.model_ids) != len(set(self.model_ids)):
            raise ValueError("Each model must have a different ID.")
        self.flavor = flavor

        self.organisms_biomass = None
        self.organisms_biomass_metabolite = None
        self.biomass = None

        self.reaction_map = None
        self.metabolite_map = None
        self.gene_map = None
        self.ext_mets = None

        self._reverse_map = None
        if abundances and any(e <= 0 for e in abundances):
            raise ValueError("All abundances need to be positive")

        self._merge_biomasses = True if abundances is not None else merge_biomasses
        self._add_compartments = add_compartments
        self._balance_exchange = balance_exchange

        # Warn if balance_exchange is enabled (deprecated feature with mass balance issues)
        if balance_exchange:
            warn(
                "balance_exchange=True is deprecated and may violate conservation of mass. "
                "Stoichiometric coefficients of exchange reactions are being modified based on "
                "organism abundances, which can lead to mass imbalance (e.g., 1 mol consumed "
                "producing only 0.3 mol if abundance=0.3). "
                "Abundance scaling is already handled through the merged biomass equation. "
                "This parameter will be removed in a future version. "
                "Set balance_exchange=False to suppress this warning.",
                DeprecationWarning,
                stacklevel=2,
            )

        if len(self.model_ids) < len(models):
            warn("Model ids are not unique, repeated models will be discarded.")

        for model in models:
            m = get_simulator(model)
            if not m.objective:
                raise ValueError(f"Model {m.id} has no objective")
            self.organisms[m.id] = deepcopy(m) if copy_models else m

        if self._merge_biomasses:
            if abundances and len(abundances) == len(self.organisms):
                self.organisms_abundance = dict(zip(self.organisms.keys(), abundances))
            else:
                self.organisms_abundance = {org_id: 1 for org_id in self.organisms.keys()}

        self._comm_model = None

    def __repr__(self):
        """Rich representation showing community model details."""
        lines = []
        lines.append("=" * 60)
        lines.append("Community Model")
        lines.append("=" * 60)

        # Number of organisms
        try:
            org_count = len(self.organisms)
            lines.append(f"{'Organisms:':<20} {org_count}")

            # List organisms (up to 5, then truncate)
            if org_count > 0:
                org_ids = list(self.organisms.keys())
                if org_count <= 5:
                    for org_id in org_ids:
                        lines.append(f"{'  -':<20} {org_id}")
                else:
                    for org_id in org_ids[:5]:
                        lines.append(f"{'  -':<20} {org_id}")
                    lines.append(f"{'  ...':<20} and {org_count - 5} more")
        except:
            pass

        # Flavor
        try:
            if self.flavor:
                lines.append(f"{'Flavor:':<20} {self.flavor}")
        except:
            pass

        # Abundances
        try:
            if hasattr(self, "organisms_abundance") and self.organisms_abundance:
                # Check if abundances are uniform
                abundances = list(self.organisms_abundance.values())
                if len(set(abundances)) == 1:
                    lines.append(f"{'Abundances:':<20} Uniform ({abundances[0]})")
                else:
                    lines.append(f"{'Abundances:':<20} Variable")
                    # Show first few
                    items = list(self.organisms_abundance.items())[:3]
                    for org_id, abundance in items:
                        lines.append(f"{'  ' + org_id + ':':<20} {abundance:.4g}")
                    if len(self.organisms_abundance) > 3:
                        lines.append(f"{'  ...':<20}")
        except:
            pass

        # Configuration
        try:
            if self._merge_biomasses:
                lines.append(f"{'Merged biomass:':<20} Yes")
            if self._add_compartments:
                lines.append(f"{'Add compartments:':<20} Yes")
        except:
            pass

        # Community model built status
        try:
            if self._comm_model is not None:
                lines.append(f"{'Status:':<20} Built")
                # Get community model stats
                if hasattr(self._comm_model, "reactions"):
                    rxn_count = len(self._comm_model.reactions)
                    lines.append(f"{'Community reactions:':<20} {rxn_count}")
                if hasattr(self._comm_model, "metabolites"):
                    met_count = len(self._comm_model.metabolites)
                    lines.append(f"{'Community metabolites:':<20} {met_count}")
            else:
                lines.append(f"{'Status:':<20} Not built (call merge() or build())")
        except:
            pass

        # Biomass reaction
        try:
            if self.biomass:
                lines.append(f"{'Community biomass:':<20} {self.biomass}")
        except:
            pass

        lines.append("=" * 60)
        return "\n".join(lines)

    def _repr_html_(self):
        """Pandas-like HTML representation for Jupyter notebooks."""
        from mewpy.util.html_repr import render_html_table

        rows = []

        org_count = len(self.organisms)
        rows.append(("Organisms", str(org_count)))

        if org_count > 0:
            org_ids = list(self.organisms.keys())
            if org_count <= 5:
                for org_id in org_ids:
                    rows.append((f"  {org_id}", ""))
            else:
                for org_id in org_ids[:5]:
                    rows.append((f"  {org_id}", ""))
                rows.append(("  ...", f"and {org_count - 5} more"))

        if self.flavor:
            rows.append(("Flavor", self.flavor))

        if hasattr(self, "organisms_abundance") and self.organisms_abundance:
            abundances = list(self.organisms_abundance.values())
            if len(set(abundances)) == 1:
                rows.append(("Abundances", f"Uniform ({abundances[0]})"))
            else:
                rows.append(("Abundances", "Variable"))
                items = list(self.organisms_abundance.items())[:3]
                for org_id, abundance in items:
                    rows.append((f"  {org_id}", f"{abundance:.4g}"))
                if len(self.organisms_abundance) > 3:
                    rows.append(("  ...", ""))

        if self._merge_biomasses:
            rows.append(("Merged biomass", "Yes"))
        if self._add_compartments:
            rows.append(("Add compartments", "Yes"))

        if self._comm_model is not None:
            rows.append(("Status", "Built"))
            if hasattr(self._comm_model, "reactions"):
                rxn_count = len(self._comm_model.reactions)
                rows.append(("Community reactions", str(rxn_count)))
            if hasattr(self._comm_model, "metabolites"):
                met_count = len(self._comm_model.metabolites)
                rows.append(("Community metabolites", str(met_count)))
        else:
            rows.append(("Status", "Not built (call merge() or build())"))

        if self.biomass:
            rows.append(("Community biomass", self.biomass))

        return render_html_table("Community Model", rows)

    def init_model(self):
        sid = " ".join(sorted(self.model_ids))
        if self.flavor == "reframed":
            from reframed.core.cbmodel import CBModel

            model = CBModel(sid)
        else:
            from cobra.core.model import Model

            model = Model(sid)
        self._comm_model = get_simulator(model)

    def clear(self):
        self.organisms_biomass = None
        self.organisms_biomass_metabolite = None
        self.biomass = None
        self.reaction_map = None
        self.metabolite_map = None
        self.gene_map = None
        self.ext_mets = None
        self._reverse_map = None
        self._comm_model = None

    @property
    def add_compartments(self):
        return self._add_compartments

    @add_compartments.setter
    def add_compartments(self, value: bool):
        if self._add_compartments == value:
            pass
        else:
            self._add_compartments = value
            self.clear()

    @property
    def balance_exchanges(self):
        return self._balance_exchange

    @balance_exchanges.setter
    def balance_exchanges(self, value: bool):
        if value == self._balance_exchange:
            return
        if value:
            warn(
                "balance_exchange=True is deprecated and may violate conservation of mass. "
                "This parameter will be removed in a future version.",
                DeprecationWarning,
                stacklevel=2,
            )
        self._balance_exchange = value
        if value:
            self._update_exchanges()
        else:
            self._update_exchanges({k: 1 for k in self.model_ids})

    @property
    def merge_biomasses(self):
        return self._merge_biomasses

    @merge_biomasses.setter
    def merge_biomasses(self, value: bool):
        if self._merge_biomasses == value:
            pass
        else:
            self._merge_biomasses = value
            self.clear()

    @property
    def reverse_map(self):
        if self._reverse_map is not None:
            return self._reverse_map
        else:
            self._reverse_map = dict()
            self._reverse_map.update({v: k for k, v in self.reaction_map.items()})
            self._reverse_map.update({v: k for k, v in self.gene_map.items()})
            return self._reverse_map

    def set_abundance(self, abundances: Dict[str, float], rebuild=False):
        if not self._merge_biomasses:
            raise ValueError("The community model has no merged biomass equation")
        # Validate organism IDs
        invalid_orgs = set(abundances.keys()) - set(self.organisms.keys())
        if invalid_orgs:
            raise ValueError(
                f"Unknown organism IDs: {invalid_orgs}. " f"Valid organisms are: {set(self.organisms.keys())}"
            )
        if any([x < 0 for x in abundances.values()]):
            raise ValueError("All abundance value need to be non negative.")
        if sum(list(abundances.values())) == 0:
            raise ValueError("At least one organism needs to have a positive abundance.")
        # update the biomass equation
        self.organisms_abundance.update(abundances)
        if rebuild:
            self.clear()
            self._merge_models()
        else:
            comm_growth = CommunityModel.GROWTH_ID
            biomass_stoichiometry = {
                met: -self.organisms_abundance[org_id]
                for org_id, met in self.organisms_biomass_metabolite.items()
                if self.organisms_abundance[org_id] > 0
            }
            self._comm_model.add_reaction(
                comm_growth,
                name="Community growth rate",
                stoichiometry=biomass_stoichiometry,
                lb=0,
                ub=inf,
                reaction_type="SINK",
            )
            self._comm_model.objective = comm_growth
            self._comm_model.solver = None

        if self._balance_exchange:
            self._update_exchanges()

    def _update_exchanges(self, abundances: dict = None):
        """
        Update exchange reaction stoichiometry based on organism abundances.

        WARNING: This method modifies stoichiometric coefficients which violates
        conservation of mass. For example, if abundance=0.3, a transport reaction
        M_org <-> M_ext with stoichiometry {M_org: -1, M_ext: 1} becomes
        {M_org: -1, M_ext: 0.3}, meaning 1 mol consumed produces only 0.3 mol.

        This feature is DEPRECATED and will be removed in a future version.
        Abundance scaling should be handled through flux constraints or is already
        addressed by the merged biomass equation.

        :param abundances: Optional dict of organism abundances to use instead of
                          self.organisms_abundance
        """
        if self.merged_model and self._merge_biomasses and self._balance_exchange:
            exchange = self.merged_model.get_exchange_reactions()
            m_r = self.merged_model.metabolite_reaction_lookup()
            for met in self.ext_mets:
                rxns = m_r[met]
                for rx, st in rxns.items():
                    if rx in exchange:
                        continue
                    org = self.reverse_map[rx][0]
                    if abundances:
                        ab = abundances[org]
                    else:
                        ab = self.organisms_abundance[org]
                    rxn = self.merged_model.get_reaction(rx)
                    stch = rxn.stoichiometry
                    new_stch = stch.copy()
                    new_stch[met] = ab if st > 0 else -ab
                    self.merged_model.update_stoichiometry(rx, new_stch)

    def get_community_model(self):
        """Returns a Simulator for the merged model"""
        return self.merged_model

    def size(self):
        return len(self.organisms)

    def get_organisms_biomass(self) -> Dict[str, str]:
        return self.organisms_biomass

    @property
    def merged_model(self):
        """Returns a community model (COBRApy or REFRAMED)"""
        if self._comm_model is None:
            self._merge_models()
        return self._comm_model

    def _merge_models(self):
        """Merges the models with optimizations for large communities."""

        self.init_model()

        old_ext_comps = []
        self.ext_mets = []
        self._reverse_map = None

        # Pre-calculate dictionary sizes for better memory efficiency
        total_reactions = sum(len(model.reactions) for model in self.organisms.values())
        total_metabolites = sum(len(model.metabolites) for model in self.organisms.values())
        total_genes = sum(len(model.genes) for model in self.organisms.values())

        # Pre-allocate dictionaries with estimated sizes (reduces reallocation)
        self.organisms_biomass = {}
        self.reaction_map = dict() if total_reactions < 1000 else {}
        self.metabolite_map = dict() if total_metabolites < 1000 else {}
        self.gene_map = dict() if total_genes < 1000 else {}

        if self._merge_biomasses:
            self.organisms_biomass_metabolite = {}

        # default IDs
        ext_comp_id = CommunityModel.EXT_COMP
        comm_growth = CommunityModel.GROWTH_ID

        # create external compartment
        self._comm_model.add_compartment(ext_comp_id, "extracellular environment", external=True)

        # community biomass
        if not self._merge_biomasses:
            biomass_id = "community_biomass"
            self._comm_model.add_metabolite(biomass_id, name="Total community biomass", compartment=ext_comp_id)

        # add each organism (with optional progress bar)
        organism_iter = tqdm(self.organisms.items(), "Organism") if self._verbose else self.organisms.items()
        for org_id, model in organism_iter:

            # Cache prefix information to avoid repeated string operations
            g_prefix_match = model._g_prefix == self._comm_model._g_prefix
            m_prefix_match = model._m_prefix == self._comm_model._m_prefix
            r_prefix_match = model._r_prefix == self._comm_model._r_prefix

            g_prefix_len = len(model._g_prefix) if not g_prefix_match else 0
            m_prefix_len = len(model._m_prefix) if not m_prefix_match else 0
            r_prefix_len = len(model._r_prefix) if not r_prefix_match else 0

            def rename(old_id):
                return f"{old_id}_{org_id}"

            def r_gene(old_id, organism=True):
                if g_prefix_match:
                    _id = old_id
                else:
                    _id = self._comm_model._g_prefix + old_id[g_prefix_len:]
                return rename(_id) if organism else _id

            def r_met(old_id, organism=True):
                if m_prefix_match:
                    _id = old_id
                else:
                    _id = self._comm_model._m_prefix + old_id[m_prefix_len:]
                return rename(_id) if organism else _id

            def r_rxn(old_id, organism=True):
                if r_prefix_match:
                    _id = old_id
                else:
                    _id = self._comm_model._r_prefix + old_id[r_prefix_len:]
                return rename(_id) if organism else _id

            # add internal compartments
            for c_id in model.compartments:
                comp = model.get_compartment(c_id)
                if comp.external:
                    old_ext_comps.append(c_id)
                    if not self._add_compartments:
                        continue
                self._comm_model.add_compartment(rename(c_id), name=f"{comp.name} ({org_id})")

            # add metabolites
            for m_id in model.metabolites:
                met = model.get_metabolite(m_id)
                if met.compartment not in old_ext_comps or self._add_compartments:
                    new_mid = r_met(m_id)
                    self._comm_model.add_metabolite(
                        new_mid,
                        formula=met.formula,
                        name=met.name,
                        compartment=rename(met.compartment),
                    )
                    self.metabolite_map[(org_id, m_id)] = new_mid

                if met.compartment in old_ext_comps and r_met(m_id, False) not in self._comm_model.metabolites:
                    new_mid = r_met(m_id, False)
                    self._comm_model.add_metabolite(
                        new_mid,
                        formula=met.formula,
                        name=met.name,
                        compartment=ext_comp_id,
                    )
                    self.ext_mets.append(new_mid)

            # add genes
            for g_id in model.genes:
                new_id = r_gene(g_id)
                self.gene_map[(org_id, g_id)] = new_id
                if self.flavor == "reframed":
                    gene = model.get_gene(g_id)
                    self._comm_model.add_gene(new_id, gene.name)

            # add reactions
            ex_rxns = model.get_exchange_reactions()

            for r_id in model.reactions:
                rxn = model.get_reaction(r_id)
                new_id = r_rxn(r_id)

                if r_id in ex_rxns:
                    mets = list(rxn.stoichiometry.keys())

                    if self._add_compartments and r_met(mets[0], False) in self.ext_mets:
                        new_stoichiometry = {
                            r_met(mets[0]): -1,
                            r_met(mets[0], False): 1,
                        }
                        self._comm_model.add_reaction(
                            new_id,
                            name=rxn.name,
                            stoichiometry=new_stoichiometry,
                            lb=-inf,
                            ub=inf,
                            reaction_type="TRP",
                        )
                        self.reaction_map[(org_id, r_id)] = new_id

                    elif len(mets) == 1 and r_met(mets[0]) in self._comm_model.metabolites:
                        # some models (e.g. AGORA models) have sink reactions (for biomass)
                        new_stoichiometry = {r_met(mets[0]): -1}
                        self._comm_model.add_reaction(
                            new_id,
                            name=rxn.name,
                            stoichiometry=new_stoichiometry,
                            lb=0,
                            ub=inf,
                            reaction_type="SINK",
                        )
                        self.reaction_map[(org_id, r_id)] = new_id

                else:
                    if self._add_compartments:
                        new_stoichiometry = {r_met(m_id): coeff for m_id, coeff in rxn.stoichiometry.items()}
                    else:
                        new_stoichiometry = {
                            r_met(m_id, False) if r_met(m_id, False) in self.ext_mets else r_met(m_id): coeff
                            for m_id, coeff in rxn.stoichiometry.items()
                        }
                    # assumes that the models' objective is the biomass
                    if r_id in [x for x, v in model.objective.items() if v > 0]:
                        if self._merge_biomasses:
                            met_id = r_met("Biomass")
                            self._comm_model.add_metabolite(
                                met_id,
                                name=f"Biomass {org_id}",
                                compartment=ext_comp_id,
                            )

                            new_stoichiometry[met_id] = 1
                            self.organisms_biomass_metabolite[org_id] = met_id

                            # add biomass sink reaction
                            self._comm_model.add_reaction(
                                r_rxn("Sink_biomass"),
                                name=f"Sink Biomass {org_id}",
                                stoichiometry={met_id: -1},
                                lb=0,
                                ub=inf,
                                reaction_type="SINK",
                            )

                        else:
                            new_stoichiometry[biomass_id] = 1

                        self.organisms_biomass[org_id] = new_id

                    if rxn.gpr:
                        t = build_tree(rxn.gpr, Boolean)
                        ren = {x: r_gene(x) for x in t.get_operands()}
                        new_gpr = t.replace(ren).to_infix()
                    else:
                        new_gpr = rxn.gpr

                    self._comm_model.add_reaction(
                        new_id,
                        name=rxn.name,
                        stoichiometry=new_stoichiometry,
                        lb=rxn.lb,
                        ub=rxn.ub,
                        gpr=new_gpr,
                        annotations=rxn.annotations,
                    )

                    self.reaction_map[(org_id, r_id)] = new_id

        # Add exchange reactions
        for m_id in self.ext_mets:
            m = m_id[len(self._comm_model._m_prefix) :] if m_id.startswith(self._comm_model._m_prefix) else m_id
            r_id = f"{self._comm_model._r_prefix}EX_{m}"
            self._comm_model.add_reaction(
                r_id,
                name=r_id,
                stoichiometry={m_id: -1},
                lb=-inf,
                ub=inf,
                reaction_type="EX",
            )

        if self._merge_biomasses:
            # if the biomasses are to be merged add
            # a new product to each organism biomass
            biomass_stoichiometry = {
                met: -1 * self.organisms_abundance[org_id] for org_id, met in self.organisms_biomass_metabolite.items()
            }
        else:
            biomass_stoichiometry = {biomass_id: -1}

        self._comm_model.add_reaction(
            comm_growth,
            name="Community growth rate",
            stoichiometry=biomass_stoichiometry,
            lb=0,
            ub=inf,
            reaction_type="SINK",
        )

        if self._balance_exchange:
            self._update_exchanges()

        self._comm_model.objective = comm_growth
        self._comm_model.biomass_reaction = comm_growth
        self.biomass = comm_growth
        setattr(self._comm_model, "organisms_biomass", self.organisms_biomass)
        setattr(self._comm_model, "community", self)
        return self._comm_model

    def copy(self, copy_models=False, flavor=None):
        models = [m.model for m in self.organisms.values()]
        f = flavor if flavor is not None else self.flavor
        return CommunityModel(models, copy_models=copy_models, flavor=f)
