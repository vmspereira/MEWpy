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
COBRA utility module

Authors: Vitor Pereira
##############################################################################
"""
import logging
from copy import copy, deepcopy
from math import inf
from typing import TYPE_CHECKING, Union

from tqdm import tqdm

from mewpy.simulation import Simulator, get_simulator
from mewpy.util.constants import ModelConstants
from mewpy.util.parsing import Boolean, build_tree, isozymes

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from cobra import Model
    from reframed.core.cbmodel import CBModel


def convert_gpr_to_dnf(model) -> None:
    """
    Convert all existing GPR (Gene-Protein-Reaction) associations to DNF (Disjunctive Normal Form).

    DNF is a standardized form where the GPR rule is expressed as an OR of ANDs,
    e.g., (geneA and geneB) or (geneC and geneD)

    :param model: A COBRApy or REFRAMED Model or an instance of Simulator
    """
    sim = get_simulator(model)
    for rxn_id in tqdm(sim.reactions, desc="Converting GPRs to DNF"):
        rxn = sim.get_reaction(rxn_id)
        if not rxn.gpr:
            continue
        try:
            tree = build_tree(rxn.gpr, Boolean)
            gpr_dnf = tree.to_dnf().to_infix()
            # Update the reaction's GPR
            rxn.gpr = gpr_dnf
        except Exception as e:
            # If conversion fails, keep original GPR
            import warnings

            warnings.warn(f"Failed to convert GPR for reaction {rxn_id}: {e}")


def convert_to_irreversible(model: Union[Simulator, "Model", "CBModel"]):
    """Split reversible reactions into two irreversible reactions
    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    :param model: A COBRApy or REFRAMED Model or an instance of
        mewpy.simulation.simulation.Simulator
    :return: A new irreversible model simulator and a reverse mapping.
    :rtype: (Simulator, dict)

    .. note::
        This function always returns a new model; the input model is not modified.
    """

    sim = get_simulator(deepcopy(model))

    objective = sim.objective.copy()
    irrev_map = dict()
    for r_id in tqdm(sim.reactions, "Converting to irreversible"):
        lb, _ = sim.get_reaction_bounds(r_id)
        if lb < 0:
            rxn = sim.get_reaction(r_id)
            rev_rxn_id = r_id + "_REV"

            rev_rxn = dict()
            rev_rxn["name"] = rxn.name + " reverse"
            rev_rxn["lb"] = 0
            rev_rxn["ub"] = -rxn.lb
            rev_rxn["gpr"] = rxn.gpr
            sth = {k: v * -1 for k, v in rxn.stoichiometry.items()}
            rev_rxn["stoichiometry"] = sth
            rev_rxn["reversible"] = False
            rev_rxn["annotations"] = copy(rxn.annotations)

            sim.add_reaction(rev_rxn_id, **rev_rxn)
            sim.set_reaction_bounds(r_id, 0, rxn.ub, False)

            irrev_map[r_id] = rev_rxn_id

            if r_id in objective:
                objective[rev_rxn_id] = -objective[r_id]

    sim.objective = objective
    return sim, irrev_map


def split_isozymes(model: Union[Simulator, "Model", "CBModel"]):
    """Splits reactions with isozymes into separated reactions

    :param model: A COBRApy or REFRAMED Model or an instance of
        mewpy.simulation.simulation.Simulator
    :return: A new simulator with split isozyme reactions and a mapping from original to splitted reactions
    :rtype: (Simulator, dict)

    .. note::
        This function always returns a new model; the input model is not modified.
    """

    sim = get_simulator(deepcopy(model))

    objective = sim.objective
    mapping = dict()
    newobjective = {}

    for r_id in tqdm(sim.reactions, "Splitting isozymes"):
        rxn = sim.get_reaction(r_id)
        gpr = rxn.gpr

        if gpr is not None and len(gpr.strip()) > 0:
            proteins = isozymes(gpr)
            mapping[r_id] = []
            for i, protein in enumerate(proteins):
                r_id_new = "{}_No{}".format(r_id, i + 1)
                mapping[r_id].append(r_id_new)

                rxn_new = dict()
                rxn_new["name"] = "{} No{}".format(rxn.name, i + 1)
                rxn_new["lb"] = rxn.lb
                rxn_new["ub"] = rxn.ub
                rxn_new["gpr"] = protein
                rxn_new["stoichiometry"] = rxn.stoichiometry.copy()
                rxn_new["annotations"] = copy(rxn.annotations)

                sim.add_reaction(r_id_new, **rxn_new)
            sim.remove_reaction(r_id)

    # set the objective
    for k, v in objective.items():
        if k in mapping.keys():
            for r in mapping[k]:
                newobjective[r] = v
        else:
            newobjective[k] = v

    sim.objective = newobjective

    return sim, mapping


def __enzime_constraints(
    model: Union[Simulator, "Model", "CBModel"],
    prot_mw=None,
    enz_kcats=None,
    c_compartment: str = "c",
    inline: bool = False,
):
    """Auxiliary method to add enzyme constraints to a model

    :param model: A model or simulator
    :type model: A COBRApy or REFRAMED Model or an instance of
        mewpy.simulation.simulation.Simulator
    :param data: Protein MW and Kcats
    :type data: None
    :param c_compartment: The compartment where gene/proteins pseudo species are to be added.
        Defaults to 'c'
    :type c_compartment: str, optional
    :param (boolean) inline: apply the modifications to the same of generate a new model.
        Default generates a new model.
    :type inline: bool, optional
    :return: a new enzyme constrained model
    :rtype: Simulator
    """

    if inline:
        sim = get_simulator(model)
    else:
        sim = deepcopy(get_simulator(model))

    objective = sim.objective

    if prot_mw is None:
        prot_mw = dict()
        for gene in sim.genes:
            prot_mw[gene] = {"protein": gene[len(sim._g_prefix) :], "mw": 1}

    if enz_kcats is None:
        enz_kcats = dict()
        for gene in sim.genes:
            enz_kcats[gene] = dict()
            rxns = sim.get_gene(gene).reactions
            for rxn in rxns:
                enz_kcats[gene][rxn] = {"protein": gene[len(sim._g_prefix) :], "kcat": 1}

    # Add protein pool and species
    common_protein_pool_id = sim._m_prefix + "prot_pool_c"
    pool_reaction = sim._r_prefix + "prot_pool_exchange"

    sim.add_metabolite(common_protein_pool_id, name="prot_pool [cytoplasm]", compartment=c_compartment)

    sim.add_reaction(
        pool_reaction,
        name="protein pool exchange",
        stoichiometry={common_protein_pool_id: 1},
        lb=0,
        ub=inf,
        reversible=False,
        reaction_type="EX",
    )

    # Add gene/protein species and draw protein pseudo-reactions
    # MW in kDa, [kDa = g/mmol]
    gene_meta = dict()
    skipped_gene = []
    for gene in tqdm(sim.genes, "Adding gene species"):
        if gene not in prot_mw.keys():
            skipped_gene.append(gene)
            continue
        mw = prot_mw[gene]
        m_prot_id = f"prot_{mw['protein']}_{c_compartment}"
        m_name = f"prot_{mw['protein']} {c_compartment}"
        sim.add_metabolite(m_prot_id, name=m_name, compartment=c_compartment)

        gene_meta[gene] = m_prot_id

        r_prot_id = f"draw_prot_{mw['protein']}"
        sim.add_reaction(
            r_prot_id,
            name=r_prot_id,
            stoichiometry={common_protein_pool_id: -1 * mw["mw"], m_prot_id: 1},
            lb=0,
            ub=inf,
            reversible=False,
            gpr=gene,
        )

    if skipped_gene:
        logger.info(
            f"{len(skipped_gene)} gene species not added (missing protein MW data). " f"First few: {skipped_gene[:5]}"
        )

    # Add enzymes to reactions stoichiometry.
    # 1/Kcats in per hour. Considering kcats in per second.
    for rxn_id in tqdm(sim.reactions, "Adding proteins usage to reactions"):
        rxn = sim.get_reaction(rxn_id)
        if rxn.gpr:
            s = rxn.stoichiometry
            genes = build_tree(rxn.gpr, Boolean).get_operands()
            for g in genes:
                if g in gene_meta:
                    # Get kcat from enz_kcats dictionary (gene -> reaction -> kcat mapping)
                    kcat = ModelConstants.DEFAULT_KCAT  # Default value
                    if g in enz_kcats and rxn_id in enz_kcats[g]:
                        kcat_data = enz_kcats[g][rxn_id]
                        if isinstance(kcat_data.get("kcat"), (int, float)):
                            kcat = kcat_data["kcat"]

                    s[gene_meta[g]] = -1 / kcat
            sim.update_stoichiometry(rxn_id, s)
    sim.objective = objective
    return sim


def add_enzyme_constraints(
    model: Union[Simulator, "Model", "CBModel"],
    prot_mw=None,
    enz_kcats=None,
    c_compartment: str = "c",
):
    """Adds enzyme constraints to a model.

    This function applies a series of transformations to prepare a model for enzyme constraints:
    1. Converts reversible reactions to irreversible
    2. Splits reactions with isozymes
    3. Adds enzyme constraints

    :param model: A model or simulator
    :type model: A COBRApy or REFRAMED Model or an instance of
        mewpy.simulation.simulation.Simulator
    :param prot_mw: Dictionary mapping gene IDs to protein molecular weight data
    :type prot_mw: dict, optional
    :param enz_kcats: Dictionary mapping gene IDs to kcat values per reaction
    :type enz_kcats: dict, optional
    :param c_compartment: The compartment where gene/proteins pseudo species are to be added.
        Defaults to 'c'
    :type c_compartment: str, optional
    :return: A new enzyme constrained model
    :rtype: Simulator

    .. note::
        This function always returns a new model; the input model is not modified.
    """
    sim, _ = convert_to_irreversible(model)
    sim, _ = split_isozymes(sim)
    sim = __enzime_constraints(sim, prot_mw=prot_mw, enz_kcats=enz_kcats, c_compartment=c_compartment, inline=True)
    return sim
