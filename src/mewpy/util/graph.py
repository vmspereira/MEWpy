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
Graph methods for metabolic models.
Generates a networkx graph based on metabolic networks

Author: Vitor Pereira
##############################################################################
"""
import math
import networkx as nx
import numpy as np
from mewpy.simulation import get_simulator
from .constants import COFACTORS

METABOLITE = 'METABOLITE'
REACTION = 'REACTION    '
REV = 'REV'
IRREV = 'IRREV'




def create_metabolic_graph(model, directed=True, carbon=True, reactions=None, remove=[], edges_labels=False, biomass=False, metabolites=False):
    """ Creates a metabolic graph

    :param model: A model or a model containter
    :param (bool) directed: Defines if the graph to be directed or undirected. Defaults to True.
    :param (bool) carbon: Only included edges for metabolites with a carbon atom. Defaults to True.
    :param (list) reactions: List of reactions to be included in the graph. Defaults to None, in which\
        all reactions are included.
    :param list remove: list os metabolites not to be included. May be used to remove cofactores such as ATP/ADP, \
        NAD(P)(H), and acetyl-CoA/CoA.
    :param (bool) edges_labels: Adds a reversabily label to edges. Defaults to False.
    :returns: A networkx graph of the metabolic network.
    """


    container = get_simulator(model)

    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    if not reactions:
        reactions = container.reactions

    reactions = list(set(reactions) - set(remove))

    for r in reactions:
        G.add_node(r, label=r, node_class=REACTION, node_id=r)

    for r in reactions:
        the_metabolites = container.get_reaction_metabolites(r)
        for m in the_metabolites:
            if m in remove or container.get_metabolite(m)['formula'] in COFACTORS.values():
                continue
            if carbon and 'C' not in container.metabolite_elements(m).keys():
                continue
            if m not in G.nodes:
                G.add_node(m, label=m, node_class=METABOLITE, node_id=m)
            # evaluating if the metabolite has been defined as a reactant or product
            if the_metabolites[m] < 0:
                (tail, head) = (m, r)
            elif the_metabolites[m] > 0:
                (tail, head) = (r, m)

            # adding an arc between a metabolite and a reactions
            G.add_edge(tail, head)
            label = IRREV
            lb, _ = container.get_reaction_bounds(r)

            if lb < 0:
                G.add_edge(head, tail)
                label = REV

            if edges_labels:
                G[tail][head]['label'] = label

            G[tail][head]['reversible'] = lb < 0

    if not metabolites:
        met_nodes = [x for x, v in dict(G.nodes(data="node_class")).items() if v == METABOLITE]
        for m in met_nodes:
            in_ = G.in_edges(m, data=True)
            out_ = G.out_edges(m, data=True)
            for s, _, r1 in in_:
                for _, t, r2 in out_:
                    try:
                        rev = r1['reversible'] and r2['reversible']
                    except:
                        rev = False
                    G.add_edge(s, t)
            G.remove_node(m)

    return G


def filter_by_degree(G, max_degree, inplace=True):
    s = list(sorted(G.degree, key=lambda item: item[1], reverse=True))
    stop = False
    while not stop:
        # find the metabolite with highest degree
        print(s[:5])
        position = 0
        found = False
        k = None
        v = None
        while not found and position < len(s):
            k, v = s[position]
            if G.nodes[k]['node_class'] == METABOLITE:
                found = True
            else:
                position += 1
        if k and v > max_degree:
            G.remove_node(k)
            print('removed ', k)
            s = list(sorted(G.degree, key=lambda item: item[1], reverse=True))
        else:
            stop = True
    return G


def shortest_distance(model, reaction, reactions=None, remove=[]):
    """ Returns the unweighted shortest path distance from a list of reactions to a reaction.
    Distances are the number of required reactions. If there is no pathway between the reactions the distance is inf·

    :param model: A model or a Simulator instance.
    :param str reaction: target reaction.
    :param list reactions: List os source reactions. Defaults to None, in which case all model reactions are considered.
    :param list remove: List os metabolites not to be included. May be used to remove path that include \
        cofactores such as ATP/ADP, NAD(P)(H), and acetyl-CoA/CoA.
    :returns: A dictionary of distances.
    """
    container = get_simulator(model)
    
    rxns = reactions if reactions else container.reactions
    if reaction not in rxns:
        rxns.append(reaction)

    G = create_metabolic_graph(container, reactions=rxns, remove=remove)
    sp = dict(nx.single_target_shortest_path_length(G, reaction))

    distances = {}
    for rxn in rxns:
        if rxn in sp:
            distances[rxn] = sp[rxn] // 2
        else:
            distances[rxn] = np.inf
    return distances


def probabilistic_reaction_targets(model, product, targets, factor=10):
    """Builds a new target list reflecting the shortest path distances from all original
    as a probability,ie, reactions closer to the product are repeated more often in the new target list.
    Moreover, reactions from which there is no path (pathway or cofactors usage) to the product are removed.

    :param model: A model or a Simulator instance.
    :param str product: Product to be optimized.
    :param str targets: EA target reactions.
    :param int factor: Maximum number of repetitions, also the distance after which all reactions are\
        considered with equal probability. Defaults to 10.
    :returns: A probabilistic target list.
    """
    distances = shortest_distance(model, product, targets)
    prob_targets = []
    for t in targets:
        if distances[t] == np.inf or distances[t] == 0:
            continue
        else:
            coef = math.ceil(1 / distances[t] * factor)
        x = [t] * coef
        prob_targets.extend(x)
    return prob_targets


def probabilistic_gene_targets(model, product, targets, factor=10):
    """Builds a new target list reflecting the shortest path distances from all original
    as a probability,ie, genes on GPRs of reactions closer to the product are repeated more
    often in the new target list.

    :param model: A model or a Simulator instance.
    :param str product: Product to be optimized.
    :param str targets: EA target genes.
    :param int factor: Maximum number of repetitions. Defaults to 10.
    :returns: A probabilistic target list.
    """
    
    container = get_simulator(model)
    
    # Reaction targets
    if not targets:
        genes = container.genes
    else:
        genes = targets

    rxns = container.get_reactions_for_genes(genes)
    rxn_distances = shortest_distance(model, product, rxns)

    # genes distances are the maximum of all reaction
    # distances that they catalyse.

    prob_targets = []

    for gene in genes:
        rxs = container.get_reactions_for_genes([gene])
        dd = [d for r, d in rxn_distances.items() if r in rxs]
        d = max(dd)
        if d == np.inf:
            coef = 1
        elif d == 0:
            continue
        else:
            coef = math.ceil(1 / d * factor)
        x = [gene] * coef
        prob_targets.extend(x)

    return prob_targets


def probabilistic_protein_targets(model, product, targets, factor=10):
    """Builds a new target list reflecting the shortest path distances from all original
    as a probability,ie, proteins used in reactions closer to the product are repeated
    more often in the new target list.

    :param model: A model or a Simulator instance.
    :param str product: Product to be optimized.
    :param str targets: EA target genes.
    :param int factor: Maximum number of repetitions. Defaults to 10.
    :returns: A probabilistic target list.
    """
    raise NotImplementedError
