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
Defines an interface for the simulators, objects that wrap metabolic phenotipe
simulation toolboxes.

Author: Vitor Pereira
##############################################################################
"""
import logging
import math
from abc import ABC, abstractmethod
from collections import OrderedDict
from copy import deepcopy
from enum import Enum
from typing import List, Union

from joblib import Parallel, delayed
from tqdm import tqdm

from ..util.parsing import evaluate_expression_tree
from ..util.process import cpu_count

logger = logging.getLogger(__name__)


class SimulationMethod(Enum):
    FBA = "FBA"
    pFBA = "pFBA"
    MOMA = "MOMA"
    lMOMA = "lMOMA"
    ROOM = "ROOM"
    NONE = "NONE"

    def __eq__(self, other):
        """Overrides equal to enable string name comparison.
        Allows to seamlessly use:
            SimulationMethod.FBA = SimulationMethod.FBA
            SimulationMethod.FBA = 'FBA'
        without requiring an additional level of comparison (SimulationMethod.FBA.name = 'FBA')
        """
        if isinstance(other, SimulationMethod):
            return super().__eq__(other)
        elif isinstance(other, str):
            return self.name == other
        else:
            return False

    def __hash__(self):
        return hash(self.name)

    def __str__(self) -> str:
        return self.name


class SStatus(Enum):
    """Enumeration of possible solution status."""

    OPTIMAL = "Optimal"
    UNKNOWN = "Unknown"
    SUBOPTIMAL = "Suboptimal"
    UNBOUNDED = "Unbounded"
    INFEASIBLE = "Infeasible"
    INF_OR_UNB = "Infeasible or Unbounded"

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name


class SimulationInterface(ABC):

    @abstractmethod
    def simulate(**kwargs):
        """Abstract method to run a simulation.

        :returns: A SimulationResult.

        """
        raise NotImplementedError


class ModelContainer(ABC):
    """Interface for Model container.
    Provides an abstraction from models implementations.
    """

    # Backend id prefixes (REFRAMED keeps the SBML R_/M_/G_ prefixes; COBRApy and GERM drop them).
    # prefix(kind) is the public form; these attributes are kept as aliases.
    _r_prefix = ""
    _m_prefix = ""
    _g_prefix = ""

    # -----------------------------------------------------------------------------
    # Identifier resolution: callers may pass ids with or without the SBML prefix
    # -----------------------------------------------------------------------------
    _ID_KINDS = {"reaction": "reactions", "metabolite": "metabolites", "gene": "genes"}
    _SBML_PREFIXES = {"reaction": "R_", "metabolite": "M_", "gene": "G_"}

    def prefix(self, kind: str) -> str:
        """The backend prefix of an id kind ('reaction', 'metabolite' or 'gene'); '' when ids are unprefixed."""
        self._check_kind(kind)
        return {"reaction": self._r_prefix, "metabolite": self._m_prefix, "gene": self._g_prefix}[kind]

    def to_native(self, base_id: str, kind: str) -> str:
        """Adds the backend prefix to a base (unprefixed) id, e.g. to name a new reaction."""
        return self.prefix(kind) + base_id

    def to_canonical(self, native_id: str, kind: str) -> str:
        """Strips the backend prefix from a backend-native id."""
        prefix = self.prefix(kind)
        return native_id[len(prefix) :] if prefix and native_id.startswith(prefix) else native_id

    def resolve(self, identifier, kind: str = "reaction"):
        """Returns the backend-native id for an id given with or without its SBML prefix.

        An exact match wins; otherwise the prefixed form, then the form with the SBML prefix stripped.
        An id that matches neither is returned unchanged, so lookups fail as before.

        :raises ValueError: if there is no exact match and both candidate forms exist.
        """
        ids = self._native_ids(kind)
        if not isinstance(identifier, str) or identifier in ids:
            return identifier
        sbml_prefix = self._SBML_PREFIXES[kind]
        candidates = {self.to_native(identifier, kind)}
        if identifier.startswith(sbml_prefix):
            candidates.add(identifier[len(sbml_prefix) :])
        found = sorted(c for c in candidates if c in ids)
        if len(found) > 1:
            raise ValueError(f"Ambiguous {kind} id '{identifier}': matches {found}")
        return found[0] if found else identifier

    def resolve_many(self, ids, kind: str = "reaction") -> list:
        """Resolves each id of an iterable (see resolve)."""
        return [self.resolve(i, kind) for i in ids]

    def resolve_keys(self, mapping, kind: str = "reaction"):
        """Returns a copy of a dict (e.g. constraints) with its keys resolved (see resolve); None stays None."""
        if mapping is None:
            return None
        resolved = OrderedDict() if isinstance(mapping, OrderedDict) else dict()
        for key, value in mapping.items():
            resolved[self.resolve(key, kind)] = value
        return resolved

    def _resolve_mapping(self, mapping, kind: str = "reaction"):
        """Like resolve_keys, but returns the mapping itself when all its keys are already native."""
        if not mapping or all(self.resolve(key, kind) == key for key in mapping):
            return mapping
        return self.resolve_keys(mapping, kind)

    def _resolve_objective(self, objective):
        """Resolves the reaction ids of an objective given as an id, a list of ids or a {id: coefficient} dict."""
        if isinstance(objective, dict):
            return self.resolve_keys(objective, "reaction")
        if isinstance(objective, str):
            return self.resolve(objective, "reaction")
        if isinstance(objective, (list, tuple)):
            return self.resolve_many(objective, "reaction")
        return objective

    def _check_kind(self, kind):
        if kind not in self._ID_KINDS:
            raise ValueError(f"Unknown id kind '{kind}'; use one of {list(self._ID_KINDS)}")

    def _id_count(self, kind):
        try:
            return len(getattr(self.model, self._ID_KINDS[kind]))
        except (AttributeError, TypeError):
            return None

    def _native_ids(self, kind):
        """Cached set of the model's ids of a kind, rebuilt when the model's count of that kind changes."""
        self._check_kind(kind)
        cache = self.__dict__.setdefault("_id_cache", {})
        count = self._id_count(kind)
        entry = cache.get(kind)
        if entry is None or count is None or entry[0] != count:
            ids = frozenset(getattr(self, self._ID_KINDS[kind]))
            entry = (count, ids)
            cache[kind] = entry
        return entry[1]

    def _invalidate_id_cache(self):
        """Drops the cached id sets; called by methods that add or remove model entities."""
        self.__dict__.pop("_id_cache", None)

    @property
    def reactions(self):
        """
        :returns: A list of reaction identifiers.

        """
        raise NotImplementedError

    @property
    def genes(self):
        """
        :returns: A list of gene identifiers.

        """
        raise NotImplementedError

    @property
    def proteins(self):
        """
        :returns: A list of proteins identifiers.

        """
        raise NotImplementedError

    @property
    def metabolites(self):
        """
        :returns: A list of metabolite identifiers.

        """
        raise NotImplementedError

    @property
    def compartments(self):
        """
        :returns: A list of compartments identifiers.

        """
        raise NotImplementedError

    @property
    def medium(self):
        raise NotImplementedError

    def get_reaction(self, r_id):
        raise NotImplementedError

    def get_gene(self, g_id):
        raise NotImplementedError

    def get_compartment(self, c_id):
        raise NotImplementedError

    def get_reaction_metabolites(self, r_id):
        rxn = self.get_reaction(r_id)
        return rxn["stoichiometry"]

    def get_substrates(self, rxn_id):
        met = self.get_reaction_metabolites(rxn_id)
        return {m_id: coeff for m_id, coeff in met.items() if coeff < 0}

    def get_products(self, rxn_id):
        met = self.get_reaction_metabolites(rxn_id)
        return {m_id: coeff for m_id, coeff in met.items() if coeff > 0}

    def get_exchange_reactions(self):
        raise NotImplementedError("Subclasses must implement get_exchange_reactions()")

    def get_gene_reactions(self):
        raise NotImplementedError

    def get_gpr(self, rxn_id):
        """
        :returns: A string representation of the reaction GPR if exists None otherwise.
        """
        rxn = self.get_reaction(rxn_id)
        return rxn["gpr"]

    def summary(self):
        logger.info(f"Metabolites: {len(self.metabolites)}")
        logger.info(f"Reactions: {len(self.reactions)}")
        logger.info(f"Genes: {len(self.genes)}")

    def set_objective(self, reaction):
        raise NotImplementedError


def current_bounds(simulator, r_id, constraints=None):
    """A reaction's bounds now: its constraint in constraints (e.g. those of a simulate call), else the simulator's
    constraint or medium, else the model's bounds.

    :param simulator: a Simulator.
    :param r_id: a native reaction id.
    :param constraints: constraints {id: value or (lb, ub)} that take precedence, or None.
    :returns: a value or a (lb, ub) tuple.
    """
    simulator_constraints = getattr(simulator, "_constraints", None)
    medium = getattr(simulator, "environmental_conditions", None)
    for source in (constraints, simulator_constraints, medium):
        if source and r_id in source:
            return source[r_id]
    return simulator.get_reaction_bounds(r_id)


class Simulator(ModelContainer, SimulationInterface):
    """
    Interface for simulators

    Results without an optimum (every backend): when a simulation is not OPTIMAL or SUBOPTIMAL, its
    SimulationResult has the solver status, objective_value None and fluxes None, and
    simulate(..., slim=True) returns None. A slim simulation otherwise returns the objective value (a float).
    """

    # methods whose objective value is a minimised quantity (total flux, distance, changed fluxes)
    _MINIMISING_METHODS = (SimulationMethod.pFBA, SimulationMethod.lMOMA, SimulationMethod.MOMA, SimulationMethod.ROOM)

    @classmethod
    def _optimised_sense(cls, method, maximize):
        """The sense of the objective a method optimises: pFBA, (l)MOMA and ROOM always minimise."""
        return False if method in cls._MINIMISING_METHODS else maximize

    def __repr__(self):
        """Rich representation showing simulator state."""
        lines = []
        lines.append("=" * 60)

        # Simulator type
        sim_type = self.__class__.__name__
        lines.append(f"Simulator: {sim_type}")
        lines.append("=" * 60)

        # Model info
        try:
            model_id = self.model.id if hasattr(self.model, "id") else str(self.model)
            lines.append(f"{'Model:':<20} {model_id}")
        except:
            lines.append(f"{'Model:':<20} <not available>")

        # Model statistics
        try:
            lines.append(f"{'Reactions:':<20} {len(self.reactions)}")
        except:
            pass

        try:
            lines.append(f"{'Metabolites:':<20} {len(self.metabolites)}")
        except:
            pass

        try:
            lines.append(f"{'Genes:':<20} {len(self.genes)}")
        except:
            pass

        # Objective
        try:
            obj = self.get_objective()
            if obj:
                if isinstance(obj, dict):
                    obj_ids = list(obj.keys())[:3]  # Show first 3
                    obj_str = ", ".join(str(o) for o in obj_ids)
                    if len(obj) > 3:
                        obj_str += f", ... ({len(obj)} total)"
                    lines.append(f"{'Objective:':<20} {obj_str}")
                else:
                    lines.append(f"{'Objective:':<20} {obj}")
        except:
            pass

        # Environmental conditions
        try:
            env = self.environmental_conditions
            if env and len(env) > 0:
                lines.append(f"{'Medium conditions:':<20} {len(env)} constraints")
        except:
            pass

        # Status
        lines.append(f"{'Status:':<20} Ready")

        lines.append("=" * 60)
        return "\n".join(lines)

    def _repr_html_(self):
        """Pandas-like HTML representation for Jupyter notebooks."""
        from mewpy.util.html_repr import render_html_table

        rows = []

        # Simulator type
        sim_type = self.__class__.__name__

        # Model info
        try:
            model_id = self.model.id if hasattr(self.model, "id") else str(self.model)
            rows.append(("Model", model_id))
        except:
            rows.append(("Model", "<not available>"))

        # Model statistics
        try:
            rows.append(("Reactions", str(len(self.reactions))))
        except:
            pass

        try:
            rows.append(("Metabolites", str(len(self.metabolites))))
        except:
            pass

        try:
            rows.append(("Genes", str(len(self.genes))))
        except:
            pass

        # Objective
        try:
            obj = self.get_objective()
            if obj:
                if isinstance(obj, dict):
                    obj_ids = list(obj.keys())[:3]
                    obj_str = ", ".join(str(o) for o in obj_ids)
                    if len(obj) > 3:
                        obj_str += f", ... ({len(obj)} total)"
                    rows.append(("Objective", obj_str))
                else:
                    rows.append(("Objective", str(obj)))
        except:
            pass

        # Environmental conditions
        try:
            env = self.environmental_conditions
            if env and len(env) > 0:
                rows.append(("Medium conditions", f"{len(env)} constraints"))
        except:
            pass

        # Status
        rows.append(("Status", "Ready"))

        return render_html_table(f"Simulator: {sim_type}", rows)

    def simulate(self, *args, **kwargs):
        """Abstract method to run a phenotype simulation.

        :returns: A SimulationResult.

        """
        method = kwargs.get("method", None)
        if method and callable(method):
            new_kwargs = {k: v for k, v in kwargs.items() if k in ["method"]}
            return self._simulate_callable(method=method, *args, **new_kwargs)
        else:
            return self._simulate(self, *args, **kwargs)

    def _simulate_callable(self, method, *args, **kwargs):
        if not callable(method):
            raise ValueError(f"The method {method} is not callable.")
        simulation_result = method(self, *args, **kwargs)
        return simulation_result

    def _simulate(self, *args, **kwargs):
        """Abstract method to run a phenotype simulation.

        :returns: A SimulationResult.

        """
        raise NotImplementedError

    @abstractmethod
    def FVA(self, reactions=None, obj_frac=0, constraints=None, loopless=False, internal=None, solver=None):
        """Abstract method to run Flux Variability Analysis (FVA).

        :returns: A dictionary of flux range values.

        """
        raise NotImplementedError

    def simulate_mp(
        self,
        objective=None,
        method=SimulationMethod.FBA,
        maximize=True,
        constraints_list=None,
        reference=None,
        solver=None,
        jobs=None,
        desc="Parallel Simulation",
        **kwargs,
    ):
        """Parallel phenotype simulations.

        :param (dict) objective: The simulations objective. If none, the model objective is considered.
        :param (SimulationMethod) method: The SimulationMethod (FBA, pFBA, lMOMA, etc ...)
        :param (boolean) maximize: The optimization direction
        :param (dict) constraints: A dictionary of constraints to be applied to the model.
        :param (dict) reference: A dictionary of reaction flux values.
        :param (Solver) solver: An instance of the solver.
        :param (int) jobs: The number of parallel jobs.
        :param (str) desc: Description to present in tdqm.
        """
        constraints_list = [None] if not constraints_list else constraints_list
        jobs = jobs if jobs else cpu_count()
        logger.debug(f"Using {jobs} jobs")
        from ..util.utilities import tqdm_joblib

        with tqdm_joblib(tqdm(desc=desc, total=len(constraints_list))):
            res = Parallel(n_jobs=jobs)(
                delayed(simulate)(
                    self.model,
                    self.environmental_conditions,
                    objective,
                    method,
                    maximize,
                    constraints,
                    reference,
                    solver,
                    **kwargs,
                )
                for constraints in constraints_list
            )
        return res

    @abstractmethod
    def get_reaction_bounds(self, r_id):
        raise NotImplementedError

    def set_environmental_conditions(self, medium):
        medium = self.resolve_keys(medium)
        for k, v in medium.items():
            if isinstance(v, tuple):
                lb, ub = v
            elif isinstance(v, (float, int)):
                lb, ub = v, v
            else:
                raise ValueError(f"{v} is an inappropiate bound.")
            self.set_reaction_bounds(k, lb, ub, False)
        self._environmental_conditions = medium

    @abstractmethod
    def metabolite_reaction_lookup(self, force_recalculate=False):
        raise NotImplementedError

    def find(self, pattern=None, sort=False, find_in="r"):
        """A user friendly method to find metabolites, reactions or genes in the model.

        :param pattern: The pattern which can be a regular expression, defaults to None in which case all entries are listed.
        :type pattern: str, optional
        :param sort: if the search results should be sorted, defaults to False
        :type sort: bool, optional
        :param find_in: The search set: 'r' for reactions, 'm' for metabolites, 'g' for genes, defaults to 'r'
        :type find_in: str, optional
        :return: the search results
        :rtype: pandas dataframe
        """
        if find_in == "m":
            values = self.metabolites
        elif find_in == "g":
            values = self.genes
        else:
            values = self.reactions
        if pattern:
            import re

            if isinstance(pattern, list):
                patt = "|".join(pattern)
                re_expr = re.compile(patt)
            else:
                re_expr = re.compile(pattern)
            values = [x for x in values if re_expr.search(x) is not None]
        if sort:
            values.sort()

        import pandas as pd

        if find_in == "m":
            data = [self.get_metabolite(x) for x in values]
        elif find_in == "g":
            data = [self.get_gene(x) for x in values]
        else:
            data = [self.get_reaction(x) for x in values]

        if data:
            df = pd.DataFrame(data)
            df = df.set_index(df.columns[0])
        else:
            df = pd.DataFrame()
        return df

    def find_genes(self, pattern=None, sort=False):
        return self.find(pattern=pattern, sort=sort, find_in="g")

    def find_metabolites(self, pattern=None, sort=False):
        return self.find(pattern=pattern, sort=sort, find_in="m")

    def find_metabolite_by_formula(self, formula: Union[str, List[str]]):
        from mewpy.util.utilities import elements

        formulas = [formula] if isinstance(formula, str) else formula
        elems = [elements(f) for f in formulas]
        mets = []
        for met in self.metabolites:
            f = self.get_metabolite(met).formula
            if elements(f) in elems:
                mets.append(met)
        if mets:
            return self.find_metabolites(mets)

    def metabolite_by_formula(self, formula: str, compartment="c"):
        from mewpy.util.utilities import elements

        elem = elements(formula)
        for met in self.metabolites:
            m = self.get_metabolite(met)
            f = m.formula
            c = m.compartment
            if elements(f) == elem and c == compartment:
                return met
        return None

    def find_reactions(self, pattern=None, sort=False):
        return self.find(pattern=pattern, sort=sort, find_in="r")

    def _knockout_growth(self, constraints):
        """Objective value under the constraints, or None when the solver fails (the candidate is non-viable).

        Solver errors (e.g. "SCIP: error in LP solver!") are logged once per essentiality scan; errors that
        point to a bug (KeyError, ValueError, TypeError, AttributeError) are raised.
        """
        try:
            return self.simulate(constraints=constraints, slim=True)
        except (KeyError, ValueError, TypeError, AttributeError):
            raise
        except Exception as e:
            if not getattr(self, "_essentiality_warned", False):
                logger.warning("A simulation failed during an essentiality check; treated as non-viable: %s", e)
                self._essentiality_warned = True
            return None

    def is_essential_reaction(self, rxn, min_growth=0.01):
        """Whether knocking out a reaction leaves a growth (objective value) below min_growth (absolute)."""
        res = self._knockout_growth({self.resolve(rxn, "reaction"): 0})
        return res is None or math.isnan(res) or res < min_growth

    def essential_reactions(self, min_growth=0.01):
        """Essential reactions are those whose knockout leaves a growth (objective value) below min_growth.

        :param float min_growth: Minimal growth, an absolute objective value (not a fraction of the wild type).
            Default 0.01.
        :returns: A list of essential reactions.
        """
        essential = getattr(self, "_essential_reactions", None)
        if essential is not None:
            return essential
        self._essentiality_warned = False
        essential = []
        for rxn in tqdm(self.reactions):
            if self.is_essential_reaction(rxn, min_growth=min_growth):
                essential.append(rxn)
        self._essential_reactions = essential
        return self._essential_reactions

    def evaluate_gprs(self, active_genes):
        """Returns the list of active reactions for a given list of active genes.

        :param list active_genes: List of genes identifiers.
        :returns: A list of active reaction identifiers.

        """
        active_genes = set(self.resolve_many(active_genes, "gene"))
        active_reactions = []
        reactions = self.reactions
        for r_id in reactions:
            gpr = self.get_gpr(r_id)
            if gpr:
                if evaluate_expression_tree(str(gpr), active_genes):
                    active_reactions.append(r_id)
            else:
                active_reactions.append(r_id)
        return active_reactions

    def _compartment_of(self, m_id):
        return self.get_metabolite(m_id)["compartment"]

    def _formula_of(self, m_id):
        return self.get_metabolite(m_id).get("formula") or ""

    def _is_proton(self, m_id):
        """A metabolite with formula H; without a formula, one whose id stem is 'h' (h_c, M_h_e)."""
        formula = self._formula_of(m_id).strip()
        if formula:
            return formula == "H"
        return self.to_canonical(m_id, "metabolite").rsplit("_", 1)[0] == "h"

    def get_transport_reactions(self):
        """Transport reactions: metabolites other than protons in two or more compartments, exchange (boundary)
        reactions excluded. Proton-pumping enzymes (e.g. ATP synthase, NADH dehydrogenase) are not transport.

        Defined here, once, so that every backend agrees.
        :returns: A list of reaction identifiers.
        """
        exchanges = set(self.get_exchange_reactions())
        compartment = {}
        transport = []
        for r_id in self.reactions:
            if r_id in exchanges:
                continue
            compartments = set()
            for m_id in self.get_reaction_metabolites(r_id):
                if m_id not in compartment:
                    compartment[m_id] = None if self._is_proton(m_id) else self._compartment_of(m_id)
                if compartment[m_id] is not None:
                    compartments.add(compartment[m_id])
            if len(compartments) > 1:
                transport.append(r_id)
        return transport

    def get_transport_genes(self):
        """Transport genes: genes all of whose reactions are transport reactions (see get_transport_reactions).

        :returns: A list of gene identifiers.
        """
        transport = set(self.get_transport_reactions())
        return [g for g, rxns in self.get_gene_reactions().items() if rxns and set(rxns) <= transport]

    def is_essential_gene(self, gene, min_growth=0.01):
        """Whether deleting a gene leaves a growth (objective value) below min_growth (absolute)."""
        gene = self.resolve(gene, "gene")
        gr = self.get_gene_reactions()
        active_genes = set(self.genes) - {gene}
        rxns = gr.get(gene, [])
        if not rxns:
            return False
        inactive_reactions = []
        for r_id in rxns:
            gpr = self.get_gpr(r_id)
            if gpr and not evaluate_expression_tree(str(gpr), active_genes):
                inactive_reactions.append(r_id)
        constraints = {rxn: 0 for rxn in inactive_reactions}
        res = self._knockout_growth(constraints)
        return res is None or math.isnan(res) or res < min_growth

    def essential_genes(self, min_growth=0.01):
        """Essential genes are those whose deletion leaves a growth (objective value) below min_growth.

        :param float min_growth: Minimal growth, an absolute objective value (not a fraction of the wild type).
            Default 0.01.
        :returns: A list of essential genes.

        """
        essential = getattr(self, "_essential_genes", None)
        if essential is not None:
            return essential
        self._essentiality_warned = False
        essential = []
        for g in tqdm(self.genes):
            if self.is_essential_gene(g, min_growth=min_growth):
                essential.append(g)
        self._essential_genes = essential
        return self._essential_genes

    @property
    def reference(self):
        """The reference wild type reaction flux values.

        :returns: A dictionary of wild type reaction flux values.

        """
        ref = getattr(self, "_reference", None)
        if ref is not None:
            return ref
        self._reference = self.simulate(method="pFBA").fluxes
        return self._reference

    def create_empty_model(self, model_id: str):
        raise NotImplementedError("Subclasses must implement create_empty_model()")

    def get_external_metabolites(self):
        external = []
        ext_com = [c_id for c_id in self.compartments if self.get_compartment(c_id).external]
        for m_id in self.metabolites:
            if self.get_metabolite(m_id).compartment in ext_com:
                external.append(m_id)
        return external

    def blocked_reactions(self, constraints=None, reactions=None, abstol=1e-9):
        """Find all blocked reactions in a model

        :param (dict) constraints: additional constraints (optional)
        :param (list) reactions: List of reactions which will be tested (default: None, test all reactions)
        :param (float) abstol: absolute tolerance (default: 1e-9)

        :returns: A list blocked reactions
        """

        variability = self.FVA(obj_frac=0, reactions=reactions, constraints=constraints)

        return [r_id for r_id, (lb, ub) in variability.items() if (abs(lb) + abs(ub)) < abstol]

    def get_metabolite_reactions(self, m_id: str) -> List[str]:
        """Returns the list or reactions that produce or consume a metabolite.

        :param m_id: the metabolite identifier
        :return: A list of reaction identifiers
        """
        m_r = self.metabolite_reaction_lookup()
        return list(m_r[self.resolve(m_id, "metabolite")].keys())

    def get_metabolite_producers(self, m_id: str) -> List[str]:
        """Returns the list or reactions that produce a metabolite.

        :param m_id: the metabolite identifier
        :return: A list of reaction identifiers
        """
        m_r = self.metabolite_reaction_lookup()
        return [k for k, v in m_r[self.resolve(m_id, "metabolite")].items() if v > 0]

    def get_metabolite_consumers(self, m_id):
        """Returns the list or reactions that consume a metabolite.

        :param m_id: the metabolite identifier
        :return: A list of reaction identifiers
        """
        m_r = self.metabolite_reaction_lookup()
        return [k for k, v in m_r[self.resolve(m_id, "metabolite")].items() if v < 0]

    def copy(self):
        """Retuns a copy of the Simulator instance."""
        return deepcopy(self)

    def __getstate__(self):
        # A cached LP solver (reset_solver=False) wraps a native object (e.g. SCIP) that cannot be
        # copied or pickled; the copy builds its own solver on the next simulation.
        state = self.__dict__.copy()
        if state.get("solver") is not None:
            state["solver"] = None
        return state


class SimulationResult(object):
    """Class that represents simulation results and performs operations over them."""

    def __init__(
        self,
        model,
        objective_value,
        fluxes=None,
        status=None,
        envcond=None,
        model_constraints=None,
        simul_constraints=None,
        maximize=True,
        method=None,
        shadow_prices=None,
    ):
        """
        :param model: A model instance.
        :param objective_value: The phenotype simulation objective value.
        :param dict fluxes: A dictionary of reaction fluxes values.
        :param status: The LP status.
        :param dict envcond: The environmental conditions of the phenotype simulation.
        :param dict model_constraints: Possible persistent additional constraints.
        :param dict simul_constraints: The simulation constraints.
        :param boolean maximize: Optimization direction.
        :param SimulationMethod method: The phenotypic methos
        :param dict shadow_prices: shadow prices
        """
        self.model = model
        self.objective_value = objective_value
        self.status = status
        self.fluxes = fluxes
        self.maximize = maximize
        # Environmental conditions
        self.envcond = envcond if envcond else OrderedDict()
        # Persistent additional constraints not included in the environmental conditions
        self.model_constraints = model_constraints if model_constraints else OrderedDict()
        # Constraints specific to the simulation
        self.simulation_constraints = simul_constraints if simul_constraints else OrderedDict()
        self.method = method
        self.shadow_prices = shadow_prices

    def get_constraints(self):
        """
        :returns: All constraints applied during the simulation both persistent and simulation specific.
        """
        constraints = OrderedDict()
        constraints.update(self.envcond)
        constraints.update(self.model_constraints)
        constraints.update(self.simulation_constraints)
        return constraints

    def __repr__(self):
        """Rich representation showing simulation result details."""
        lines = []
        lines.append("=" * 60)
        lines.append("Simulation Result")
        lines.append("=" * 60)

        # Status
        try:
            if self.status:
                lines.append(f"{'Status:':<20} {self.status}")
        except:
            pass

        # Objective value with direction
        try:
            if self.objective_value is not None:
                direction = "maximize" if self.maximize else "minimize"
                label = f"Objective ({direction}):"
                lines.append(f"{label:<20} {self.objective_value:.6g}")
        except:
            pass

        # Method
        try:
            if self.method:
                if callable(self.method):
                    method_name = getattr(self.method, "__name__", repr(self.method))
                else:
                    method_name = str(self.method)
                lines.append(f"{'Method:':<20} {method_name}")
        except:
            pass

        # Model info
        try:
            if self.model:
                if hasattr(self.model, "id"):
                    model_id = self.model.id
                else:
                    model_id = str(self.model)[:30]
                lines.append(f"{'Model:':<20} {model_id}")
        except:
            pass

        # Fluxes summary
        try:
            if self.fluxes:
                total_fluxes = len(self.fluxes)
                non_zero_fluxes = sum(1 for v in self.fluxes.values() if abs(v) > 1e-9)
                lines.append(f"{'Fluxes:':<20} {non_zero_fluxes} non-zero / {total_fluxes} total")
        except:
            pass

        # Constraints summary
        try:
            all_constraints = self.get_constraints()
            if all_constraints:
                constraint_count = len(all_constraints)
                lines.append(f"{'Constraints:':<20} {constraint_count}")

                # Break down by type
                env_count = len(self.envcond) if self.envcond else 0
                model_count = len(self.model_constraints) if self.model_constraints else 0
                simul_count = len(self.simulation_constraints) if self.simulation_constraints else 0

                if env_count > 0:
                    lines.append(f"{'  Environment:':<20} {env_count}")
                if model_count > 0:
                    lines.append(f"{'  Model:':<20} {model_count}")
                if simul_count > 0:
                    lines.append(f"{'  Simulation:':<20} {simul_count}")
        except:
            pass

        # Shadow prices
        try:
            if self.shadow_prices:
                sp_count = len(self.shadow_prices)
                lines.append(f"{'Shadow prices:':<20} {sp_count} available")
        except:
            pass

        lines.append("=" * 60)
        return "\n".join(lines)

    def _repr_html_(self):
        """Pandas-like HTML representation for Jupyter notebooks."""
        from mewpy.util.html_repr import render_html_table

        rows = []

        # Status
        if self.status:
            rows.append(("Status", str(self.status)))

        # Objective value with direction
        if self.objective_value is not None:
            direction = "maximize" if self.maximize else "minimize"
            rows.append((f"Objective ({direction})", f"{self.objective_value:.6g}"))

        # Method
        if self.method:
            if callable(self.method):
                method_name = getattr(self.method, "__name__", repr(self.method))
            else:
                method_name = str(self.method)
            rows.append(("Method", method_name))

        # Model info
        if self.model:
            if hasattr(self.model, "id"):
                model_id = self.model.id
            else:
                model_id = str(self.model)[:30]
            rows.append(("Model", model_id))

        # Fluxes summary
        if self.fluxes:
            total_fluxes = len(self.fluxes)
            non_zero_fluxes = sum(1 for v in self.fluxes.values() if abs(v) > 1e-9)
            rows.append(("Fluxes", f"{non_zero_fluxes} non-zero / {total_fluxes} total"))

        # Constraints summary
        all_constraints = self.get_constraints()
        if all_constraints:
            constraint_count = len(all_constraints)
            rows.append(("Constraints", str(constraint_count)))

            env_count = len(self.envcond) if self.envcond else 0
            model_count = len(self.model_constraints) if self.model_constraints else 0
            simul_count = len(self.simulation_constraints) if self.simulation_constraints else 0

            if env_count > 0:
                rows.append(("  Environment", str(env_count)))
            if model_count > 0:
                rows.append(("  Model", str(model_count)))
            if simul_count > 0:
                rows.append(("  Simulation", str(simul_count)))

        # Shadow prices
        if self.shadow_prices:
            sp_count = len(self.shadow_prices)
            rows.append(("Shadow prices", f"{sp_count} available"))

        return render_html_table("Simulation Result", rows)

    def __str__(self):
        return self.__repr__()

    def find(self, pattern=None, sort=False, shadow_prices=False, show_nulls=False):
        """Returns a dataframe of reactions and their fluxes matching a pattern or a list of patterns.

        :param pattern: a string or a list of strings. defaults to None
        :param sort: If the dataframe is to be sorted by flux rates. Defaults to False
        :type sort: bool, optional
        :return: returns a dataframe.
        :rtype: pandas.DataFrame
        """
        if shadow_prices:
            try:
                values = [(key, value) for key, value in self.shadow_prices.items()]
                columns = ["Metabolite", "Shadow Price"]
            except Exception:
                raise ValueError("No shadow prices")
        else:
            try:
                if show_nulls:
                    values = [(key, value) for key, value in self.fluxes.items()]
                else:
                    values = [(key, value) for key, value in self.fluxes.items() if value != 0.0]
                columns = ["Reaction ID", "Flux rate"]
            except Exception:
                raise ValueError("No fluxes")
        if pattern:
            import re

            if isinstance(pattern, list):
                patt = "|".join(pattern)
                re_expr = re.compile(patt)
            else:
                re_expr = re.compile(pattern)
            values = [x for x in values if re_expr.search(x[0]) is not None]
        if sort:
            values.sort(key=lambda x: x[1])
        import pandas as pd

        df = pd.DataFrame(values, columns=columns)
        df = df.set_index(columns[0])
        return df

    @property
    def dataframe(self):
        return self.find()

    def get_net_conversion(self, biomassId=None):
        """Returns a string representation of the net conversion.

        :params str biosmassId: Biomass identifier (optional)
        :returns: A string representation of the net conversion.

        """

        left = ""
        right = ""
        firstLeft, firstRight = True, True
        from . import get_simulator

        sim = get_simulator(self.model)
        ssFluxes = self.fluxes
        for r_id in sim.reactions:
            fluxValue = ssFluxes[r_id]
            sub = list(sim.get_substrates(r_id).keys())
            prod = list(sim.get_products(r_id).keys())
            # if rId is a drain reaction
            if fluxValue != 0.0 and (len(sub) == 0 or len(prod) == 0):
                m = sub + prod
                if fluxValue < 0:
                    if firstLeft:
                        firstLeft = False
                    else:
                        left = left + " + "
                    left = left + str(-1 * fluxValue)
                    left = left + " " + m[0]
                else:
                    if firstRight:
                        firstRight = False
                    else:
                        right = right + " + "
                    right = right + str(fluxValue)
                    right = right + " " + m[0]

        biomassId = sim.resolve(biomassId, "reaction") if biomassId else biomassId
        if biomassId and biomassId in ssFluxes.keys():
            biomassFlux = ssFluxes[biomassId]
            if biomassFlux > 0:
                if firstRight:
                    firstRight = False
                else:
                    right = right + " + "
                right = right + str(biomassFlux)
                right = right + " " + biomassId

        return left + " --> " + right

    def get_metabolites_turnover(self, pattern=None, format="df"):
        """Calculate metabolite turnovers.

         :param str format: the display format (pandas.DataFrame or dict).
             Default 'df' pandas.DataFrame

        Returns:
             dict or pandas.DataFrame: metabolite turnover rates
        """
        from . import get_simulator

        sim = get_simulator(self.model)

        if not self.fluxes:
            return None

        mets = None
        if pattern:
            import re

            if isinstance(pattern, list):
                patt = "|".join(pattern)
                re_expr = re.compile(patt)
            else:
                re_expr = re.compile(pattern)
            mets = [x for x in sim.metabolites if re_expr.search(x) is not None]

        m_r_table = sim.metabolite_reaction_lookup()

        data = {
            m_id: 0.5 * sum([abs(coeff * self.fluxes[r_id]) for r_id, coeff in neighbours.items()])
            for m_id, neighbours in m_r_table.items()
        }
        if mets is not None:
            data = {k: v for k, v in data.items() if k in mets}

        if format == "df":
            import pandas as pd

            df = pd.DataFrame(list(data.values()), index=list(data.keys()), columns=["Turnover"])
            df.index.name = "Metabolite"
            return df
        return data

    def get_metabolite(self, met_id, format="df"):
        """Displays the consumption/production of a metabolite in the reactions
        it participates.

        :param str met_id: the metabolite identifier.
        :param str format: the display format (pandas.DataFrame or dict).
             Default 'df' pandas.DataFrame

        Returns:
            dict or pandas.DataFrame: metabolite turnover rates
        """
        from . import get_simulator

        sim = get_simulator(self.model)

        if not self.fluxes:
            return None

        m_r = sim.metabolite_reaction_lookup()[sim.resolve(met_id, "metabolite")]
        data = {r_id: coeff * self.fluxes[r_id] for r_id, coeff in m_r.items()}
        if format == "df":
            import pandas as pd

            df = pd.DataFrame(list(data.values()), index=list(data.keys()), columns=["Value"])
            df.index.name = "Reaction"
            return df
        return data

    @classmethod
    def from_linear_solver(cls, solution):
        """Converts a solver Solution object to a SolutionResult object.

        :param solution: solution to be converted
        :type solution: mewpy.solver.solution.Solution
        :raises ValueError: if solution is not an instance of mewpy.solver.solution.Solution
        :return: an instance of SolutionResult
        :rtype: SolutionResult
        """
        from mewpy.solvers import Solution, Status

        smap = {
            Status.OPTIMAL: SStatus.OPTIMAL,
            Status.UNKNOWN: SStatus.UNKNOWN,
            Status.SUBOPTIMAL: SStatus.SUBOPTIMAL,
            Status.UNBOUNDED: SStatus.UNBOUNDED,
            Status.INFEASIBLE: SStatus.INFEASIBLE,
            Status.INF_OR_UNB: SStatus.INF_OR_UNB,
        }
        if not isinstance(solution, Solution):
            raise ValueError("solution should be and instance of mewpy.solvers.solution.Solution")
        return cls(None, solution.fobj, fluxes=solution.values, status=smap[solution.status])


def simulate(
    model,
    envcond=None,
    objective=None,
    method=SimulationMethod.FBA,
    maximize=True,
    constraints=None,
    reference=None,
    solver=None,
    **kwargs,
):
    """Runs an FBA phenotype simulation.

    :param model: cobrapy, reframed, GERM constraint-base model
    :param (dict ) envcond : Environmental conditions, defaults to None
    :param (dict) objective: the FBA objective , defaults to None
    :param method: The FBA method, defaults to SimulationMethod.FBA
    :param maximize: optimization sense, defaults to True
    :param (dict) constraints: additional constraints , defaults to None
    :param (dict) reference: reference fluxes, defaults to None
    :param solver: a solver instance, defaults to None

    :returns: SimultationResult
    """
    from . import get_simulator

    sim = get_simulator(model, envcond=envcond)
    res = sim.simulate(
        objective=objective,
        method=method,
        maximize=maximize,
        constraints=constraints,
        reference=reference,
        solver=solver,
        **kwargs,
    )
    return res
