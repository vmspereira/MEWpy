"""
Comprehensive validation tests for RFBA and SRFBA implementations.

Tests verify that:
1. RFBA correctly applies regulatory constraints
2. SRFBA integrates Boolean logic into MILP
3. Results are consistent with expected behavior
4. Regulatory network properly constrains metabolic fluxes

Oracles: COBRApy FBA on the same SBML file, with the reactions whose GPR (COBRApy) is false for the
decoded gene state switched off; the documented RFBA growth on E. coli core (0.8514, CHANGELOG); and the
textbook anaerobic growth when ArcA or Fnr are forced active.
"""

import sys
from functools import lru_cache
from pathlib import Path

import pytest

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

MODEL_PATH = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core.xml"
REG_PATH = Path(__file__).parent.parent / "examples" / "models" / "germ" / "e_coli_core_trn.csv"
# e_coli_core_trn.csv: header row, columns bNum, Gene (alias), Rule
TRN_CSV = dict(regulatory_format="csv", sep=",", id_col=0, rule_col=2, aliases_cols=[1], header=0)
# RFBA growth on E. coli core with the default regulatory state (CHANGELOG, "as in the documentation")
RFBA_GROWTH = 0.8514
# Orth et al. (2010): anaerobic growth of e_coli_core on glucose
ANAEROBIC_GROWTH = 0.2117


def _integrated_model():
    from mewpy.germ.models import RegulatoryExtension

    return RegulatoryExtension.from_sbml(str(MODEL_PATH), str(REG_PATH), flavor="reframed", **TRN_CSV)


@lru_cache(maxsize=None)
def _cobra_model():
    import cobra

    return cobra.io.read_sbml_model(str(MODEL_PATH))


def _cobra_growth(genes_off=(), anaerobic=False):
    """COBRApy FBA with the reactions whose GPR is false for the switched-off genes closed"""
    model = _cobra_model()
    knockouts = {g[2:] if g.startswith("G_") else g for g in genes_off}
    with model:
        for rxn in model.reactions:
            if knockouts & {g.id for g in rxn.genes} and not rxn.gpr.eval(knockouts=knockouts):
                rxn.bounds = (0, 0)
        if anaerobic:
            model.reactions.EX_o2_e.lower_bound = 0
        return model.slim_optimize(error_value=0.0)


def _rfba_genes_off(rfba, model, initial_state=None):
    """Metabolic genes off in RFBA: regulators that are targets follow their rules, then the targets are
    decoded from that regulatory state; explicit states win"""
    initial_state = initial_state or {}
    state = rfba.initial_state(initial_state)
    regulators = set(model.regulators)
    regulatory = {k: v for k, v in rfba.decode_regulatory_state(state).items() if k in regulators}
    metabolic = {**rfba.decode_metabolic_state({**state, **regulatory, **initial_state}), **initial_state}
    return {g for g, v in metabolic.items() if g in model.genes and v == 0}


class TestRFBAValidation:
    """Validation tests for RFBA implementation."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        return _integrated_model()

    def test_rfba_basic_functionality(self, integrated_model):
        """RFBA with the default regulatory state: the documented 0.8514, FBA without the reactions it closes"""
        from mewpy.germ.analysis import RFBA
        from mewpy.solvers.solution import Status

        # Create RFBA instance
        rfba = RFBA(integrated_model)

        # Test that it builds without errors
        rfba.build()

        # Optimize with the default initial state
        solution = rfba.optimize()

        assert solution.status == Status.OPTIMAL
        assert solution.objective_value == pytest.approx(RFBA_GROWTH, abs=1e-4)
        genes_off = _rfba_genes_off(rfba, integrated_model)
        assert genes_off
        assert solution.objective_value == pytest.approx(_cobra_growth(genes_off), abs=1e-6)

    def test_rfba_with_inactive_regulators(self, integrated_model):
        """RFBA with a regulator inactive is FBA without the reactions the decoded gene state closes"""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()

        regulators = list(integrated_model.regulators.keys())
        initial_state = {regulators[0]: 0.0}

        solution = rfba.optimize(initial_state=initial_state)

        expected = _cobra_growth(_rfba_genes_off(rfba, integrated_model, initial_state))
        assert solution.objective_value == pytest.approx(expected, abs=1e-6)

    def test_rfba_dynamic_mode(self, integrated_model):
        """Dynamic RFBA iterates until the regulatory state no longer changes: the last two steps agree and
        start from the steady-state RFBA solution"""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()

        solution = rfba.optimize(dynamic=True)

        growth = [sol.objective_value for sol in solution.solutions.values()]
        assert 2 <= len(growth) < 100
        assert growth[0] == pytest.approx(rfba.optimize().objective_value, abs=1e-9)
        assert growth[-1] == pytest.approx(growth[-2], abs=1e-9)
        assert growth[-1] == pytest.approx(RFBA_GROWTH, abs=1e-4)

    def test_rfba_regulatory_constraints_applied(self, integrated_model):
        """With every regulator inactive, RFBA is FBA without the reactions the decoded gene state closes"""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()

        regulators = list(integrated_model.regulators.keys())
        all_inactive_state = {reg: 0.0 for reg in regulators}

        solution = rfba.optimize(initial_state=all_inactive_state)

        expected = _cobra_growth(_rfba_genes_off(rfba, integrated_model, all_inactive_state))
        assert solution.objective_value == pytest.approx(expected, abs=1e-6)
        assert solution.objective_value <= _cobra_growth() + 1e-6

    def test_rfba_gene_knockout_is_lethal(self, integrated_model):
        """An explicit gene state in initial_state is applied: G_b1779 is the only gene of the essential GAPD."""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()
        assert rfba.optimize().objective_value > 0.1
        knockout = rfba.optimize(initial_state={"G_b1779": 0.0})
        assert (knockout.objective_value or 0.0) == pytest.approx(0.0, abs=1e-6)

    def test_rfba_decode_methods(self, integrated_model):
        """The decode methods: a state for every regulator, a coefficient for every target, and closed bounds
        for exactly the reactions whose GPR (COBRApy) is false for the decoded gene state"""
        from mewpy.germ.analysis import RFBA

        rfba = RFBA(integrated_model).build()

        initial_state = rfba.initial_state()
        assert set(initial_state) == set(integrated_model.regulators)

        regulatory_state = rfba.decode_regulatory_state(initial_state)
        assert set(regulatory_state) == set(integrated_model.targets)
        assert set(regulatory_state.values()) <= {0.0, 1.0}

        metabolic_state = rfba.decode_metabolic_state(initial_state)
        assert metabolic_state == regulatory_state

        constraints = rfba.decode_constraints(metabolic_state)
        off = {g[2:] for g, v in metabolic_state.items() if g in integrated_model.genes and v == 0}
        expected = {
            "R_" + r.id
            for r in _cobra_model().reactions
            if off & {g.id for g in r.genes} and not r.gpr.eval(knockouts=off)
        }
        assert set(constraints) == expected
        assert set(constraints.values()) <= {(0.0, 0.0)}


class TestSRFBAValidation:
    """Validation tests for SRFBA implementation."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        return _integrated_model()

    def test_srfba_basic_functionality(self, integrated_model):
        """SRFBA: one boolean per reaction with a GPR; with free regulators the optimum is FBA's"""
        from mewpy.germ.analysis import SRFBA
        from mewpy.solvers.solution import Status

        # Create SRFBA instance
        srfba = SRFBA(integrated_model)

        # Test that it builds without errors
        srfba.build()

        gpr_reactions = {"R_" + r.id for r in _cobra_model().reactions if r.gpr.body is not None}
        assert srfba._boolean_variables == {"bool_" + r: r for r in gpr_reactions}

        # Optimize
        solution = srfba.optimize()

        assert solution.status == Status.OPTIMAL
        assert solution.objective_value == pytest.approx(_cobra_growth(), abs=1e-6)

    def test_srfba_builds_gpr_constraints(self, integrated_model):
        """The solver holds the mass balances plus, for every reaction with a GPR, a row v <= y * ub (ub != 0) and
        a row v >= y * lb (lb != 0) linking the flux to the reaction boolean y"""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model).build()

        constraints = set(srfba.solver.list_constraints())
        model = _cobra_model()
        assert {"M_" + m.id for m in model.metabolites} <= constraints
        rows = 0
        for rxn in model.reactions:
            if rxn.gpr.body is None:
                continue
            expected = {f"gpr_upper_R_{rxn.id}"} if rxn.upper_bound != 0 else set()
            expected |= {f"gpr_lower_R_{rxn.id}"} if rxn.lower_bound != 0 else set()
            assert expected and expected <= constraints, rxn.id
            rows += len(expected)
        assert len(constraints) > len(model.metabolites) + rows

    def test_srfba_builds_regulatory_constraints(self, integrated_model):
        """Every gene, regulator and reaction boolean of the integrated model is a solver variable"""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model).build()

        variables = set(srfba.solver.list_variables())
        assert len(list(integrated_model.yield_interactions())) == 159
        assert set(integrated_model.genes) <= variables
        assert set(integrated_model.regulators) <= variables
        assert set(srfba._boolean_variables) <= variables

    def test_srfba_with_initial_state(self, integrated_model):
        """Fixing regulators restricts SRFBA; forcing ArcA (b4401) active represses aerobic respiration, so the
        optimum is the anaerobic one (COBRApy, textbook 0.2117)"""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model).build()
        free = srfba.optimize().objective_value

        regulators = list(integrated_model.regulators.keys())
        initial_state = {regulators[0]: 0.0}

        solution = srfba.optimize(initial_state=initial_state)
        assert solution.objective_value <= free + 1e-6

        arca = srfba.optimize(initial_state={"b4401": 1.0})
        assert arca.objective_value == pytest.approx(_cobra_growth(anaerobic=True), abs=1e-6)
        assert arca.objective_value == pytest.approx(ANAEROBIC_GROWTH, abs=5e-5)

    def test_srfba_essential_gene_knockout_is_lethal(self, integrated_model):
        """The SRFBA MILP reaches the solver: wild type equals FBA, the GAPD gene knockout is lethal."""
        from mewpy.germ.analysis import SRFBA

        fba = integrated_model.simulator.simulate().objective_value
        wild_type = SRFBA(integrated_model).build().optimize()
        assert wild_type.objective_value == pytest.approx(fba, abs=1e-6)
        knockout = SRFBA(integrated_model).build().optimize(initial_state={"G_b1779": 0.0})
        assert (knockout.objective_value or 0.0) == pytest.approx(0.0, abs=1e-6)

    def test_srfba_integer_variables(self, integrated_model):
        """The SRFBA solver is a MILP: the reaction booleans are binary variables"""
        from mewpy.germ.analysis import SRFBA

        srfba = SRFBA(integrated_model).build()

        problem = srfba.solver.problem
        binary = {v.name for v in problem.getVars() if v.vtype() == "BINARY"}
        assert set(srfba._boolean_variables) <= binary
        assert len(srfba._boolean_variables) == 69


class TestRFBAvsSRFBA:
    """Compare RFBA and SRFBA results."""

    @pytest.fixture
    def integrated_model(self):
        """Create an integrated model for testing."""
        return _integrated_model()

    def test_compare_basic_results(self, integrated_model):
        """SRFBA chooses the regulatory state, RFBA fixes it: RFBA (0.8514) <= SRFBA (FBA)"""
        from mewpy.germ.analysis import RFBA, SRFBA
        from mewpy.solvers.solution import Status

        # Run RFBA
        rfba = RFBA(integrated_model)
        rfba_solution = rfba.optimize()

        # Run SRFBA
        srfba = SRFBA(integrated_model)
        srfba_solution = srfba.optimize()

        assert rfba_solution.status == Status.OPTIMAL
        assert srfba_solution.status == Status.OPTIMAL
        assert rfba_solution.objective_value == pytest.approx(RFBA_GROWTH, abs=1e-4)
        assert srfba_solution.objective_value == pytest.approx(_cobra_growth(), abs=1e-6)
        assert rfba_solution.objective_value < srfba_solution.objective_value - 1e-3

    def test_compare_with_metabolic_only(self, integrated_model):
        """Regulation only restricts: RFBA and SRFBA never exceed FBA"""
        from mewpy.germ.analysis import RFBA, SRFBA

        # Get FBA result (no regulatory constraints)
        fba_obj = integrated_model.simulator.simulate().objective_value

        rfba_obj = RFBA(integrated_model).optimize().objective_value
        srfba_obj = SRFBA(integrated_model).optimize().objective_value

        assert fba_obj == pytest.approx(_cobra_growth(), abs=1e-6)
        assert rfba_obj <= fba_obj + 1e-6
        assert srfba_obj <= fba_obj + 1e-6


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
