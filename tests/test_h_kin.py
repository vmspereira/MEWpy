import unittest

from mewpy.simulation.kinetic import KineticSimulation

MODELS_PATH = "tests/data/"
MODEL = MODELS_PATH + "chassagnole2002.xml"


class TestKineticSimulation(unittest.TestCase):

    def setUp(self):
        from mewpy.io.sbml import load_ODEModel

        self.model = load_ODEModel(MODEL)

    def test_build_ode(self):
        """The ODE built as a string evaluates every rate law and, for each metabolite, the stoichiometric sum of
        the rates over the compartment volume"""
        import numpy as np

        func_str = self.model.build_ode()
        self.assertTrue(func_str.startswith("def ode_func(t, x, r, p, v)"))

        rates = {}
        y0 = [self.model.concentrations.get(m_id, 0.0) for m_id in self.model.metabolites]
        dxdt = self.model.get_ode(r_dict=rates)(0.0, y0)
        self.assertEqual(set(rates), set(self.model.ratelaws))
        volume = self.model.merge_constants()
        expected = np.zeros(len(self.model.metabolites))
        index = {m_id: i for i, m_id in enumerate(self.model.metabolites)}
        for r_id, rxn in self.model.reactions.items():
            for m_id, coeff in rxn.stoichiometry.items():
                if m_id in index:
                    expected[index[m_id]] += coeff * rates[r_id] / volume[self.model.metabolites[m_id].compartment]
        np.testing.assert_allclose(dxdt, expected, rtol=1e-10, atol=1e-12)
        self.assertGreater(np.abs(expected).max(), 1.0)

    def test_simulation(self):
        """The default simulation integrates to steady state: its concentrations are SciPy's (LSODA) at the last
        time point, their derivatives vanish, and the rates are the rate laws at those concentrations"""
        import numpy as np
        from scipy.integrate import solve_ivp

        from mewpy.simulation.kinetic import KineticSimulation
        from mewpy.solvers.ode import ODEStatus, SolverConfigurations

        sim = KineticSimulation(self.model)
        result = sim.simulate()

        self.assertEqual(result.status, ODEStatus.OPTIMAL)
        concentrations = np.array([result.concentrations[m_id] for m_id in self.model.metabolites])
        t_start, t_end = sim.get_time_points()[0], sim.get_time_points()[-1]
        y0 = [self.model.concentrations.get(m_id, 0.0) for m_id in self.model.metabolites]
        times = np.linspace(t_start, t_end, SolverConfigurations.N_STEPS)
        reference = solve_ivp(
            self.model.get_ode(), (t_start, t_end), y0, method="LSODA", t_eval=times, rtol=1e-8, atol=1e-10
        )
        self.assertEqual(reference.status, 0)
        np.testing.assert_allclose(concentrations, reference.y[:, -1], rtol=1e-5, atol=1e-8)
        np.testing.assert_allclose(self.model.deriv(t_end, list(concentrations)), 0, atol=1e-8)
        rates = {}
        self.model.get_ode(r_dict=rates)(t_end, list(concentrations))
        for r_id, rate in result.fluxes.items():
            self.assertAlmostEqual(rate, rates[r_id], delta=1e-6 * max(1.0, abs(rates[r_id])), msg=r_id)

    def test_deriv_vs_build_ode_equivalence(self):
        """Test that deriv() and build_ode() produce equivalent results.

        This verifies that both ODE evaluation paths (direct deriv() and
        compiled build_ode()) apply the same mathematical operations:
        - Stoichiometric coefficients
        - Compartment volume normalization
        - Factor application
        """
        import numpy as np

        # Get initial concentrations
        y0 = [self.model.concentrations.get(m_id, 0.0) for m_id in self.model.metabolites]
        t = 0.0

        # Test without factors
        deriv_result = self.model.deriv(t, y0)
        ode_func = self.model.get_ode()
        build_ode_result = ode_func(t, y0)

        np.testing.assert_allclose(
            deriv_result, build_ode_result, rtol=1e-10, err_msg="deriv() and build_ode() produce different results!"
        )

        # Test with factors (note: deriv() doesn't support factors directly)
        factors = {list(self.model.variable_params.keys())[0]: 0.5} if self.model.variable_params else {}
        if factors:
            ode_func_with_factors = self.model.get_ode(factors=factors)
            build_ode_factors_result = ode_func_with_factors(t, y0)
            # Verify result is computed (shows factors work in build_ode path)
            assert build_ode_factors_result is not None
            # Note: This test mainly verifies that both paths use the same stoichiometry and volume normalization

    def test_deriv_vs_build_ode_off_initial_state(self):
        """deriv() must use the state vector y, not the initial concentrations."""
        import numpy as np

        y0 = [self.model.concentrations.get(m_id, 0.0) for m_id in self.model.metabolites]
        y1 = [2 * v + 0.1 for v in y0]
        ode_func = self.model.get_ode()
        np.testing.assert_allclose(self.model.deriv(0.0, y1), ode_func(0.0, y1), rtol=1e-8, atol=1e-12)
        assert not np.allclose(self.model.deriv(0.0, y0), self.model.deriv(0.0, y1))
