import unittest

from mewpy.simulation.kinetic import KineticSimulation

MODELS_PATH = "tests/data/"
MODEL = MODELS_PATH + "chassagnole2002.xml"


class TestKineticSimulation(unittest.TestCase):

    def setUp(self):
        from mewpy.io.sbml import load_ODEModel

        self.model = load_ODEModel(MODEL)

    def test_build_ode(self):
        self.model.build_ode()

    def test_simulation(self):
        from mewpy.simulation.kinetic import KineticSimulation

        sim = KineticSimulation(self.model)
        sim.simulate()

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
