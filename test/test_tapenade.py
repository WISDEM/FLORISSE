import unittest
from openmdao.api import pyOptSparseDriver, Problem, Group, IndepVarComp
from florisse.floris import Floris
import numpy as np


class TotalDerivTestsFlorisAEP(unittest.TestCase):

    def setUp(self):
        print 'Test FLORIS TAPENADE derivatives'
        nTurbines = 3
        self.rtol = 1E-6
        self.atol = 1E-6

        np.random.seed(seed=10)

        turbineX = np.random.rand(nTurbines)*3000.
        turbineY = np.random.rand(nTurbines)*50.
        hubHeight = np.random.rand(nTurbines)*100.+50.

        # initialize input variable arrays
        rotorDiameter = np.ones(nTurbines)*np.random.random()*150.
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        Ct = np.ones(nTurbines)*np.random.random()
        yaw = np.random.rand(nTurbines)*60. - 30.

        # Define flow properties
        nDirections = 1
        wind_speed = float(np.random.rand(1)*10.)       # m/s
        air_density = 1.1716    # kg/m^3
        windDirections = np.random.rand(nDirections)*360.0
        windFrequencies = np.random.rand(nDirections)

        shearExp = 0.15
        z0 = 0.
        zref = 50.

        prob = Problem()
        root = prob.root = Group()

        root.add('turbineXw', IndepVarComp('turbineXw', turbineX), promotes=['*'])
        root.add('turbineYw', IndepVarComp('turbineYw', turbineY), promotes=['*'])
        root.add('yaw0', IndepVarComp('yaw0', yaw), promotes=['*'])
        root.add('hubHeight', IndepVarComp('hubHeight', hubHeight), promotes=['*'])
        root.add('rotorDiameter', IndepVarComp('rotorDiameter', rotorDiameter), promotes=['*'])


        root.add('floris', Floris(nTurbines, differentiable=True, use_rotor_components=False, nSamples=0,
                 verbose=False),promotes=['*'])

        # set up optimizer
        prob.driver = pyOptSparseDriver()
        prob.driver.options['optimizer'] = 'SNOPT'
        prob.driver.add_objective('wtVelocity0', scaler=1.0)

        # select design variables
        prob.driver.add_desvar('turbineXw', scaler=1.0)
        prob.driver.add_desvar('turbineYw', scaler=1.0)
        prob.driver.add_desvar('hubHeight', scaler=1.0)
        prob.driver.add_desvar('yaw0', scaler=1.0)
        prob.driver.add_desvar('rotorDiameter', scaler=1.0)

        prob.root.ln_solver.options['single_voi_relevance_reduction'] = True

        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
        prob['turbineXw'] = turbineX
        prob['turbineYw'] = turbineY
        prob['hubHeight'] = hubHeight
        prob['yaw0'] = yaw
        prob['rotorDiameter'] = rotorDiameter
        prob['axialInduction'] = axialInduction
        prob['Ct'] = Ct
        prob['wind_speed'] = wind_speed
        prob['axialInduction'] = axialInduction
        prob['floris_params:shearExp'] = shearExp
        prob['floris_params:z_ref'] = zref
        prob['floris_params:z0'] = z0

        # run problem
        prob.run_once()

        print prob['wtVelocity0']

        # pass results to self for use with unit test
        self.J = prob.check_total_derivatives(out_stream=None)
        self.nDirections = nDirections

        print 'Check derivatives'
        print 'wrt turbineXw'
        print 'FD: ', self.J[('wtVelocity0', 'turbineXw')]['J_fd']
        print 'FWD: ', self.J[('wtVelocity0', 'turbineXw')]['J_rev']

        print 'wrt turbineYw'
        print 'FD: ', self.J[('wtVelocity0', 'turbineYw')]['J_fd']
        print 'FWD: ', self.J[('wtVelocity0', 'turbineYw')]['J_rev']

        print 'wrt hubHeight'
        print 'FD: ', self.J[('wtVelocity0', 'hubHeight')]['J_fd']
        print 'FWD: ', self.J[('wtVelocity0', 'hubHeight')]['J_rev']

        print 'wrt yaw0'
        print 'FD: ', self.J[('wtVelocity0', 'yaw0')]['J_fd']
        print 'FWD: ', self.J[('wtVelocity0', 'yaw0')]['J_rev']

        print 'wrt rotorDiameter'
        print 'FD: ', self.J[('wtVelocity0', 'rotorDiameter')]['J_fd']
        print 'FWD: ', self.J[('wtVelocity0', 'rotorDiameter')]['J_rev']


    def testObj(self):

        np.testing.assert_allclose(self.J[('wtVelocity0', 'turbineXw')]['J_rev'], self.J[('wtVelocity0', 'turbineXw')]['J_fd'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('wtVelocity0', 'turbineYw')]['J_rev'], self.J[('wtVelocity0', 'turbineYw')]['J_fd'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('wtVelocity0', 'hubHeight')]['J_rev'], self.J[('wtVelocity0', 'hubHeight')]['J_fd'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('wtVelocity0', 'yaw0')]['J_rev'], self.J[('wtVelocity0', 'yaw0')]['J_fd'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('wtVelocity0', 'rotorDiameter')]['J_rev'], self.J[('wtVelocity0', 'rotorDiameter')]['J_fd'], self.rtol, self.atol)


if __name__ == "__main__":
    unittest.main()
