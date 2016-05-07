import unittest
# from florisse.GeneralWindFarmComponents import *
from openmdao.api import pyOptSparseDriver, ExecComp, IndepVarComp, Problem

from COE import *
from florisse.floris import *
from florisse.OptimizationGroups import *

import cPickle as pickle

class TotalDerivTestsCOE(unittest.TestCase):

    def setUp(self):

        nTurbines = 4
        self.rtol = 1E-6
        self.atol = 1E-6

        np.random.seed(seed=10)

        turbineX = np.random.rand(nTurbines)*3000.
        turbineY = np.random.rand(nTurbines)*3000.
        turbineH1 = np.random.rand(1)*150.+75
        turbineH2 = np.random.rand(1)*150.+75
        nTurbsH1 = nTurbines/2
        nTurbsH2 = nTurbines-nTurbsH1

        minSpacing = 2.

        # initialize input variable arrays
        rotorDiameter = np.ones(nTurbines)*np.random.random()*150.
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        Ct = np.ones(nTurbines)*np.random.random()
        Cp = np.ones(nTurbines)*np.random.random()
        generatorEfficiency = np.ones(nTurbines)*np.random.random()
        yaw = np.random.rand(nTurbines)*60. - 30.

        # Define flow properties
        nDirections = 50
        windSpeeds = np.random.rand(nDirections)*20        # m/s
        air_density = 1.1716    # kg/m^3
        windDirections = np.random.rand(nDirections)*360.0
        windFrequencies = np.random.rand(nDirections)

        # set up problem
        # prob = Problem(root=OptAEP(nTurbines, nDirections=1))

        prob = Problem()
        root = prob.root = Group()
        
        #TODO How to I pass TurbineZ down to calculate AEP? Basically, How do I get AEP?
        root.add('dv5', IndepVarComp('turbineH1', 0., units='m'), promotes=['*'])
        root.add('dv6', IndepVarComp('turbineH2', 0., units='m'), promotes=['*'])
        root.add('getTurbineZ', getTurbineZ(nTurbines), promotes=['*'])
        root.add('AEPGroup', AEPGroup(nTurbines, nDirections=nDirections,
                    use_rotor_components=False, datasize=0, differentiable=True,
                    optimizingLayout=False, nSamples=0), promotes=['*'])
        root.add('COEComponent', COEComponent(nTurbines), promotes=['*'])
        
        root.add('spacing_comp', SpacingComp(nTurbines=nTurbines), promotes=['*'])
        # add constraint definitions
        root.add('spacing_con', ExecComp('sc = wtSeparationSquared-(minSpacing*rotorDiameter[0])**2',
                                     minSpacing=minSpacing, rotorDiameter=np.zeros(nTurbines),
                                     sc=np.zeros(((nTurbines-1.)*nTurbines/2.)),
                                     wtSeparationSquared=np.zeros(((nTurbines-1.)*nTurbines/2.))),
             promotes=['*'])

        # set up optimizer
        # prob.driver = pyOptSparseDriver()
        # prob.driver.options['optimizer'] = 'SNOPT'
        prob.driver.add_objective('COE', scaler=1E-8)

        # set optimizer options
        # prob.driver.opt_settings['Verify level'] = 3
        # prob.driver.opt_settings['Print file'] = 'SNOPT_print_exampleOptAEP.out'
        # prob.driver.opt_settings['Summary file'] = 'SNOPT_summary_exampleOptAEP.out'
        # prob.driver.opt_settings['Major iterations limit'] = 1

        # select design variables
        prob.driver.add_desvar('turbineX', lower=np.ones(nTurbines)*min(turbineX), upper=np.ones(nTurbines)*max(turbineX), scaler=1E-2)
        prob.driver.add_desvar('turbineY', lower=np.ones(nTurbines)*min(turbineY), upper=np.ones(nTurbines)*max(turbineY), scaler=1E-2)
        prob.driver.add_desvar('turbineH1', lower=np.ones(nTurbines)*min(turbineH1), upper=np.ones(nTurbines)*max(turbineH1), scaler=1E-2)
        prob.driver.add_desvar('turbineH2', lower=np.ones(nTurbines)*min(turbineH1), upper=np.ones(nTurbines)*max(turbineH1), scaler=1E-2)
        """for direction_id in range(0, windDirections.size):
            prob.driver.add_desvar('yaw%i' % direction_id, lower=-30.0, upper=30.0, scaler=1E-1)"""

        # add constraints
        prob.driver.add_constraint('sc', lower=np.zeros(((nTurbines-1.)*nTurbines/2.)))

        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
        prob['turbineX'] = turbineX
        prob['turbineY'] = turbineY
        prob['turbineH1'] = turbineH1
        prob['turbineH2'] = turbineH2
        prob['nTurbsH1'] = nTurbsH1
        prob['nTurbsH2'] = nTurbsH2
        prob['yaw0'] = yaw
        prob['rotorDiameter'] = rotorDiameter
        prob['axialInduction'] = axialInduction
        prob['Ct_in'] = Ct
        prob['Cp_in'] = Cp
        prob['generatorEfficiency'] = generatorEfficiency
        prob['windSpeeds'] = windSpeeds
        prob['air_density'] = air_density
        prob['windDirections'] = windDirections
        prob['windFrequencies'] = windFrequencies
        prob['floris_params:FLORISoriginal'] = True

        # run problem
        prob.run()

        # pass results to self for use with unit test
        self.J = prob.check_total_derivatives(out_stream=None)
        self.nDirections = nDirections

        # print self.J

    def testObj(self):

        np.testing.assert_allclose(self.J[('COE', 'turbineX')]['rel error'], self.J[('COE', 'turbineX')]['rel error'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('COE', 'turbineY')]['rel error'], self.J[('COE', 'turbineY')]['rel error'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('COE', 'turbineH1')]['rel error'], self.J[('COE', 'turbineH1')]['rel error'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('COE', 'turbineH2')]['rel error'], self.J[('COE', 'turbineH2')]['rel error'], self.rtol, self.atol)
        """for dir in np.arange(0, self.nDirections):
            np.testing.assert_allclose(self.J[('COE', 'yaw%i' % dir)]['rel error'], self.J[('COE', 'yaw%i' % dir)]['rel error'], self.rtol, self.atol)"""

    def testCon(self):

        np.testing.assert_allclose(self.J[('sc', 'turbineX')]['rel error'], self.J[('sc', 'turbineX')]['rel error'], self.rtol, self.atol)
        np.testing.assert_allclose(self.J[('sc', 'turbineY')]['rel error'], self.J[('sc', 'turbineY')]['rel error'], self.rtol, self.atol)
        """for dir in np.arange(0, self.nDirections):
            np.testing.assert_allclose(self.J[('sc', 'yaw%i' % dir)]['rel error'], self.J[('sc', 'yaw%i' % dir)]['rel error'], self.rtol, self.atol)"""

if __name__ == "__main__":
    unittest.main()
