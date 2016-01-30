import unittest
# from florisse.GeneralWindFarmComponents import *
from florisse.floris import *
# from florisse.OptimizationGroups import *
from _floris import *

import cPickle as pickle


class GradientTestsFLORIS(unittest.TestCase):

    def setUp(self):

        nTurbines = 4
        self.atol = 1E-3
        self.rtol = 1E-6

        np.random.seed(seed=5)

        turbineX = np.random.rand(nTurbines)*3000.
        turbineY = np.random.rand(nTurbines)*3000.

        # initialize input variable arrays
        rotorDiameter = np.ones(nTurbines)*np.random.random()*150.
        # rotorDiameter = np.ones(nTurbines)*126.4
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        Ct = np.ones(nTurbines)*np.random.random()
        Cp = np.ones(nTurbines)*np.random.random()
        generator_efficiency = np.ones(nTurbines)*np.random.random()
        yaw = np.random.rand(nTurbines)*60. - 30.

        # Define flow properties
        wind_speed = np.random.random()*20        # m/s
        air_density = 1.1716    # kg/m^3
        wind_direction = np.random.random()*360    # deg (N = 0 deg., using direction FROM, as in met-mast data)

        # set up problem
        prob = Problem(root=Group())
        prob.root.add('myFloris', FLORIS(nTurbines, resolution=0))
        prob.root.add('CtCp', AdjustCtCpYaw(nTurbines))
        prob.root.add('v1', IndepVarComp('rotorDiameter', rotorDiameter, units='m'), promotes=['*'])
        prob.root.add('v2', IndepVarComp('yaw', yaw, units='deg'), promotes=['*'])
        prob.root.add('v3', IndepVarComp('axialInduction', axialInduction), promotes=['*'])
        prob.root.add('v4', IndepVarComp('turbineX', turbineX), promotes=['*'])
        prob.root.add('v5', IndepVarComp('turbineY', turbineY), promotes=['*'])
        prob.root.add('v6', IndepVarComp('Ct', Ct), promotes=['*'])
        prob.root.add('v7', IndepVarComp('Cp', Cp), promotes=['*'])

        # connect components
        prob.root.connect('CtCp.Ct_out', 'myFloris.Ct')
        prob.root.connect('CtCp.Cp_out', 'myFloris.Cp')
        prob.root.connect('CtCp.yaw', 'myFloris.yaw')
        prob.root.connect('rotorDiameter', 'myFloris.rotorDiameter')
        prob.root.connect('yaw', 'myFloris.yaw')
        prob.root.connect('axialInduction', 'myFloris.axialInduction')
        prob.root.connect('turbineX', 'myFloris.turbineX')
        prob.root.connect('turbineY', 'myFloris.turbineY')
        prob.root.connect('Ct', 'CtCp.Ct_in')
        prob.root.connect('Cp', 'CtCp.Cp_in')

        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
        prob['myFloris.generator_efficiency'] = generator_efficiency
        prob['myFloris.wind_speed'] = wind_speed
        prob['myFloris.air_density'] = air_density
        prob['myFloris.wind_direction'] = wind_direction
        prob['myFloris.floris_params:FLORISoriginal'] = True
        prob['CtCp.floris_params:FLORISoriginal'] = True

        # run problem
        prob.run()

        # indep_list = ['turbineX', 'turbineY', 'yaw', 'rotorDiameter']
        # unknown_list = ['power0']
        # self.J = prob.calc_gradient(indep_list, unknown_list, return_format='array')
        # print self.J
        # pass results to self for use with unit test
        self.J = prob.check_partial_derivatives(out_stream=None)

        # print self.J
    def testWindFrameGrads_turbineXw(self):

        np.testing.assert_allclose(self.J['myFloris.f_1'][('turbineXw', 'turbineX')]['J_fwd'], self.J['myFloris.f_1'][('turbineXw', 'turbineX')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_1'][('turbineXw', 'turbineY')]['J_fwd'], self.J['myFloris.f_1'][('turbineXw', 'turbineY')]['J_fd'], self.atol, self.rtol)

    def testWindFrameGrads_turbineYw(self):

        np.testing.assert_allclose(self.J['myFloris.f_1'][('turbineYw', 'turbineX')]['J_fwd'], self.J['myFloris.f_1'][('turbineYw', 'turbineX')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_1'][('turbineYw', 'turbineY')]['J_fwd'], self.J['myFloris.f_1'][('turbineYw', 'turbineY')]['J_fd'], self.atol, self.rtol)

    def testFlorisCentDiamGrads_wakeCentersYT(self):
        atol = self.atol
        rtol = 1.0e-0
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeCentersYT', 'yaw')]['J_fwd'], self.J['myFloris.f_2'][('wakeCentersYT', 'yaw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeCentersYT', 'Ct')]['J_fwd'], self.J['myFloris.f_2'][('wakeCentersYT', 'Ct')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeCentersYT', 'turbineXw')]['J_fwd'], self.J['myFloris.f_2'][('wakeCentersYT', 'turbineXw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeCentersYT', 'turbineYw')]['J_fwd'], self.J['myFloris.f_2'][('wakeCentersYT', 'turbineYw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2']['wakeCentersYT', 'rotorDiameter']['J_fwd'], self.J['myFloris.f_2']['wakeCentersYT', 'rotorDiameter']['J_fd'], atol, rtol)
        return

    def testFlorisCentDiamGrads_wakeDiametersT(self):
        atol = self.atol
        rtol = 1.0e-3
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeDiametersT', 'yaw')]['J_fwd'], self.J['myFloris.f_2'][('wakeDiametersT', 'yaw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeDiametersT', 'Ct')]['J_fwd'], self.J['myFloris.f_2'][('wakeDiametersT', 'Ct')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeDiametersT', 'turbineXw')]['J_fwd'], self.J['myFloris.f_2'][('wakeDiametersT', 'turbineXw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeDiametersT', 'turbineYw')]['J_fwd'], self.J['myFloris.f_2'][('wakeDiametersT', 'turbineYw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeDiametersT', 'rotorDiameter')]['J_fwd'], self.J['myFloris.f_2'][('wakeDiametersT', 'rotorDiameter')]['J_fd'], atol, rtol)

    def testFlorisOverlapGrads_wakeOverlapTRel(self):

        np.testing.assert_allclose(self.J['myFloris.f_3'][('wakeOverlapTRel', 'turbineYw')]['J_fwd'], self.J['myFloris.f_3'][('wakeOverlapTRel', 'turbineYw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_3'][('wakeOverlapTRel', 'rotorDiameter')]['J_fwd'], self.J['myFloris.f_3'][('wakeOverlapTRel', 'rotorDiameter')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_3'][('wakeOverlapTRel', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_3'][('wakeOverlapTRel', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_3'][('wakeOverlapTRel', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_3'][('wakeOverlapTRel', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)

    def testFlorisOverlapGrads_cosFac(self):

        np.testing.assert_allclose(self.J['myFloris.f_3'][('cosFac', 'turbineYw')]['J_fwd'], self.J['myFloris.f_3'][('cosFac', 'turbineYw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_3'][('cosFac', 'rotorDiameter')]['J_fwd'], self.J['myFloris.f_3'][('cosFac', 'rotorDiameter')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_3'][('cosFac', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_3'][('cosFac', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_3'][('cosFac', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_3'][('cosFac', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)

    def testFlorisPowerGrads_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines0', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines0', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines0', 'cosFac')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines0', 'cosFac')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines0', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines0', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines0', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines0', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines0', 'axialInduction')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines0', 'axialInduction')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines0', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines0', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines0', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines0', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines0', 'rotorDiameter')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines0', 'rotorDiameter')]['J_fd'], self.atol, self.rtol)
        # return

    def testFlorisPowerGrads_wt_power(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power0', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('wt_power0', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power0', 'cosFac')]['J_fwd'], self.J['myFloris.f_4'][('wt_power0', 'cosFac')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power0', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('wt_power0', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power0', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('wt_power0', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power0', 'axialInduction')]['J_fwd'], self.J['myFloris.f_4'][('wt_power0', 'axialInduction')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power0', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('wt_power0', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power0', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('wt_power0', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power0', 'rotorDiameter')]['J_fwd'], self.J['myFloris.f_4'][('wt_power0', 'rotorDiameter')]['J_fd'], self.atol, self.rtol)

    def testFlorisPowerGrads_power(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('power0', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('power0', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power0', 'cosFac')]['J_fwd'], self.J['myFloris.f_4'][('power0', 'cosFac')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power0', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('power0', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power0', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('power0', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power0', 'axialInduction')]['J_fwd'], self.J['myFloris.f_4'][('power0', 'axialInduction')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power0', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('power0', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power0', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('power0', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4']['power0', 'rotorDiameter']['J_fwd'], self.J['myFloris.f_4']['power0', 'rotorDiameter']['J_fd'], self.atol, self.rtol)


class GradientTestsCtCp(unittest.TestCase):

    def setUp(self):

        nTurbines = 4
        self.atol = 1E-6
        self.rtol = 1E-6

        np.random.seed(seed=10)

        turbineX = np.random.rand(nTurbines)*3000.
        turbineY = np.random.rand(nTurbines)*3000.

        # initialize input variable arrays
        rotorDiameter = np.ones(nTurbines)*np.random.random()*150.
        # rotorDiameter = np.ones(nTurbines)*126.4
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        Ct = np.ones(nTurbines)*np.random.random()
        Cp = np.ones(nTurbines)*np.random.random()
        generator_efficiency = np.ones(nTurbines)*np.random.random()
        yaw = np.random.rand(nTurbines)*60. - 30.

        # Define flow properties
        wind_speed = np.random.random()*20        # m/s
        air_density = 1.1716    # kg/m^3
        wind_direction = np.random.random()*360    # deg (N = 0 deg., using direction FROM, as in met-mast data)

        # set up problem
        prob = Problem(root=DirectionGroupFLORIS(nTurbines=nTurbines))
        prob.root.add('v1', IndepVarComp('rotorDiameter', rotorDiameter, units='m'), promotes=['*'])
        prob.root.add('v2', IndepVarComp('yaw', yaw, units='deg'), promotes=['*'])
        prob.root.add('v3', IndepVarComp('axialInduction', axialInduction), promotes=['*'])
        prob.root.add('v4', IndepVarComp('turbineX', turbineX), promotes=['*'])
        prob.root.add('v5', IndepVarComp('turbineY', turbineY), promotes=['*'])
        prob.root.add('v6', IndepVarComp('Ct_in', np.zeros(nTurbines)), promotes=['*'])
        prob.root.add('v7', IndepVarComp('Cp_in', np.zeros(nTurbines)), promotes=['*'])


        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
        prob['generator_efficiency'] = generator_efficiency
        prob['wind_speed'] = wind_speed
        prob['air_density'] = air_density
        prob['wind_direction'] = wind_direction
        prob['floris_params:FLORISoriginal'] = False

        # run problem
        prob.run()

        # pass gradient test results to self for use with unit tests
        self.J = prob.check_partial_derivatives(out_stream=None)

    def testCtCp_Ct_out(self):
        np.testing.assert_allclose(self.J['CtCp'][('Ct_out', 'Ct_in')]['J_fwd'], self.J['CtCp'][('Ct_out', 'Ct_in')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['CtCp'][('Ct_out', 'Cp_in')]['J_fwd'], self.J['CtCp'][('Ct_out', 'Cp_in')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['CtCp'][('Ct_out', 'yaw')]['J_fwd'], self.J['CtCp'][('Ct_out', 'yaw')]['J_fd'], self.atol, self.rtol)

    def testCtCp_Cp_out(self):
        np.testing.assert_allclose(self.J['CtCp'][('Cp_out', 'Ct_in')]['J_fwd'], self.J['CtCp'][('Cp_out', 'Ct_in')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['CtCp'][('Cp_out', 'Cp_in')]['J_fwd'], self.J['CtCp'][('Cp_out', 'Cp_in')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['CtCp'][('Cp_out', 'yaw')]['J_fwd'], self.J['CtCp'][('Cp_out', 'yaw')]['J_fd'], self.atol, self.rtol)


class GradientTestsCtCpRotor(unittest.TestCase):

    def setUp(self):

        nTurbines = 4
        self.atol = 1E-6
        self.rtol = 1E-6

        np.random.seed(seed=10)

        turbineX = np.random.rand(nTurbines)*3000.
        turbineY = np.random.rand(nTurbines)*3000.

        # initialize input variable arrays
        rotorDiameter = np.ones(nTurbines)*np.random.random()*150.
        # rotorDiameter = np.ones(nTurbines)*126.4
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        Ct = np.ones(nTurbines)*np.random.random()
        Cp = np.ones(nTurbines)*np.random.random()
        generator_efficiency = np.ones(nTurbines)*np.random.random()
        yaw = np.random.rand(nTurbines)*60. - 30.

        # Define flow properties
        wind_speed = np.random.random()*20        # m/s
        air_density = 1.1716    # kg/m^3
        wind_direction = np.random.random()*360    # deg (N = 0 deg., using direction FROM, as in met-mast data)

        NREL5MWCPCT = pickle.load(open('../NREL5MWCPCT_dict.p'))
        datasize = NREL5MWCPCT['CP'].size

        # set up problem
        prob = Problem(root=DirectionGroupFLORIS(nTurbines=nTurbines, use_rotor_components=True, datasize=datasize))
        prob.root.add('v1', IndepVarComp('rotorDiameter', rotorDiameter, units='m'), promotes=['*'])
        prob.root.add('v2', IndepVarComp('yaw', yaw, units='deg'), promotes=['*'])
        prob.root.add('v3', IndepVarComp('axialInduction', axialInduction), promotes=['*'])
        prob.root.add('v4', IndepVarComp('turbineX', turbineX), promotes=['*'])
        prob.root.add('v5', IndepVarComp('turbineY', turbineY), promotes=['*'])

        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
        prob['generator_efficiency'] = generator_efficiency
        prob['wind_speed'] = wind_speed
        prob['air_density'] = air_density
        prob['wind_direction'] = wind_direction
        prob['floris_params:FLORISoriginal'] = True
        prob['floris_params:FLORISoriginal'] = True

        # values for rotor coupling
        prob['params:windSpeedToCPCT:CP'] = NREL5MWCPCT['CP']
        prob['params:windSpeedToCPCT:CT'] = NREL5MWCPCT['CT']
        prob['params:windSpeedToCPCT:wind_speed'] = NREL5MWCPCT['wind_speed']
        prob['floris_params:ke'] = 0.05
        prob['floris_params:kd'] = 0.17
        prob['floris_params:aU'] = 12.0
        prob['floris_params:bU'] = 1.3
        prob['floris_params:initialWakeAngle'] = 3.0
        prob['floris_params:useaUbU'] = True
        prob['floris_params:useWakeAngle'] = True
        prob['floris_params:adjustInitialWakeDiamToYaw'] = False

        # run problem
        prob.run()

        # pass gradient test results to self for use with unit tests
        self.J = prob.check_partial_derivatives(out_stream=None)

    def testCtCpRotor_Cp_out(self):
        np.testing.assert_allclose(self.J['CtCp'][('Cp_out', 'yaw')]['J_fwd'], self.J['CtCp'][('Cp_out', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['CtCp'][('Cp_out', 'velocitiesTurbines0')]['J_fwd'], self.J['CtCp'][('Cp_out', 'velocitiesTurbines0')]['J_fd'], self.atol, self.rtol)

    def testCtCpRotor_Ct_out(self):
        np.testing.assert_allclose(self.J['CtCp'][('Ct_out', 'yaw')]['J_fwd'], self.J['CtCp'][('Ct_out', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['CtCp'][('Ct_out', 'velocitiesTurbines0')]['J_fwd'], self.J['CtCp'][('Ct_out', 'velocitiesTurbines0')]['J_fd'], self.atol, self.rtol)


if __name__ == "__main__":
    unittest.main()