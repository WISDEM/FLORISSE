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
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        Ct = np.ones(nTurbines)*np.random.random()
        Cp = np.ones(nTurbines)*np.random.random()
        generator_efficiency = np.ones(nTurbines)*np.random.random()
        yaw = np.random.rand(nTurbines)*60. - 30.

        # Define flow properties
        wind_speed = np.random.random()*20        # m/s
        air_density = 1.1716    # kg/m^3
        wind_direction = np.random.random()*360    # deg (N = 0 deg., using direction FROM, as in met-mast data)
        wind_frequency = np.random.random()    # probability of wind in given direction

        # set up problem
        prob = Problem(root=AEPGroupFLORIS(nTurbines, nDirections=1, resolution=0))

        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
        prob['turbineX'] = turbineX
        prob['turbineY'] = turbineY
        prob['yaw0'] = yaw
        prob['rotorDiameter'] = rotorDiameter
        prob['axialInduction'] = axialInduction
        prob['Ct_in'] = Ct
        prob['Cp_in'] = Cp
        prob['generator_efficiency'] = generator_efficiency
        prob['windSpeeds'] = np.array([wind_speed])
        prob['air_density'] = air_density
        prob['windDirections'] = np.array([wind_direction])
        prob['windrose_frequencies'] = np.array([wind_frequency])
        prob['floris_params:FLORISoriginal'] = True

        # run problem
        prob.run()

        # pass results to self for use with unit test
        self.J = prob.check_partial_derivatives(out_stream=None)

        # print self.J

    def testWindFrameGrads_turbineXw(self):

        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_1'][('turbineXw', 'turbineX')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_1'][('turbineXw', 'turbineX')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_1'][('turbineXw', 'turbineY')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_1'][('turbineXw', 'turbineY')]['J_fd'], self.atol, self.rtol)

    def testWindFrameGrads_turbineYw(self):

        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_1'][('turbineYw', 'turbineX')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_1'][('turbineYw', 'turbineX')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_1'][('turbineYw', 'turbineY')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_1'][('turbineYw', 'turbineY')]['J_fd'], self.atol, self.rtol)

    def testFlorisCentDiamGrads_wakeCentersYT(self):
        atol = self.atol
        rtol = 1.0e-0
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeCentersYT', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeCentersYT', 'yaw0')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeCentersYT', 'Ct')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeCentersYT', 'Ct')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeCentersYT', 'turbineXw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeCentersYT', 'turbineXw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeCentersYT', 'turbineYw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeCentersYT', 'turbineYw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2']['wakeCentersYT', 'rotorDiameter']['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2']['wakeCentersYT', 'rotorDiameter']['J_fd'], atol, rtol)
        return

    def testFlorisCentDiamGrads_wakeDiametersT(self):
        atol = self.atol
        rtol = 1.0e-3
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'yaw0')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'Ct')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'Ct')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'turbineXw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'turbineXw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'turbineYw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'turbineYw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'rotorDiameter')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_2'][('wakeDiametersT', 'rotorDiameter')]['J_fd'], atol, rtol)

    def testFlorisOverlapGrads_wakeOverlapTRel(self):

        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_3'][('wakeOverlapTRel', 'turbineYw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_3'][('wakeOverlapTRel', 'turbineYw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_3'][('wakeOverlapTRel', 'rotorDiameter')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_3'][('wakeOverlapTRel', 'rotorDiameter')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_3'][('wakeOverlapTRel', 'wakeDiametersT')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_3'][('wakeOverlapTRel', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_3'][('wakeOverlapTRel', 'wakeCentersYT')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_3'][('wakeOverlapTRel', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)

    def testFlorisOverlapGrads_cosFac(self):

        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_3'][('cosFac', 'turbineYw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_3'][('cosFac', 'turbineYw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_3'][('cosFac', 'rotorDiameter')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_3'][('cosFac', 'rotorDiameter')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_3'][('cosFac', 'wakeDiametersT')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_3'][('cosFac', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_3'][('cosFac', 'wakeCentersYT')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_3'][('cosFac', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)

    def testFlorisPowerGrads_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'wakeOverlapTRel')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'cosFac')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'cosFac')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'Ct')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'Cp')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'axialInduction')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'axialInduction')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'turbineXw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'yaw0')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'rotorDiameter')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('velocitiesTurbines0', 'rotorDiameter')]['J_fd'], self.atol, self.rtol)
        # return

    def testFlorisPowerGrads_wt_power(self):

        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'wakeOverlapTRel')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'cosFac')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'cosFac')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'Ct')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'Cp')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'axialInduction')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'axialInduction')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'turbineXw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'yaw0')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'rotorDiameter')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('wt_power0', 'rotorDiameter')]['J_fd'], self.atol, self.rtol)

    def testFlorisPowerGrads_power(self):

        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'wakeOverlapTRel')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'cosFac')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'cosFac')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'Ct')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'Cp')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'axialInduction')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'axialInduction')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'turbineXw')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4'][('power0', 'yaw0')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.myFloris.f_4']['power0', 'rotorDiameter']['J_fwd'], self.J['all_directions.direction_group0.myFloris.f_4']['power0', 'rotorDiameter']['J_fd'], self.atol, self.rtol)


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
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        Ct = np.ones(nTurbines)*np.random.random()
        Cp = np.ones(nTurbines)*np.random.random()
        generator_efficiency = np.ones(nTurbines)*np.random.random()
        yaw = np.random.rand(nTurbines)*60. - 30.

        # Define flow properties
        wind_speed = np.random.random()*20        # m/s
        air_density = 1.1716    # kg/m^3
        wind_direction = np.random.random()*360    # deg (N = 0 deg., using direction FROM, as in met-mast data)
        wind_frequency = np.random.random()    # probability of wind in given direction

        # set up problem
        prob = Problem(root=AEPGroupFLORIS(nTurbines=nTurbines))

        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
                # assign values to constant inputs (not design variables)
        prob['turbineX'] = turbineX
        prob['turbineY'] = turbineY
        prob['yaw0'] = yaw
        prob['rotorDiameter'] = rotorDiameter
        prob['axialInduction'] = axialInduction
        prob['Ct_in'] = Ct
        prob['Cp_in'] = Cp
        prob['generator_efficiency'] = generator_efficiency
        prob['windSpeeds'] = np.array([wind_speed])
        prob['windrose_frequencies'] = np.array([wind_frequency])
        prob['air_density'] = air_density
        prob['windDirections'] = np.array([wind_direction])
        prob['floris_params:FLORISoriginal'] = False

        # run problem
        prob.run()

        # pass gradient test results to self for use with unit tests
        self.J = prob.check_partial_derivatives(out_stream=None)

    def testCtCp_Ct_out(self):
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'Ct_in')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'Ct_in')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'Cp_in')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'Cp_in')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'yaw0')]['J_fd'], self.atol, self.rtol)

    def testCtCp_Cp_out(self):
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'Ct_in')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'Ct_in')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'Cp_in')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'Cp_in')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'yaw0')]['J_fd'], self.atol, self.rtol)


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
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        generator_efficiency = np.ones(nTurbines)*np.random.random()
        yaw = np.random.rand(nTurbines)*60. - 30.

        # Define flow properties
        wind_speed = np.random.random()*20        # m/s
        air_density = 1.1716    # kg/m^3
        wind_direction = np.random.random()*360    # deg (N = 0 deg., using direction FROM, as in met-mast data)
        wind_frequency = np.random.random()    # probability of wind in given direction



        NREL5MWCPCT = pickle.load(open('NREL5MWCPCT_dict.p'))
        datasize = NREL5MWCPCT['CP'].size

        # set up problem
        prob = Problem(root=AEPGroupFLORIS(nTurbines=nTurbines, use_rotor_components=True, datasize=datasize))

        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
        prob['turbineX'] = turbineX
        prob['turbineY'] = turbineY
        prob['yaw0'] = yaw
        prob['rotorDiameter'] = rotorDiameter
        prob['axialInduction'] = axialInduction
        prob['generator_efficiency'] = generator_efficiency
        prob['windSpeeds'] = np.array([wind_speed])
        prob['windrose_frequencies'] = np.array([wind_frequency])
        prob['air_density'] = air_density
        prob['windDirections'] = np.array([wind_direction])
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
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'yaw0')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'velocitiesTurbines0')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Cp_out', 'velocitiesTurbines0')]['J_fd'], self.atol, self.rtol)

    def testCtCpRotor_Ct_out(self):
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'yaw0')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'yaw0')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'velocitiesTurbines0')]['J_fwd'], self.J['all_directions.direction_group0.CtCp'][('Ct_out', 'velocitiesTurbines0')]['J_fd'], self.atol, self.rtol)


if __name__ == "__main__":
    unittest.main()


# indep_list = ['turbineX', 'turbineY', 'yaw', 'rotorDiameter']
# unknown_list = ['power0']
# self.J = prob.calc_gradient(indep_list, unknown_list, return_format='array')
# print self.J