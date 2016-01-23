import unittest
from GeneralWindfarmComponents import *
from floris_openmdao1 import *
# from OptimizationGroups import *
from _floris import *


#
# class GeneralWindFarmComponentsTest(unittest.TestCase):
#
#     def setUp(self):
#
#     def test(self):
# class MPITests(unittest.TestCase):
#
#     def setUp(self):
#


class GradientTests(unittest.TestCase):

    def setUp(self):

        nTurbines = 4
        self.atol = 1E-3
        self.rtol = 1E-6
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
        atol = 1.0e-6
        rtol = 1.0e-0
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeCentersYT', 'yaw')]['J_fwd'], self.J['myFloris.f_2'][('wakeCentersYT', 'yaw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeCentersYT', 'Ct')]['J_fwd'], self.J['myFloris.f_2'][('wakeCentersYT', 'Ct')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeCentersYT', 'turbineXw')]['J_fwd'], self.J['myFloris.f_2'][('wakeCentersYT', 'turbineXw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2'][('wakeCentersYT', 'turbineYw')]['J_fwd'], self.J['myFloris.f_2'][('wakeCentersYT', 'turbineYw')]['J_fd'], atol, rtol)
        np.testing.assert_allclose(self.J['myFloris.f_2']['wakeCentersYT', 'rotorDiameter']['J_fwd'], self.J['myFloris.f_2']['wakeCentersYT', 'rotorDiameter']['J_fd'], atol, rtol)
        return

    def testFlorisCentDiamGrads_wakeDiametersT(self):
        atol = 1.0e-6
        rtol = 1.0e-0
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

if __name__ == "__main__":
    unittest.main()