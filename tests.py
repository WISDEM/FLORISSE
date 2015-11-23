import unittest
from GeneralWindfarmComponents import *
from floris_openmdao1 import *
from OptimizationGroups import *
from _floris import *


#
# class GeneralWindFarmComponentsTest(unittest.TestCase):
#
#     def setUp(self):
#
#     def test(self):

class FlorisFortranTest(unittest.TestCase):

    def setUp(self):

        nTurbines = 4
        self.atol = 1E-3
        self.rtol = 1E-6
        turbineX = np.random.rand(nTurbines)*3000.
        turbineY = np.random.rand(nTurbines)*3000.

        # initialize input variable arrays
        rotorDiameter = np.ones(nTurbines)*np.random.random()*150.
        axialInduction = np.ones(nTurbines)*np.random.random()*(1./3.)
        Ct = np.ones(nTurbines)*np.random.random()
        Cp = np.ones(nTurbines)*np.random.random()
        generator_efficiency = np.ones(nTurbines)*np.random.random()*150
        yaw = np.random.rand(nTurbines)*3000.

        # Define flow properties
        wind_speed = np.random.random()*20        # m/s
        air_density = 1.1716    # kg/m^3
        wind_direction = np.random.random()*360    # deg (N = 0 deg., using direction FROM, as in met-mast data)

        # set up problem
        prob = Problem(root=Group())
        prob.root.add('myFloris', FLORIS(nTurbines, resolution=0))
        prob.root.add('CtCp', AdjustCtCpYaw(nTurbines))

        # connect components
        prob.root.connect('CtCp.Ct_out', 'myFloris.Ct')
        prob.root.connect('CtCp.Cp_out', 'myFloris.Cp')
        prob.root.connect('CtCp.yaw', 'myFloris.yaw')

        # initialize problem
        prob.setup()

        # assign values to constant inputs (not design variables)
        prob['myFloris.turbineX'] = turbineX
        prob['myFloris.turbineY'] = turbineY
        prob['myFloris.rotorDiameter'] = rotorDiameter
        prob['myFloris.axialInduction'] = axialInduction
        prob['myFloris.generator_efficiency'] = generator_efficiency
        prob['myFloris.wind_speed'] = wind_speed
        prob['myFloris.air_density'] = air_density
        prob['myFloris.wind_direction'] = wind_direction
        prob['CtCp.yaw'] = yaw
        prob['CtCp.Ct_in'] = Ct
        prob['CtCp.Cp_in'] = Cp
        prob['myFloris.floris_params:FLORISoriginal'] = True
        prob['CtCp.floris_params:FLORISoriginal'] = True

        prob.run()

        self.prob = prob

        self.J = self.prob.check_partial_derivatives()

        # print self.J

    def testFlorisPowerGrads_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)
        # np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)

    def testFlorisPowerGrads_wt_power(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)
        # np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)

    def testFlorisPowerGrads_power(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('power', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('power', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('power', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)
        # np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('power', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('power', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('power', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('power', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)

    def testFlorisCentDiamGrads_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)
        # np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)

    def testFlorisCentDiamGrads_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)
        # np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('wt_power', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('wt_power', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)

    def testFlorisCentDiamGrads_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('power', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('power', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('power', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)
        # np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('power', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('power', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('power', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('power', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('power', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)

    def testFlorisOverlap_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'yaw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'yaw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)

    def testFlorisOverlap_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        # np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)

    def testFlorisOverlap_velocitiesTurbines(self):

        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'turbineXw')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeDiametersT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeCentersYT')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Cp')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'Ct')]['J_fd'], self.atol, self.rtol)
        np.testing.assert_allclose(self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fwd'], self.J['myFloris.f_4'][('velocitiesTurbines', 'wakeOverlapTRel')]['J_fd'], self.atol, self.rtol)



if __name__ == "__main__":
    unittest.main()