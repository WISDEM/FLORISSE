from openmdao.main.api import Component, VariableTree
from openmdao.lib.datatypes.api import Array, Bool, Float, VarTree
from Parameters import FLORISParameters
import numpy as np
import _floris


class floris_wcent_wdiam(Component):
    """ Calculates the center and diameter of each turbine wake at each other turbine """

    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')
    turbineXw = Array(iotype='in', desc='x coordinates of turbines in wind dir. ref. frame')
    turbineYw = Array(iotype='in', desc='y coordinates of turbines in wind dir. ref. frame')
    yaw = Array(iotype='in', desc='yaw of each turbine')
    rotorDiameter = Array(dtype='float', iotype='in', desc='rotor diameter of each turbine')
    Ct = Array(iotype='in', dtype='float', desc='thrust coefficient of each turbine')

    wakeCentersYT = Array(iotype='out', dtype='float', desc='wake center y position at each turbine')
    wakeDiametersT = Array(iotype='out', dtype='float', desc='wake diameter of each zone of each wake at each turbine')
    # p_near0 = Float(iotyp='out', desc='upwind location of diameter spline in rotor diameters')

    def execute(self):

        # print 'entering wcent_wdiam - tapenade'

        # rename inputs and outputs
        # pP = self.parameters.pP
        kd = self.parameters.kd
        ke = self.parameters.ke
        initialWakeDisplacement = self.parameters.initialWakeDisplacement
        initialWakeAngle = self.parameters.initialWakeAngle
        rotorDiameter = self.rotorDiameter
        Ct = self.Ct
        keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        me = self.parameters.me

        turbineXw = self.turbineXw
        turbineYw = self.turbineYw

        yaw_deg = self.yaw
        # yaw = self.yaw*np.pi/180.0
        # print 'yaw = ', yaw_deg

        wakeCentersYT_vec, wakeDiametersT_vec = _floris.floris_wcent_wdiam(kd, initialWakeDisplacement, \
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me)



        self.wakeCentersYT = wakeCentersYT_vec
        self.wakeDiametersT = wakeDiametersT_vec

    def list_deriv_vars(self):

        return ('yaw', 'Ct', 'turbineXw', 'turbineYw', 'rotorDiameter'), ('wakeCentersYT', 'wakeDiametersT')

    def provideJ(self):

        kd = self.parameters.kd
        ke = self.parameters.ke
        initialWakeDisplacement = self.parameters.initialWakeDisplacement
        initialWakeAngle = self.parameters.initialWakeAngle
        rotorDiameter = self.rotorDiameter
        Ct = self.Ct
        keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        me = self.parameters.me

        turbineXw = self.turbineXw
        turbineYw = self.turbineYw
        # print self.turbineXw
        yaw_deg = self.yaw
        # yaw = self.yaw*np.pi/180.0

        nTurbines = np.size(turbineXw)
        nbdirs = 3*nTurbines*nTurbines

        wakeDiametersT_vecb = np.zeros((nbdirs, 3*nTurbines*nTurbines))

        wakeCentersYT_vecb = np.eye(nbdirs, nTurbines*nTurbines)

        # print nTurbines, turbineXw, me, wakeCentersYT_vecb

        yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb, _, _ = _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement, \
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me, wakeCentersYT_vecb, wakeDiametersT_vecb)

        dwc = np.hstack((yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb))
        dwc = dwc[:nTurbines*nTurbines, :]

        wakeCentersYT_vecb[:, :] = 0.0

        wakeDiametersT_vecb = np.eye(nbdirs, nbdirs)

        yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb, _, _ = _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement, \
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me, wakeCentersYT_vecb, wakeDiametersT_vecb)


        dwd = np.hstack((yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb))

        J = np.vstack((dwc, dwd))

        # print J.shape

        return J


class floris_overlap(Component):
    """ Calculates the overlap between each turbine rotor and the existing turbine wakes """
    turbineXw = Array(iotype='in', units='m', deriv_ignore=True, desc='X positions of turbines wrt the wind direction')
    turbineYw = Array(iotype='in', units='m', desc='Y positions of turbines wrt the wind direction')
    rotorDiameter = Array(iotype='in', units='m', desc='diameters of all turbine rotors')
    wakeDiametersT = Array(iotype='in', units='m', desc='diameters of all turbines wake zones')
    wakeCentersYT = Array(iotype='in', units='m', desc='Y positions of all wakes at each turbine')
    rotorArea = Array(iotype='in', units='m*m', desc='Area of each turbine rotor')

    wakeOverlapTRel = Array(iotype='out', desc='relative wake zone overlap to rotor area')

    def execute(self):
        # print 'entering overlap - tapenade'

        nTurbines = self.turbineYw.size
        # p_near0 = self.p_near0
        # calculate overlap areas at rotors
        # wakeOverlapT(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake
        # of turbine TURB with rotor of turbine TURBI
        # wakeOverlapT = calcOverlapAreas(self.turbineXw, self.turbineYw, self.rotorDiameter, self.wakeDiametersT, self.wakeCentersYT, p_near0)

        # make overlap relative to rotor area (maximum value should be 1)
        # wakeOverlapTRel = wakeOverlapT
        # for turb in range(0, nTurbines): # Jared: I think it would make more sense to use turbI for consistency
        #     wakeOverlapTRel[turb] = wakeOverlapTRel[turb]/self.rotorArea[turb]

        # print 'wakeDiametersT = ', self.wakeDiametersT
        # print 'wakeCentersYT = ', self.wakeCentersYT
        # print 'turbineYw = ', self.turbineYw

        # wakeOverlapTRel = calcOverlapAreas(self.turbineXw, self.turbineYw, self.rotorDiameter, self.wakeDiametersT, self.wakeCentersYT, p_near0)
        wakeOverlapTRel_vec = _floris.floris_overlap(self.turbineXw, self.turbineYw, self.rotorDiameter, \
                                                     self.wakeDiametersT, self.wakeCentersYT)

        # self.wakeOverlapTRel = wakeOverlapTRel
        # print self.wakeOverlapTRel
        # print '_'

        # self.wakeOverlapTRel = wakeOverlapTRel
        self.wakeOverlapTRel = wakeOverlapTRel_vec

    def list_deriv_vars(self):
        """specifies the inputs and outputs where derivatives are defined"""

        # return ('turbineYw', 'rotorDiameter', 'wakeDiametersT', 'wakeCentersYT'), 'wakeOverlapTRel'
        return ('turbineYw', 'rotorDiameter', 'wakeDiametersT', 'wakeCentersYT'), ('wakeOverlapTRel',)

    def provideJ(self):

        turbineXw = self.turbineXw
        turbineYw = self.turbineYw
        rotorDiameter = self.rotorDiameter
        wakeDiametersT_vec = self.wakeDiametersT
        wakeCentersYT_vec = self.wakeCentersYT
        wakeOverlapTRel_vec = self.wakeOverlapTRel

        nTurbines = np.size(turbineXw)
        nbdirs = 3*nTurbines*nTurbines

        wakeOverlapTRel_vecb = np.eye(nbdirs, 3*nTurbines*nTurbines)

        turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb \
            = _floris.floris_overlap_bv(turbineXw, turbineYw, rotorDiameter, wakeDiametersT_vec, \
                                         wakeCentersYT_vec, wakeOverlapTRel_vec,  wakeOverlapTRel_vecb)

        J = np.hstack((turbineYwb[:, :], rotorDiameterb[:, :], wakeDiametersT_vecb[:, :], \
                             wakeCentersYT_vecb[:, :]))

        return J


class floris_power(Component):
    """ Calculates the turbine power and effective wind speed for each turbine """

    # original variables in Pieter's OpenMDAO stand-alone version of FLORIS
    parameters = VarTree(FLORISParameters(), iotype='in')
    velocitiesTurbines = Array(iotype='out', units='m/s')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    # input variables added so I don't have to use WISDEM while developing gradients
    rotorDiameter = Array(dtype='float', iotype='in', units='m', desc='rotor diameters of all turbine')
    rotorArea = Array(iotype='in', dtype='float', units='m*m', desc='rotor area of all turbines')
    axialInduction = Array(iotype='in', dtype='float', desc='axial induction of all turbines')
    Ct = Array(iotype='in', dtype='float', desc='Thrust coefficient for all turbines')
    Cp = Array(iotype='in', dtype='float', desc='power coefficient for all turbines')
    generator_efficiency = Array(iotype='in', dtype='float', desc='generator efficiency of all turbines')
    # yaw = Array(iotype='in', missing_deriv_policy='assume_zero', desc='yaw of each turbine')
    turbineXw = Array(iotype='in', dtype='float', units='m', desc='X positions of turbines in the wind direction reference frame')
    wakeCentersYT = Array(iotype='in', units='m', desc='centers of the wakes at each turbine')
    wakeDiametersT = Array(iotype='in', units='m', desc='diameters of each of the wake zones for each of the wakes at each turbine')
    wakeOverlapTRel = Array(iotype='in', units='m', desc='ratios of wake overlap area per zone to rotor area')

    # Flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    air_density = Float(iotype='in', units='kg/(m*m*m)', desc='air density in free stream')

    # output variables added so I don't have to use WISDEM while developing gradients
    wt_power = Array(iotype='out', units='kW')
    power = Float(iotype='out', units='kW', desc='total power output of the wind farm')

    def execute(self):
        # print 'entering power - tapenade'

        wakeOverlapTRel_v = self.wakeOverlapTRel

        ke = self.parameters.ke
        keCorrArray = self.parameters.keCorrArray
        keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        Ct = self.Ct
        Vinf = self.wind_speed
        turbineXw = self.turbineXw
        axialInduction = self.axialInduction
        rotorDiameter = self.rotorDiameter
        rotorArea = self.rotorArea
        rho = self.air_density
        generator_efficiency = self.generator_efficiency
        Cp = self.Cp
        MU = self.parameters.MU
        axialIndProvided = self.parameters.axialIndProvided

        p_near0 = 1.0
        self.p_near0 = p_near0

        # print 'turbineXw = ', turbineXw
        # print 'wakeOverlapTRel = ', wakeOverlapTRel_v

        # print 'nTurbines = ', nTurbines
        # print 'wakeOverlapTRel  = ', wakeOverlapTRel_v
        # print 'Ct  = ', Ct
        # print 'axialInduction  = ', axialInduction
        # # print 'yaw  = ', yaw
        # print 'axialIndProvided type = ', type(axialIndProvided)
        # print 'keCorrCT type = ', type(keCorrCT)
        # print 'Region2CT type = ', type(Region2CT)
        # print 'ke type = ', type(ke)
        # print 'Vinf type = ', type(Vinf)
        # print 'keCorrArray type = ', type(keCorrArray)
        # print 'turbineXw  = ', turbineXw
        # print 'p_near0 type = ', type(p_near0)
        # print 'rotorDiameter  = ', rotorDiameter
        # print 'MU  = ', MU
        # print 'rho type = ', type(rho)
        # print 'Cp  = ', Cp
        # print 'generator_efficiency  = ', generator_efficiency

        # print 'nTurbines = ', nTurbines
        # print 'wakeOverlapTRel shape  = ', wakeOverlapTRel_v.shape
        # print 'Ct  = shape ', Ct.shape
        # print 'axialInduction  shape = ', axialInduction.shape
        # print 'yaw  = ', yaw
        # print 'axialIndProvided type = ', type(axialIndProvided)
        # print 'keCorrCT type = ', type(keCorrCT)
        # print 'Region2CT type = ', type(Region2CT)
        # print 'ke type = ', type(ke)
        # print 'Vinf type = ', type(Vinf)
        # print 'keCorrArray type = ', type(keCorrArray)
        # print 'turbineXw  shape = ', turbineXw.shape
        # print 'p_near0 type = ', type(p_near0)
        # print 'rotorDiameter shape = ', rotorDiameter.shape
        # print 'MU shape = ', MU.shape
        # print 'rho type = ', type(rho)
        # print 'Cp  shape = ', Cp.shape
        # print 'generator_efficiency  = ', generator_efficiency

        velocitiesTurbines, wt_power, power = _floris.floris_power(wakeOverlapTRel_v, Ct, axialInduction, \
                                                            axialIndProvided, keCorrCT, Region2CT, ke, \
                                                            Vinf, keCorrArray, turbineXw, p_near0, rotorDiameter, MU, \
                                                            rho, Cp, generator_efficiency)

        # print 'fortran call complete'
        # print velocitiesTurbines, wt_power, power
        if self.verbose:
            print "wind speed at turbines %s [m/s]" % velocitiesTurbines
            print "rotor area %s" % rotorArea
            print "rho %s" % rho
            print "generator_efficiency %s" % generator_efficiency
            print "powers turbines %s [kW]" % wt_power

        self.velocitiesTurbines = velocitiesTurbines
        self.wt_power = wt_power
        self.power = power
        # print power



    def list_deriv_vars(self):
        """specifies the inputs and outputs where derivatives are defined"""

        return ('wakeOverlapTRel', 'Ct', 'axialInduction', 'turbineXw', \
                'rotorDiameter', 'Cp'), ('velocitiesTurbines', 'wt_power', 'power')
        # return ('turbineX', 'turbineY'), ('turbineXw', 'turbineYw', 'wsw_position')

    def provideJ(self):

        nTurbines = np.size(self.turbineXw)
        nbdirs = nTurbines
        # nbdirs = 3*nTurbines*nTurbines

        wakeOverlapTRel_v = self.wakeOverlapTRel

        ke = self.parameters.ke
        keCorrArray = self.parameters.keCorrArray
        keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        Ct = self.Ct
        Vinf = self.wind_speed
        turbineXw = self.turbineXw
        axialInduction = self.axialInduction
        rotorDiameter = self.rotorDiameter
        rotorArea = self.rotorArea
        rho = self.air_density
        generator_efficiency = self.generator_efficiency
        Cp = self.Cp
        MU = self.parameters.MU
        axialIndProvided = self.parameters.axialIndProvided

        p_near0 = self.p_near0

        ############### BACKWARD MODE #############
        velocitiesTurbinesb = np.eye(nbdirs, nTurbines)
        wt_powerb = np.zeros((nbdirs, nTurbines))
        powerb = np.zeros(nbdirs)

        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct,axialInduction,axialIndProvided, keCorrCT, Region2CT, ke, \
                                      Vinf, keCorrArray, turbineXw, p_near0, rotorDiameter, MU, rho, Cp, \
                                      generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        dvelocitiesTurbines = np.hstack((wakeOverlapTRel_vb[:, :], Ctb[:, :], axialInductionb[:, :], turbineXwb[:, :], \
                                         rotorDiameterb[:, :], Cpb[:, :]))

        velocitiesTurbinesb[:, :] = 0.0
        wt_powerb = np.eye(nbdirs, nTurbines)

        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct,axialInduction,axialIndProvided, keCorrCT, Region2CT, ke, \
                                      Vinf, keCorrArray, turbineXw, p_near0, rotorDiameter, MU, rho, Cp, \
                                      generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        dwt_power = np.hstack((wakeOverlapTRel_vb[:, :], Ctb[:, :], axialInductionb[:, :], turbineXwb[:, :], \
                               rotorDiameterb[:, :], Cpb[:, :]))

        wt_powerb[:, :] = 0.0
        powerb[0] = 1.0

        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct,axialInduction,axialIndProvided, keCorrCT, Region2CT, ke, \
                                      Vinf, keCorrArray, turbineXw, p_near0, rotorDiameter, MU, rho, Cp, \
                                      generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        dpower = np.hstack((wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, rotorDiameterb, Cpb))
        dpower = dpower[0, :]

        # print wakeOverlapTRel_vb.shape, Ctb.shape, axialInductionb.shape, turbineXwb.shape, rotorDiameterb.shape, Cpb.shape

        J = np.vstack((dvelocitiesTurbines, dwt_power, dpower))

        # print J.shape

        return J