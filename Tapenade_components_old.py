from openmdao.main.api import Component, VariableTree
from openmdao.lib.datatypes.api import Array, Bool, Float, VarTree
from Parameters import FLORISParameters
import numpy as np
import _floris


class floris_wcent_wdiam(Component):
    """ Calculates the center and diameter of each turbine wake at each other turbine """

    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    def __init__(self, nTurbines):

        # print 'entering wcent_wdiam __init__ - Tapenade'


        super(floris_wcent_wdiam, self).__init__()

        # Explicitly size input arrays
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='in', \
                                    desc='x coordinates of turbines in wind dir. ref. frame'))
        self.add('turbineYw', Array(np.zeros(nTurbines), iotype='in', \
                                    desc='y coordinates of turbines in wind dir. ref. frame'))
        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', desc='yaw of each turbine'))
        self.add('rotorDiameter', Array(np.zeros(nTurbines), dtype='float', iotype='in', \
                                        desc='rotor diameter of each turbine'))
        self.add('Ct', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                             desc='thrust coefficient of each turbine'))

        # Explicitly size output arrays
        self.add('wakeCentersYT', Array(np.zeros(nTurbines*nTurbines), iotype='out', dtype='float', \
                                        desc='wake center y position at each turbine'))
        self.add('wakeDiametersT', Array(np.zeros(3*nTurbines*nTurbines), iotype='out', dtype='float', \
                                         desc='wake diameter of each zone of each wake at each turbine'))

    def execute(self):

        # print 'entering wcent_wdiam - tapenade'

        # rename inputs and outputs
        # pP = self.parameters.pP
        # kd = self.parameters.kd
        # ke = self.parameters.ke
        # initialWakeDisplacement = self.parameters.initialWakeDisplacement
        # initialWakeAngle = self.parameters.initialWakeAngle
        rotorDiameter = self.rotorDiameter
        Ct = self.Ct
        # keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        # me = self.parameters.me

        if self.parameters.FLORISoriginal:

            ke = 0.065
            keCorrCT = 0.0
            keCorrTI = 0.0
            keCorrHR = 0.0
            keCorrHRTI = 0.0
            keCorrArray = 0.0

            kd = 0.15
            kdCorrDirection = 0.0

            me = np.array([-0.5, 0.22, 1.0])
            MU = np.array([0.5, 1.0, 5.5])
            pP = 1.88
            useWakeAngle = False
            initialWakeDisplacement = 4.5
            initialWakeAngle = self.parameters.initialWakeAngle
            bd = -0.01
            useaUbU = True
            aU = 5.0
            bU = 1.66
            adjustInitialWakeDiamToYaw = False

        else:
            # rename inputs and outputs
            ke = self.parameters.ke

            keCorrCT = self.parameters.keCorrCT
            keCorrTI = self.parameters.keCorrTI
            keCorrHR = self.parameters.keCorrHR
            keCorrHRTI = self.parameters.keCorrHRTI
            keCorrArray = self.parameters.keCorrArray

            kd = self.parameters.kd
            kdCorrYawDirection = self.parameters.kdCorrYawDirection

            me = self.parameters.me
            MU = self.parameters.MU

            initialWakeDisplacement = self.parameters.initialWakeDisplacement
            useWakeAngle = self.parameters.useWakeAngle
            initialWakeAngle = self.parameters.initialWakeAngle

            bd = self.parameters.bd

            useaUbU = self.parameters.useaUbU
            aU = self.parameters.aU
            bU = self.parameters.bU

            adjustInitialWakeDiamToYaw = self.parameters.adjustInitialWakeDiamToYaw

            pP = self.parameters.pP

        baselineCT = self.parameters.baselineCT
        baselineTI = self.parameters.baselineTI
        keSaturation = self.parameters.keSaturation

        CTcorrected = self.parameters.CTcorrected
        CPcorrected = self.parameters.CPcorrected
        axialIndProvided = self.parameters.axialIndProvided

        # x and y positions w.r.t. the wind direction (wind = +x)
        turbineXw = self.turbineXw
        turbineYw = self.turbineYw

        # yaw in degrees
        yaw_deg = self.yaw
        # print yaw_deg, Ct
        wakeCentersYT_vec, wakeDiametersT_vec = _floris.floris_wcent_wdiam(kd, initialWakeDisplacement, \
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me, bd, useWakeAngle, adjustInitialWakeDiamToYaw)

        # Outputs in vector form so they can be used in Jacobian creation
        self.wakeCentersYT = wakeCentersYT_vec
        self.wakeDiametersT = wakeDiametersT_vec

    def list_deriv_vars(self):

        return ('yaw', 'Ct', 'turbineXw', 'turbineYw', 'rotorDiameter'), ('wakeCentersYT', 'wakeDiametersT')

    def provideJ(self):

        # rename inputs and outputs
        # pP = self.parameters.pP
        # kd = self.parameters.kd
        # ke = self.parameters.ke
        # initialWakeDisplacement = self.parameters.initialWakeDisplacement
        # initialWakeAngle = self.parameters.initialWakeAngle
        rotorDiameter = self.rotorDiameter
        Ct = self.Ct
        # keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        # me = self.parameters.me

        if self.parameters.FLORISoriginal:

            ke = 0.065
            keCorrCT = 0.0
            keCorrTI = 0.0
            keCorrHR = 0.0
            keCorrHRTI = 0.0
            keCorrArray = 0.0

            kd = 0.15
            kdCorrDirection = 0.0

            me = np.array([-0.5, 0.22, 1.0])
            MU = np.array([0.5, 1.0, 5.5])
            pP = 1.88
            useWakeAngle = False
            initialWakeDisplacement = 4.5
            initialWakeAngle = self.parameters.initialWakeAngle
            bd = -0.01
            useaUbU = True
            aU = 5.0
            bU = 1.66
            adjustInitialWakeDiamToYaw = False

        else:
            # rename inputs and outputs
            ke = self.parameters.ke

            keCorrCT = self.parameters.keCorrCT
            keCorrTI = self.parameters.keCorrTI
            keCorrHR = self.parameters.keCorrHR
            keCorrHRTI = self.parameters.keCorrHRTI
            keCorrArray = self.parameters.keCorrArray

            kd = self.parameters.kd
            kdCorrYawDirection = self.parameters.kdCorrYawDirection

            me = self.parameters.me
            MU = self.parameters.MU

            initialWakeDisplacement = self.parameters.initialWakeDisplacement
            useWakeAngle = self.parameters.useWakeAngle
            initialWakeAngle = self.parameters.initialWakeAngle

            bd = self.parameters.bd

            useaUbU = self.parameters.useaUbU
            aU = self.parameters.aU
            bU = self.parameters.bU

            adjustInitialWakeDiamToYaw = self.parameters.adjustInitialWakeDiamToYaw

            pP = self.parameters.pP

        baselineCT = self.parameters.baselineCT
        baselineTI = self.parameters.baselineTI
        keSaturation = self.parameters.keSaturation

        CTcorrected = self.parameters.CTcorrected
        CPcorrected = self.parameters.CPcorrected
        axialIndProvided = self.parameters.axialIndProvided

        # x and y positions w.r.t. the wind direction (wind = +x)
        turbineXw = self.turbineXw
        turbineYw = self.turbineYw

        # turbine yaw w.r.t. wind direction
        yaw_deg = self.yaw

        # number of turbines
        nTurbines = np.size(turbineXw)

        # number of directions being differentiated in the Jacobian
        nbdirs = 3*nTurbines*nTurbines

        # input arrays to direct differentiation
        wakeDiametersT_vecb = np.zeros((nbdirs, 3*nTurbines*nTurbines))
        wakeCentersYT_vecb = np.eye(nbdirs, nTurbines*nTurbines)

        # function call to extract gradients of wakeCentersYT w.r.t. all design vars

        yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb, _, _ = _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement,
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me, bd, useWakeAngle, adjustInitialWakeDiamToYaw, wakeCentersYT_vecb, \
                              wakeDiametersT_vecb)

        # construct Jacobian of wakeCentersYT
        dwc = np.hstack((yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb))
        dwc = dwc[:nTurbines*nTurbines, :]

        # input arrays to direct differentiation
        wakeCentersYT_vecb[:, :] = 0.0
        wakeDiametersT_vecb = np.eye(nbdirs, nbdirs)

        # function call to extract gradients of wakeDiametersT w.r.t. all design vars
        yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb, _, _ = _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement,
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me, bd, useWakeAngle, adjustInitialWakeDiamToYaw, wakeCentersYT_vecb, \
                              wakeDiametersT_vecb)

        # construct Jacobian of wakeDiametersT
        dwd = np.hstack((yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb))

        # construct total Jacobian of floris_wcent_wdiam
        J = np.vstack((dwc, dwd))

        return J


class floris_overlap(Component):
    """ Calculates the overlap between each turbine rotor and the existing turbine wakes """

    def __init__(self, nTurbines):

        # print 'entering overlap __init__ - Tapenade'

        super(floris_overlap, self).__init__()

        # Explicitly size input arrays
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='in', units='m', ignore_deriv=True, \
                                    desc='X positions of turbines wrt the wind direction'))
        self.add('turbineYw', Array(np.zeros(nTurbines), iotype='in', units='m', \
                                          desc='Y positions of turbines wrt the wind direction'))
        self.add('rotorDiameter', Array(np.zeros(nTurbines), iotype='in', units='m', \
                                              desc='diameters of all turbine rotors'))
        self.add('wakeCentersYT', Array(np.zeros(nTurbines*nTurbines), iotype='in', units='m', \
                                              desc='Y positions of all wakes at each turbine'))
        self.add('wakeDiametersT', Array(3*nTurbines*np.zeros(nTurbines), iotype='in', units='m',\
                                               desc='diameters of all turbines wake zones'))

        # Explicitly size output arrays
        self.add('wakeOverlapTRel', Array(np.zeros(3*nTurbines*nTurbines), iotype='out', \
                                                desc='relative wake zone overlap to rotor area'))


    def execute(self):

        # print 'entering overlap - Tapenade'

        # call to fortran code to obtain relative wake overlap values
        wakeOverlapTRel_vec = _floris.floris_overlap(self.turbineXw, self.turbineYw, self.rotorDiameter, \
                                                     self.wakeDiametersT, self.wakeCentersYT)

        # pass results to self in the form of a vector for use in Jacobian creation
        self.wakeOverlapTRel = wakeOverlapTRel_vec

    def list_deriv_vars(self):
        """specifies the inputs and outputs where derivatives are defined"""

        return ('turbineYw', 'rotorDiameter', 'wakeDiametersT', 'wakeCentersYT'), ('wakeOverlapTRel',)

    def provideJ(self):

        # rename input variables
        turbineXw = self.turbineXw
        turbineYw = self.turbineYw
        rotorDiameter = self.rotorDiameter
        wakeDiametersT_vec = self.wakeDiametersT
        wakeCentersYT_vec = self.wakeCentersYT
        wakeOverlapTRel_vec = self.wakeOverlapTRel

        # number of turbines
        nTurbines = np.size(turbineXw)

        # number of directions being differentiated
        nbdirs = 3*nTurbines*nTurbines

        # input array to direct differentiation
        wakeOverlapTRel_vecb = np.eye(nbdirs, 3*nTurbines*nTurbines)

        # function call to fortran to obtain gradients
        turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb \
            = _floris.floris_overlap_bv(turbineXw, turbineYw, rotorDiameter, wakeDiametersT_vec, \
                                         wakeCentersYT_vec, wakeOverlapTRel_vec,  wakeOverlapTRel_vecb)

        # construct Jacobian of floris_overlap
        J = np.hstack((turbineYwb[:, :], rotorDiameterb[:, :], wakeDiametersT_vecb[:, :], \
                             wakeCentersYT_vecb[:, :]))

        return J


class floris_power(Component):
    """ Calculates the turbine power and effective wind speed for each turbine """

    # inputs
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    # Flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    air_density = Float(iotype='in', units='kg/(m*m*m)', desc='air density in free stream')

    # output
    power = Float(iotype='out', units='kW', desc='total power output of the wind farm')

    def __init__(self, nTurbines):

        # print 'entering power __init__ - Tapenade'

        super(floris_power, self).__init__()



        # Explicitly size input arrays
        # input variables added so I don't have to use WISDEM while developing gradients
        self.add('rotorDiameter', Array(np.zeros(nTurbines), dtype='float', iotype='in', units='m', \
                                        desc='rotor diameters of all turbine'))
        self.add('axialInduction', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                                         desc='axial induction of all turbines'))
        self.add('Ct', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                             desc='Thrust coefficient for all turbines'))
        self.add('Cp', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                             desc='power coefficient for all turbines'))
        self.add('generator_efficiency', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                                               desc='generator efficiency of all turbines'))
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='in', dtype='float', units='m', \
                                    desc='X positions of turbines in the wind direction reference frame'))
        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', dtype='float', units='deg', \
                                    desc='yaw angle of turbines wrt the wind direction'))
        self.add('wakeCentersYT',  Array(np.zeros(nTurbines*nTurbines), iotype='in', units='m', \
                                         desc='centers of the wakes at each turbine'))
        self.add('wakeDiametersT', Array(np.zeros(3*nTurbines*nTurbines), iotype='in', units='m', \
                                         desc='diameters of each of the wake zones for each of the wakes \
                                         at each turbine'))
        self.add('wakeOverlapTRel', Array(np.zeros(3*nTurbines*nTurbines), iotype='in', units='m', \
                                          desc='ratios of wake overlap area per zone to rotor area'))

        # Explicitly size output arrays
        self.add('velocitiesTurbines', Array(np.zeros(nTurbines), iotype='out', units='m/s'))
        self.add('wt_power', Array(np.zeros(nTurbines), iotype='out', units='kW'))

    def execute(self):
        # print 'entering power - tapenade'

        # reassign input variables
        wakeOverlapTRel_v = self.wakeOverlapTRel
        # ke = self.parameters.ke
        # keCorrArray = self.parameters.keCorrArray
        # keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        Ct = self.Ct
        Vinf = self.wind_speed
        turbineXw = self.turbineXw
        axialInduction = self.axialInduction
        rotorDiameter = self.rotorDiameter
        rho = self.air_density
        generator_efficiency = self.generator_efficiency
        Cp = self.Cp
        yaw = self.yaw
        # MU = self.parameters.MU
        # axialIndProvided = self.parameters.axialIndProvided

        if self.parameters.FLORISoriginal:

            ke = 0.065
            keCorrCT = 0.0
            # ignore these for now
            # keCorrTI = 0.0
            # keCorrHR = 0.0
            # keCorrHRTI = 0.0
            keCorrArray = 0.0

            kd = 0.15
            # kdCorrDirection = 0.0

            me = np.array([-0.5, 0.22, 1.0])
            MU = np.array([0.5, 1.0, 5.5])
            pP = 1.88
            useWakeAngle = False
            initialWakeDisplacement = 4.5
            bd = -0.01
            useaUbU = True
            aU = 5.0
            bU = 1.66
            adjustInitialWakeDiamToYaw = False

        else:
            # rename inputs and outputs
            ke = self.parameters.ke

            keCorrCT = self.parameters.keCorrCT
            keCorrTI = self.parameters.keCorrTI
            keCorrHR = self.parameters.keCorrHR
            keCorrHRTI = self.parameters.keCorrHRTI
            keCorrArray = self.parameters.keCorrArray

            kd = self.parameters.kd
            kdCorrYawDirection = self.parameters.kdCorrYawDirection

            me = self.parameters.me
            MU = self.parameters.MU

            initialWakeDisplacement = self.parameters.initialWakeDisplacement
            useWakeAngle = self.parameters.useWakeAngle
            initialWakeAngle = self.parameters.initialWakeAngle

            bd = self.parameters.bd

            useaUbU = self.parameters.useaUbU
            aU = self.parameters.aU
            bU = self.parameters.bU

            adjustInitialWakeDiamToYaw = self.parameters.adjustInitialWakeDiamToYaw

            pP = self.parameters.pP
        # print 'in floris_power aU = %s, bU = %s, useaUbU = %s, ke = %s, kd = %s' % (aU, bU, useaUbU, ke, kd)
        baselineCT = self.parameters.baselineCT
        baselineTI = self.parameters.baselineTI
        keSaturation = self.parameters.keSaturation

        CTcorrected = self.parameters.CTcorrected
        CPcorrected = self.parameters.CPcorrected
        axialIndProvided = self.parameters.axialIndProvided

        # how far in front of turbines to use overlap power calculations (in rotor diameters). This must match the
        # value used in floris_wcent_wdiam (hardcoded in fortran as 1)
        # TODO hard code this parameter in the fortran code and remove the specifier from all functions of this component
        p_near0 = 1.0

        # pass p_near0 to self for use in gradient calculations
        self.p_near0 = p_near0

        # call to fortran code to obtain output values
        velocitiesTurbines, wt_power, power = _floris.floris_power(wakeOverlapTRel_v, Ct, axialInduction, \
                                                            axialIndProvided, useaUbU, keCorrCT, Region2CT, ke, \
                                                            Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter, MU, \
                                                            rho, aU, bU, Cp, generator_efficiency)

        # optional print statements
        if self.verbose:
            print "wind speed at turbines %s [m/s]" % velocitiesTurbines
            print "rotor area %d" % (np.pi*rotorDiameter[0]*rotorDiameter[0]/4.0)
            print "rho %s" % rho
            print "generator_efficiency %s" % generator_efficiency
            print "powers turbines %s [kW]" % wt_power

        # pass outputs to self
        self.velocitiesTurbines = velocitiesTurbines
        self.wt_power = wt_power
        self.power = power




    def list_deriv_vars(self):
        """specifies the inputs and outputs where derivatives are defined"""

        return ('wakeOverlapTRel', 'Ct', 'axialInduction', 'turbineXw', 'yaw', \
                'rotorDiameter', 'Cp'), ('velocitiesTurbines', 'wt_power', 'power')
        # return ('turbineX', 'turbineY'), ('turbineXw', 'turbineYw', 'wsw_position')

    def provideJ(self):

        # number of turbines
        nTurbines = np.size(self.turbineXw)

        # number of directions to differentiate
        nbdirs = nTurbines

        # reassign input variables
        wakeOverlapTRel_v = self.wakeOverlapTRel
        # ke = self.parameters.ke
        # keCorrArray = self.parameters.keCorrArray
        # keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        Ct = self.Ct
        Vinf = self.wind_speed
        turbineXw = self.turbineXw
        axialInduction = self.axialInduction
        rotorDiameter = self.rotorDiameter
        rho = self.air_density
        generator_efficiency = self.generator_efficiency
        Cp = self.Cp
        yaw = self.yaw
        # MU = self.parameters.MU
        # axialIndProvided = self.parameters.axialIndProvided

        if self.parameters.FLORISoriginal:

            ke = 0.065
            keCorrCT = 0.0
            # ignore these for now
            # keCorrTI = 0.0
            # keCorrHR = 0.0
            # keCorrHRTI = 0.0
            keCorrArray = 0.0

            kd = 0.15
            # kdCorrDirection = 0.0

            me = np.array([-0.5, 0.22, 1.0])
            MU = np.array([0.5, 1.0, 5.5])
            pP = 1.88
            useWakeAngle = False
            initialWakeDisplacement = 4.5
            bd = -0.01
            useaUbU = True
            aU = 5.0
            bU = 1.66
            adjustInitialWakeDiamToYaw = False

        else:
            # rename inputs and outputs
            ke = self.parameters.ke

            keCorrCT = self.parameters.keCorrCT
            keCorrTI = self.parameters.keCorrTI
            keCorrHR = self.parameters.keCorrHR
            keCorrHRTI = self.parameters.keCorrHRTI
            keCorrArray = self.parameters.keCorrArray

            kd = self.parameters.kd
            kdCorrYawDirection = self.parameters.kdCorrYawDirection

            me = self.parameters.me
            MU = self.parameters.MU

            initialWakeDisplacement = self.parameters.initialWakeDisplacement
            useWakeAngle = self.parameters.useWakeAngle
            initialWakeAngle = self.parameters.initialWakeAngle

            bd = self.parameters.bd

            useaUbU = self.parameters.useaUbU
            aU = self.parameters.aU
            bU = self.parameters.bU

            adjustInitialWakeDiamToYaw = self.parameters.adjustInitialWakeDiamToYaw

            pP = self.parameters.pP

        baselineCT = self.parameters.baselineCT
        baselineTI = self.parameters.baselineTI
        keSaturation = self.parameters.keSaturation

        CTcorrected = self.parameters.CTcorrected
        CPcorrected = self.parameters.CPcorrected
        axialIndProvided = self.parameters.axialIndProvided

        # see execute(self) for explanation
        p_near0 = self.p_near0

        # input arrays to direct differentiation
        velocitiesTurbinesb = np.eye(nbdirs, nTurbines)
        wt_powerb = np.zeros((nbdirs, nTurbines))
        powerb = np.zeros(nbdirs)

        # call to fortran to obtain gradients of velocitiesTurbines
        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct, axialInduction, axialIndProvided, useaUbU, keCorrCT, \
                                      Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter, MU, \
                                      rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        # construct Jacobian of velocitiesTurbines
        dvelocitiesTurbines = np.hstack((wakeOverlapTRel_vb[:, :], Ctb[:, :], axialInductionb[:, :], turbineXwb[:, :], \
                                         yawb[:, :], rotorDiameterb[:, :], Cpb[:, :]))

        # input arrays to direct differentiation
        velocitiesTurbinesb[:, :] = 0.0
        wt_powerb = np.eye(nbdirs, nTurbines)

        # call to fortran to obtain gradients wt_power
        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct, axialInduction, axialIndProvided, useaUbU, keCorrCT, \
                                      Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter, MU, \
                                      rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        # construct Jacobian of wt_power
        dwt_power = np.hstack((wakeOverlapTRel_vb[:, :], Ctb[:, :], axialInductionb[:, :], turbineXwb[:, :], \
                               yawb[:, :], rotorDiameterb[:, :], Cpb[:, :]))

        # input arrays to direct differentiation
        wt_powerb[:, :] = 0.0
        powerb[0] = 1.0

        # call to fortran to obtain gradients of power
        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct, axialInduction, axialIndProvided, useaUbU, keCorrCT, \
                                      Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter, MU, \
                                      rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        # construct Jacobian of power
        dpower = np.hstack((wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb))
        dpower = dpower[0., :]

        # construct total Jacobian of floris_power
        J = np.vstack((dvelocitiesTurbines, dwt_power, dpower))
        # print 'power J = ', J
        # print 'in J'
        return J