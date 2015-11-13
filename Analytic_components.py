from openmdao.main.api import Component, VariableTree
from openmdao.lib.datatypes.api import Array, Bool, Float, VarTree
import numpy as np
from Parameters import FLORISParameters
import time


class floris_adjustCtCp(Component):
    """ Adjust Cp and Ct to yaw if they are not already adjusted """

    parameters = VarTree(FLORISParameters(), iotype='in')

    def __init__(self, nTurbines):

        print 'entering adjustCtCp __init__ - analytic'

        super(floris_adjustCtCp, self).__init__()

        # Explicitly size input arrays
        self.add('Ct_in', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                                desc='Thrust coefficient for all turbines'))
        self.add('Cp_in', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                                desc='power coefficient for all turbines'))
        self.add('generator_efficiency', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                desc='generator efficiency of all turbines'))
        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', desc='yaw of each turbine'))

        # Explicitly size output arrays
        self.add('Ct_out', Array(np.zeros(nTurbines), iotype='out', dtype='float', \
                                 desc='Thrust coefficient for all turbines'))
        self.add('Cp_out', Array(np.zeros(nTurbines), iotype='out', dtype='float', \
                                 desc='power coefficient for all turbines'))

    def execute(self):

        print 'entering adjustCtCP - analytic'

        # print 'CTcorrected is', self.parameters.CTcorrected
        # print 'CPcorrected is', self.parameters.CPcorrected

        Ct = self.Ct_in
        Cp = self.Cp_in
        nTurbines = np.size(Ct)
        yaw = self.yaw*np.pi/180.
        # CTcorrected = self.parameters.CTcorrected
        # CPcorrected = self.parameters.CPcorrected
        # pP = self.parameters.pP

        # print 'before', Ct, Cp
        # print 'yaw in adjust = ', yaw
        # print 'Ct in adjust = ', Ct

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

        if not CTcorrected:
            # print Ct.size, yaw.size
            self.Ct_out = Ct*np.cos(yaw)*np.cos(yaw)
            dCt_dCt = np.eye(nTurbines)*np.cos(yaw)*np.cos(yaw)
            dCt_dyaw = np.eye(nTurbines)*(-2.*Ct*np.sin(yaw)*np.cos(yaw))*np.pi/180.
            dCt_dCp = np.zeros((nTurbines, nTurbines))
            dCt = np.hstack((dCt_dCt, dCt_dCp, dCt_dyaw))

        else:

            self.Ct_out = Ct
            dCt_dCt = np.eye(nTurbines, nTurbines)
            dCt_dCp = np.zeros((nTurbines, nTurbines))
            dCt_dyaw = np.zeros((nTurbines, nTurbines))
            dCt = np.hstack((dCt_dCt, dCt_dCp, dCt_dyaw))

        if not CPcorrected:

            self.Cp_out = Cp*np.cos(yaw)**pP

            dCp_dCp = np.eye(nTurbines, nTurbines)*np.cos(yaw)**pP
            dCp_dyaw = np.eye(nTurbines, nTurbines)*(-Cp*pP*np.sin(yaw)*np.cos(yaw)**(pP-1.0))*np.pi/180.
            dCp_dCt = np.zeros((nTurbines, nTurbines))
            dCp = np.hstack((dCp_dCt, dCp_dCp, dCp_dyaw))

        else:

            self.Cp_out = Cp
            dCp_dCp = np.eye(nTurbines, nTurbines)
            dCp_dCt = np.zeros((nTurbines, nTurbines))
            dCp_dyaw = np.zeros((nTurbines, nTurbines))
            dCp = np.hstack((dCp_dCt, dCp_dCp, dCp_dyaw))

        self.J = np.vstack((dCt, dCp))

    def list_deriv_vars(self):
        return ('Ct_in', 'Cp_in', 'yaw'), ('Ct_out', 'Cp_out')

    def provideJ(self):
        return self.J


class floris_windframe(Component):
    """ Calculates the locations of each turbine in the wind direction reference frame """

    # original variables
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    # flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    wind_direction = Float(iotype='in', units='deg', desc='wind direction using direction to, in deg. ccw from east')

    def __init__(self, nTurbines, resolution):

        print 'entering windframe __init__ - analytic'

        super(floris_windframe, self).__init__()

        # Explicitly size input arrays
        self.add('turbineX', Array(np.zeros(nTurbines), iotype='in', \
                                   desc='x positions of turbines in original ref. frame'))
        self.add('turbineY', Array(np.zeros(nTurbines), iotype='in', \
                                   desc='y positions of turbines in original ref. frame'))

        # variables for verbosity
        self.add('Ct', Array(np.zeros(nTurbines), iotype='in'))
        self.add('Cp', Array(np.zeros(nTurbines), iotype='in', \
                             desc='power coefficient for all turbines'))
        self.add('axialInduction', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                                         desc='axial induction of all turbines'))
        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', \
                              desc='yaw of each turbine'))

        # variables for testing wind speed at various locations
        self.add('ws_position', Array(np.zeros(resolution*resolution), iotype='in', units='m', \
                                      desc='position of desired measurements in original ref. frame'))

        # Explicitly size output arrays
        self.add('wsw_position', Array(np.zeros(resolution*resolution), iotype='out', units='m', \
                                       desc='position of desired measurements in wind ref. frame'))

        # for testing purposes only
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='out', units='m', \
                                    desc='x coordinates of turbines in wind dir. ref. frame'))
        self.add('turbineYw', Array(np.zeros(nTurbines), iotype='out', units='m', \
                                    desc='y coordinates of turbines in wind dir. ref. frame'))

    def execute(self):

        print 'entering windframe - analytic'

        if self.parameters.FLORISoriginal:

            ke = 0.065
            keCorrCT = 0.0
            # ignore these for now
            # keCorrTI = 0.0
            # keCorrHR = 0.0
            # keCorrHRTI = 0.0
            keCorrArray = 0.0

            kd = 0.15
            # kdCorrYawDirection = 0.0  # ignored for now

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

        Vinf = self.wind_speed
        windDirection = self.wind_direction*np.pi/180.0

        #variables to satisfy verbosity
        axialInd = self.axialInduction
        Cp = self.Cp
        Ct = self.Ct
        yaw = self.yaw*np.pi/180

        if self.verbose:
            np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print "wind direction %s deg" % [windDirection*180.0/np.pi]
            print "free-stream wind speed %s" % Vinf
            print "axial induction turbines %s" % axialInd
            print "C_P turbines %s" % Cp
            print "C_T turbines %s" % Ct
            print "yaw turbines %s" % yaw

        # get turbine positions and velocity sampling positions
        # position = self.position
        # turbineX = position[:, 0]
        # turbineY = position[:, 1]
        turbineX = self.turbineX
        turbineY = self.turbineY
        # print turbineX, turbineY

        if self.ws_position.any():
            velX = self.ws_position[:, 0]
            velY = self.ws_position[:, 1]
        else:
            velX = np.zeros([0, 0])
            velY = np.zeros([0, 0])

        # convert to downwind-crosswind coordinates
        rotationMatrix = np.array([(np.cos(-windDirection), -np.sin(-windDirection)),
                                   (np.sin(-windDirection), np.cos(-windDirection))])
        # print 'rotation matrix = ', rotationMatrix
        turbineLocations = np.dot(rotationMatrix, np.array([turbineX, turbineY]))
        # print turbineLocations
        self.turbineXw = np.zeros(turbineX.size)
        self.turbineYw = np.zeros(turbineX.size)
        # print self.turbineXw
        self.turbineXw = turbineLocations[0]
        self.turbineYw = turbineLocations[1]
        #print 'windframe.turbineX = %s' %self.turbineX
        if velX.size > 0:
            locations = np.dot(rotationMatrix, np.array([velX, velY]))
            velX = locations[0]
            velY = locations[1]

        self.wsw_position = np.array([velX, velY])
        #print 'wsw_position in windframe is:', self.wsw_position
        #print 'ws_position in windframe is:', self.ws_position

        # print self.turbineXw

    def list_deriv_vars(self):
        """specifies the inputs and outputs where derivatives are defined"""

        return('turbineX', 'turbineY'), ('turbineXw', 'turbineYw')

    def provideJ(self):

        #print 'entering windframe - provideJ'

        n = np.size(self.turbineX)

        windDirection = self.wind_direction*np.pi/180

        dturbineXw_dturbineX = np.zeros([n, n])
        dturbineXw_dturbineY = np.zeros([n, n])
        dturbineYw_dturbineX = np.zeros([n, n])
        dturbineYw_dturbineY = np.zeros([n, n])

        for i in range(0, n):
            dturbineXw_dturbineX[i, i] = np.cos(-windDirection)
            dturbineXw_dturbineY[i, i] = -np.sin(-windDirection)
            dturbineYw_dturbineX[i, i] = np.sin(-windDirection)
            dturbineYw_dturbineY[i, i] = np.cos(-windDirection)

        JturbineXw = np.concatenate((dturbineXw_dturbineX, dturbineXw_dturbineY), 1)
        JturbineYw = np.concatenate((dturbineYw_dturbineX, dturbineYw_dturbineY), 1)
        J = np.concatenate((JturbineXw, JturbineYw), 0)

        return J


class AEP(Component):

    AEP = Float(iotype='out', units='kW', desc='total annual energy output of wind farm')

    def __init__(self, nDirections):

        super(AEP, self).__init__()

        self.add('power_directions', Array(np.zeros(nDirections), iotype='in', units='kW', desc='vector containing \
                                           the power production at each wind direction ccw from north'))
        self.add('windrose_frequencies', Array(np.zeros(nDirections), iotype='in', desc='vector containing \
                                               the weighted frequency of wind at each direction ccw from east using \
                                               direction too'))
        # do not use these for any gradient calculations, only for output
        self.add('power_directions_out', Array(np.zeros(nDirections), iotype='out', units='kW', desc='vector containing \
                                           the power production at each wind direction ccw from north', deriv_ignore=True))

    def execute(self):

        #print 'in AEP'

        # locally name input values
        power_directions = self.power_directions
        windrose_frequencies = self.windrose_frequencies

        # number of hours in a year
        hours = 8760.0

        # calculate approximate AEP
        AEP = sum(power_directions*windrose_frequencies)*hours

        # promote AEP result to class attribute
        self.AEP = AEP
        self.power_directions_out = power_directions

        #print 'AEP %s' % self.AEP

    def list_deriv_vars(self):

        # return ('power_directions',), ('AEP', 'power_directions_out')
        return ('power_directions',), ('AEP',)

    def provideJ(self):

        #print 'entering AEP - provideJ'

        # create local variables
        windrose_frequencies = self.windrose_frequencies
        ndirs = np.size(windrose_frequencies)

        # number of hours in a year
        hours = 8760.0

        # calculate the derivative of outputs w.r.t. each wind direction
        dAEP_dpower = np.ones(ndirs)*windrose_frequencies*hours

        J = np.array([dAEP_dpower])

        return J


class dist_const(Component):

    parameters = VarTree(FLORISParameters(), iotype='in')

    def __init__(self, nTurbines):

        #print 'entering dist_const __init__'

        super(dist_const, self).__init__()

        # Explicitly size input arrays
        self.add('turbineX', Array(np.zeros(nTurbines), iotype='in', \
                                    desc='x coordinates of turbines in wind dir. ref. frame'))
        self.add('turbineY', Array(np.zeros(nTurbines), iotype='in', \
                                    desc='y coordinates of turbines in wind dir. ref. frame'))

        # Explicitly size output array
        self.add('separation', Array(np.zeros((nTurbines-1.)*nTurbines/2.), iotype='out', dtype='float', \
                                        desc='spacing of all turbines in the wind farm'))

    def execute(self):

        #print 'in dist const'

        turbineX = self.turbineX
        turbineY = self.turbineY
        nTurbines = turbineX.size
        separation = np.zeros((nTurbines-1.)*nTurbines/2.)

        k = 0
        for i in range(0, nTurbines):
            for j in range(i+1, nTurbines):
                separation[k] = np.sqrt((turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2)
                k += 1
        self.separation = separation

    def list_deriv_vars(self):
        return ('turbineX', 'turbineY'), ('separation',)

    def provideJ(self):

        #print 'entering dist const - provideJ'
	tictot = time.time()
        turbineX = self.turbineX
        turbineY = self.turbineY
        nTurbines = turbineX.size
        J = np.zeros(((nTurbines-1.)*nTurbines/2., 2*nTurbines))
	
        k = 0
        for i in range(0, nTurbines):
            for j in range(i+1, nTurbines):
                J[k, j] = (turbineX[j]-turbineX[i])*((turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2)**(-0.5)
                J[k, i] = (turbineX[i]-turbineX[j])*((turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2)**(-0.5)
                J[k, j+nTurbines] = (turbineY[j]-turbineY[i])*((turbineX[j]-turbineX[i])**2 +
                                                               (turbineY[j]-turbineY[i])**2)**(-0.5)
                J[k, i+nTurbines] = (turbineY[i]-turbineY[j])*((turbineX[j]-turbineX[i])**2 +
                                                               (turbineY[j]-turbineY[i])**2)**(-0.5)
                k += 1
	toctot = time.time()
	#print 'done %s' % (toctot-tictot)
        return J

class hull_const(Component):

    def __init__(self, nVertices, nTurbines):

        super(hull_const, self).__init__()

        # Explicitly size input arrays
        self.add('AX', Array(np.zeros(nVertices), iotype='in'))
        self.add('AY', Array(np.zeros(nVertices), iotype='in'))
        self.add('b', Array(np.zeros(nVertices), iotype='in'))

        self.add('turbineX', Array(np.zeros(nTurbines), iotype='in', \
                                    desc='x coordinates of turbines in wind dir. ref. frame'))
        self.add('turbineY', Array(np.zeros(nTurbines), iotype='in', \
                                    desc='y coordinates of turbines in wind dir. ref. frame'))

        # Explicitly size output array 
        # (vector with positive elements if turbines outside of hull)
        self.add('inout', Array(np.zeros(nVertices*nTurbines), iotype='out'))


    def execute(self):

        #print 'in hull const'
        tictot = time.time()

        AX = self.AX
        AY = self.AY
        b = self.b
        turbineX = self.turbineX
        turbineY = self.turbineY
        nTurbines = turbineX.size

        J = hull_const_J(AX, AY, nTurbines)

        self.inout = (np.dot(J, np.concatenate((turbineX,turbineY))) - np.tile(b, (1,nTurbines))).flatten()

        toctot = time.time()
        #print 'done %s' % (toctot-tictot)

    def list_deriv_vars(self):
        return ('turbineX', 'turbineY',), ('inout',)

    def provideJ(self):

        #print 'in hull const - provide J'
        tictot = time.time()

        AX = self.AX
        AY = self.AY
        b = self.b
        nTurbines = self.turbineX.size

        J = hull_const_J(AX, AY, nTurbines)

        toctot = time.time()
        #print 'done %s' % (toctot-tictot)

        return J

def hull_const_J(AX, AY, nTurbines):
    J = np.concatenate((np.kron(np.eye(nTurbines),AX).transpose(),np.kron(np.eye(nTurbines),AY).transpose()),1)
    return J

