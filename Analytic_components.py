from openmdao.main.api import Component, VariableTree
from openmdao.lib.datatypes.api import Array, Bool, Float, VarTree
import numpy as np
from Parameters import FLORISParameters


class floris_adjustCtCp(Component):
    """ Adjust Cp and Ct to yaw if they are not already adjusted """

    parameters = VarTree(FLORISParameters(), iotype='in')
    Ct_in = Array(iotype='in', dtype='float', desc='Thrust coefficient for all turbines')
    Cp_in = Array(iotype='in', dtype='float', desc='power coefficient for all turbines')
    Ct_out = Array(iotype='out', dtype='float', desc='Thrust coefficient for all turbines')
    Cp_out = Array(iotype='out', dtype='float', desc='power coefficient for all turbines')
    generator_efficiency = Array(iotype='in', dtype='float', desc='generator efficiency of all turbines')
    yaw = Array(iotype='in', desc='yaw of each turbine')

    def execute(self):

        # print 'entering adjustCtCP - analytic'

        # print 'CTcorrected is', self.parameters.CTcorrected
        # print 'CPcorrected is', self.parameters.CPcorrected

        Ct = self.Ct_in
        Cp = self.Cp_in
        nTurbines = np.size(Ct)
        yaw = self.yaw*np.pi/180.
        CTcorrected = self.parameters.CTcorrected
        CPcorrected = self.parameters.CPcorrected
        pP = self.parameters.pP

        # print 'before', Ct, Cp
        # print 'yaw in adjust = ', yaw
        # print 'Ct in adjust = ', Ct

        if not CTcorrected:

            self.Ct_out = Ct*np.cos(yaw)*np.cos(yaw)
            #
            # dCt_dCt = np.eye(nTurbines, nTurbines)*np.cos(yaw)*np.cos(yaw)
            # dCt_dyaw = np.eye(nTurbines, nTurbines)*(-2*np.sin(yaw)*np.cos(yaw))
            # dCt_dCp = np.zeros((nTurbines, nTurbines))
            # dCt = np.hstack((dCt_dCt, dCt_dCp, dCt_dyaw))

        else:

            self.Ct_out = Ct
            # dCt_dCt = np.eye(nTurbines, nTurbines)*np.cos(yaw)*np.cos(yaw)
            # dCt_dCp = np.zeros((nTurbines, nTurbines))
            # dCt_dyaw = np.zeros((nTurbines, nTurbines))
            # dCt = np.hstack((dCt_dCt, dCt_dCp, dCt_dyaw))

        if not CPcorrected:

            self.Cp_out = Cp*np.cos(yaw)**pP

            # dCp_dCp = np.eye(nTurbines, nTurbines)*np.cos(yaw)**pP
            # dCp_dyaw = np.eye(nTurbines, nTurbines)*(-Cp*pP*np.sin(yaw)*np.cos(yaw)**(pP-1.0))
            # dCp_dCt = np.zeros((nTurbines, nTurbines))
            # dCp = np.hstack((dCp_dCt, dCp_dCp, dCp_dyaw))

        else:

            self.Cp_out = Cp
            # dCp_dCp = np.eye(nTurbines, nTurbines)*np.cos(yaw)**pP
            # dCp_dCt = np.zeros((nTurbines, nTurbines))
            # dCp_dyaw = np.zeros((nTurbines, nTurbines))
            # dCp = np.hstack((dCp_dCt, dCp_dCp, dCp_dyaw))
        #
    #     self.J = np.vstack((dCt, dCp))
    #
    # def list_deriv_vars(self):
    #
    #     return ('Ct_in', 'Cp_in', 'yaw'), ('Ct_out', 'Cp_out')
    #
    #
    # def provideJ(self):
    #
    #     return self.J


class floris_windframe(Component):
    """ Calculates the locations of each turbine in the wind direction reference frame """

    # original variables
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')
    turbineX = Array(iotype='in', desc='x positions of turbines in original ref. frame')
    turbineY = Array(iotype='in', desc='y positions of turbines in original ref. frame')

    # variables for verbosity
    Ct = Array(iotype='in', deriv_ignore='true', dtype='float')
    Cp = Array(iotype='in', deriv_ignore='true', dtype='float', desc='power coefficient for all turbines')
    axialInduction = Array(iotype='in', dtype='float', desc='axial induction of all turbines')
    yaw = Array(iotype='in', deriv_ignore=True, desc='yaw of each turbine')

    # variables for testing wind speed at various locations
    ws_position = Array(iotype='in', units='m', deriv_ignore=True, desc='position of desired measurements in original ref. frame')
    wsw_position = Array(iotype='out', units='m', deriv_ignore=True, desc='position of desired measurements in wind ref. frame')

    # flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    wind_direction = Float(iotype='in', units='deg', desc='overall wind direction for wind farm')

    # for testing purposes only
    turbineXw = Array(iotype='out', units='m', desc='x coordinates of turbines in wind dir. ref. frame')
    turbineYw = Array(iotype='out', units='m', desc='y coordinates of turbines in wind dir. ref. frame')

    def execute(self):

        # print 'entering windframe - analytic'

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


class floris_windrose(Component):
    """ Adjust Cp and Ct to yaw if they are not already adjusted """

    parameters = VarTree(FLORISParameters(), iotype='in')
    Ct_in = Array(iotype='in', dtype='float', desc='Thrust coefficient for all turbines')
    Cp_in = Array(iotype='in', dtype='float', desc='power coefficient for all turbines')
    Ct_out = Array(iotype='out', dtype='float', desc='Thrust coefficient for all turbines')
    Cp_out = Array(iotype='out', dtype='float', desc='power coefficient for all turbines')
    generator_efficiency = Array(iotype='in', dtype='float', desc='generator efficiency of all turbines')
    yaw = Array(iotype='in', desc='yaw of each turbine')

    def execute(self):

        # print 'entering adjustCtCP - analytic'

        # print 'CTcorrected is', self.parameters.CTcorrected
        # print 'CPcorrected is', self.parameters.CPcorrected

        Ct = self.Ct_in
        Cp = self.Cp_in
        nTurbines = np.size(Ct)
        yaw = self.yaw*np.pi/180.
        CTcorrected = self.parameters.CTcorrected
        CPcorrected = self.parameters.CPcorrected
        pP = self.parameters.pP

        # print 'before', Ct, Cp
        # print 'yaw in adjust = ', yaw
        # print 'Ct in adjust = ', Ct

        if not CTcorrected:

            self.Ct_out = Ct*np.cos(yaw)*np.cos(yaw)
            #
            # dCt_dCt = np.eye(nTurbines, nTurbines)*np.cos(yaw)*np.cos(yaw)
            # dCt_dyaw = np.eye(nTurbines, nTurbines)*(-2*np.sin(yaw)*np.cos(yaw))
            # dCt_dCp = np.zeros((nTurbines, nTurbines))
            # dCt = np.hstack((dCt_dCt, dCt_dCp, dCt_dyaw))

        else:

            self.Ct_out = Ct
            # dCt_dCt = np.eye(nTurbines, nTurbines)*np.cos(yaw)*np.cos(yaw)
            # dCt_dCp = np.zeros((nTurbines, nTurbines))
            # dCt_dyaw = np.zeros((nTurbines, nTurbines))
            # dCt = np.hstack((dCt_dCt, dCt_dCp, dCt_dyaw))

        if not CPcorrected:

            self.Cp_out = Cp*np.cos(yaw)**pP

            # dCp_dCp = np.eye(nTurbines, nTurbines)*np.cos(yaw)**pP
            # dCp_dyaw = np.eye(nTurbines, nTurbines)*(-Cp*pP*np.sin(yaw)*np.cos(yaw)**(pP-1.0))
            # dCp_dCt = np.zeros((nTurbines, nTurbines))
            # dCp = np.hstack((dCp_dCt, dCp_dCp, dCp_dyaw))

        else:

            self.Cp_out = Cp
            # dCp_dCp = np.eye(nTurbines, nTurbines)*np.cos(yaw)**pP
            # dCp_dCt = np.zeros((nTurbines, nTurbines))
            # dCp_dyaw = np.zeros((nTurbines, nTurbines))
            # dCp = np.hstack((dCp_dCt, dCp_dCp, dCp_dyaw))
        #
    #     self.J = np.vstack((dCt, dCp))
    #
    # def list_deriv_vars(self):
    #
    #     return ('Ct_in', 'Cp_in', 'yaw'), ('Ct_out', 'Cp_out')
    #
    #
    # def provideJ(self):
    #
    #     return self.J