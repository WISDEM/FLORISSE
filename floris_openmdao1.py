import numpy as np

from openmdao.api import Group, Component, Problem, IndepVarComp

import _floris


class FLORIS(Group):
    """ Group containing all necessary components of the floris model """

    def __init__(self, nTurbines, resolution):
        super(FLORIS, self).__init__()
        self.add('f_1', floris_windframe(nTurbines, resolution), promotes=['*'])
        self.add('f_2', floris_wcent_wdiam(nTurbines), promotes=['*'])
        self.add('f_3', floris_overlap(nTurbines), promotes=['*'])
        self.add('f_4', floris_power(nTurbines), promotes=['*'])


class adjustCtCp_yaw(Component):
    """ Adjust Cp and Ct to yaw if they are not already adjusted """

    def __init__(self, nTurbines):

        print 'entering adjustCtCp __init__ - analytic'
        super(adjustCtCp_yaw, self).__init__()

        # Explicitly size input arrays
        self.add_param('Ct_in', val=np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_param('Cp_in', val=np.zeros(nTurbines), desc='power coefficient for all turbines')
        self.add_param('yaw', val=np.zeros(nTurbines), desc='yaw of each turbine')

        # Explicitly size output arrays
        self.add_output('Ct_out', val=np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_output('Cp_out', val=np.zeros(nTurbines), desc='power coefficient for all turbines')

        # parameters since var trees are not supports
        self.add_param('floris_params:pP', 1.88)
        self.add_param('floris_params:CTcorrected', False,
                       desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('floris_params:CPcorrected', False,
                       desc='CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)')
        self.add_param('floris_params:FLORISoriginal', False,
                       desc='override all parameters and use FLORIS as original in first Wind Energy paper')

    def solve_nonlinear(self, params, unknowns, resids):

        print 'entering adjustCtCP - analytic'

        # collect inputs
        Ct = params['Ct_in']
        Cp = params['Cp_in']
        yaw = params['yaw'] * np.pi / 180.

        # determine floris_parameter values
        if params['floris_params:FLORISoriginal']:
            pP = 1.88
        else:
            pP = params['floris_params:pP']

        CTcorrected = params['floris_params:CTcorrected']
        CPcorrected = params['floris_params:CPcorrected']

        # calculate new CT values, if desired
        if not CTcorrected:
            unknowns['Ct_out'] = Ct * np.cos(yaw) * np.cos(yaw)
        else:
            unknowns['Ct_out'] = Ct

        # calculate new CP values, if desired
        if not CPcorrected:
            unknowns['Cp_out'] = Cp * np.cos(yaw) ** pP
        else:
            unknowns['Cp_out'] = Cp

    def linearize(self, params, unknowns, resids):

        # collect inputs
        Ct = params['Ct_in']
        Cp = params['Cp_in']
        nTurbines = np.size(Ct)
        yaw = params['yaw'] * np.pi / 180.

        # determine floris_parameter values
        if params['floris_params:FLORISoriginal']:
            pP = 1.88
        else:
            pP = params['floris_params:pP']

        CTcorrected = params['floris_params:CTcorrected']
        CPcorrected = params['floris_params:CPcorrected']

        # calculate gradients
        J = {}

        if not CTcorrected:
            J[('Ct_out', 'Ct_in')] = np.eye(nTurbines) * np.cos(yaw) * np.cos(yaw)
            J[('Ct_out', 'Cp_in')] = np.zeros((nTurbines, nTurbines))
            J[('Ct_out', 'yaw')] = np.eye(nTurbines) * (-2. * Ct * np.sin(yaw) * np.cos(yaw)) * np.pi / 180.
        else:
            J[('Ct_out', 'Ct_in')] = np.eye(nTurbines, nTurbines)
            J[('Ct_out', 'Cp_in')] = np.zeros((nTurbines, nTurbines))
            J[('Ct_out', 'yaw')] = np.zeros((nTurbines, nTurbines))

        if not CPcorrected:
            J[('Cp_out', 'Cp_in')] = np.eye(nTurbines, nTurbines) * np.cos(yaw) ** pP
            J[('Cp_out', 'Ct_in')] = np.zeros((nTurbines, nTurbines))
            J[('Cp_out', 'yaw')] = np.eye(nTurbines, nTurbines) * (
                -Cp * pP * np.sin(yaw) * np.cos(yaw) ** (pP - 1.0)) * np.pi / 180.
        else:
            J[('Cp_out', 'Cp_in')] = np.eye(nTurbines, nTurbines)
            J[('Cp_out', 'Ct_in')] = np.zeros((nTurbines, nTurbines))
            J[('Cp_out', 'yaw')] = np.zeros((nTurbines, nTurbines))

        return J


class floris_windframe(Component):
    """ Calculates the locations of each turbine in the wind direction reference frame """

    def __init__(self, nTurbines, resolution):

        print 'entering windframe __init__ - analytic'

        super(floris_windframe, self).__init__()

        self.nTurbines = nTurbines

        # flow property variables
        self.add_param('wind_speed', val=8.0, units='m/s', desc='free stream wind velocity')
        self.add_param('wind_direction', val=270, units='deg',
                       desc='wind direction using direction from, in deg. cw from north as in meteorological data')

        # Explicitly size input arrays
        self.add_param('turbineX', val=np.zeros(nTurbines), desc='x positions of turbines in original ref. frame')
        self.add_param('turbineY', val=np.zeros(nTurbines), desc='y positions of turbines in original ref. frame')

        # variables for testing wind speed at various locations
        self.add_param('ws_position', val=np.zeros(resolution * resolution), units='m',
                       desc='position of desired measurements in original ref. frame')

        # Explicitly size output arrays
        self.add_param('wsw_position', val=np.zeros(resolution * resolution), units='m',
                       desc='position of desired measurements in wind ref. frame')

        # add output
        self.add_output('turbineXw', val=np.zeros(nTurbines), units='m', desc='downwind coordinates of turbines')
        self.add_output('turbineYw', val=np.zeros(nTurbines), units='m', desc='crosswind coordinates of turbines')

        # define floris parameters
        self.add_param('floris_params:pP', 1.88)
        self.add_param('floris_params:ke', 0.065)
        self.add_param('floris_params:keCorrArray', 0.0)
        self.add_param('floris_params:keCorrCT', 0.0)
        self.add_param('floris_params:Region2CT', 4.0*(1.0/3.0)*(1.0-(1.0/3.0)))
        self.add_param('floris_params:kd', 0.15)
        self.add_param('floris_params:me', np.array([-0.5, 0.22, 1.0]))
        self.add_param('floris_params:initialWakeDisplacement', -4.5)
        self.add_param('floris_params:initialWakeAngle', 3.0)
        self.add_param('floris_params:baselineCT', 4./3.*(1.-1./3.))
        self.add_param('floris_params:keCorrTI', 0.0)
        self.add_param('floris_params:baselineTI', 0.045)
        self.add_param('floris_params:keCorrHR', 0.0) # neutral, with heating rate 0, is baseline
        self.add_param('floris_params:keCorrHRTI', 0.0)
        self.add_param('floris_params:keSaturation', 0.0)
        self.add_param('floris_params:kdCorrYawDirection', 0.0)
        self.add_param('floris_params:MU', np.array([0.5, 1.0, 10]))
        self.add_param('floris_params:CTcorrected', True, desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('floris_params:CPcorrected', True, desc = 'CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)')
        self.add_param('floris_params:axialIndProvided', False, desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('floris_params:useWakeAngle', True)
        self.add_param('floris_params:bd', -0.01)
        self.add_param('floris_params:useaUbU', False)
        self.add_param('floris_params:aU', 5.0, units='deg')
        self.add_param('floris_params:bU', 1.66)
        self.add_param('floris_params:adjustInitialWakeDiamToYaw', True)
        self.add_param('floris_params:FLORISoriginal', False, desc='override all parameters and use FLORIS as original in first Wind Energy paper')

    def solve_nonlinear(self, params, unknowns, resids):

        windDirectionDeg = params['wind_direction']

        # get turbine positions and velocity sampling positions
        turbineX = params['turbineX']
        turbineY = params['turbineY']

        # if self.ws_position.any():
        #     velX = self.ws_position[:, 0]
        #     velY = self.ws_position[:, 1]
        # else:
        #     velX = np.zeros([0, 0])
        #     velY = np.zeros([0, 0])

        # convert to downwind(x)-crosswind(y) coordinates
        windDirectionDeg = 270. - windDirectionDeg
        if windDirectionDeg < 0.:
            windDirectionDeg += 360.
        windDirectionRad = np.pi*windDirectionDeg/180.0             # inflow wind direction in radians
        # rotationMatrix = np.array([(np.cos(-windDirectionRad), -np.sin(-windDirectionRad)),
        #                            (np.sin(-windDirectionRad), np.cos(-windDirectionRad))])
        turbineXw = turbineX*np.cos(-windDirectionRad)-turbineY*np.sin(-windDirectionRad)
        turbineYw = turbineX*np.sin(-windDirectionRad)+turbineY*np.cos(-windDirectionRad)
        # turbineLocations = np.dot(rotationMatrix, np.array([turbineX, turbineY]))
        # turbineLocations = rotationMatrix*np.array([turbineX, turbineY])
        # print turbineLocations
        # unknowns['turbineXw'] = turbineLocations[0]
        # unknowns['turbineXw'] = turbineLocations[0]
        unknowns['turbineXw'] = turbineXw
        unknowns['turbineYw'] = turbineYw


        #print 'windframe.turbineX = %s' %self.turbineX
        # if velX.size > 0:
        #     locations = np.dot(rotationMatrix, np.array([velX, velY]))
        #     velX = locations[0]
        #     velY = locations[1]

        # self.wsw_position = np.array([velX, velY])
        #print 'wsw_position in windframe is:', self.wsw_position
        #print 'ws_position in windframe is:', self.ws_position

        # print self.turbineXw

    def linearize(self, params, unknowns, resids):

        print 'entering windframe - provideJ'

        nTurbines = self.nTurbines

        windDirectionDeg = params['wind_direction']

        windDirectionDeg = 270. - windDirectionDeg
        if windDirectionDeg < 0.:
            windDirectionDeg += 360.

        windDirectionRad = np.pi*windDirectionDeg/180.0             # inflow wind direction in radians

        dturbineXw_dturbineX = np.ones([nTurbines, nTurbines])*np.cos(-windDirectionRad)
        dturbineXw_dturbineY = np.ones([nTurbines, nTurbines])*(-np.sin(-windDirectionRad))
        dturbineYw_dturbineX = np.ones([nTurbines, nTurbines])*np.sin(-windDirectionRad)
        dturbineYw_dturbineY = np.ones([nTurbines, nTurbines])*np.cos(-windDirectionRad)

        print dturbineXw_dturbineY.shape
        J = {}

        J[('turbineXw', 'turbineX')] = dturbineXw_dturbineX
        J[('turbineXw', 'turbineY')] = dturbineXw_dturbineY
        J[('turbineYw', 'turbineX')] = dturbineYw_dturbineX
        J[('turbineYw', 'turbineY')] = dturbineYw_dturbineY

        return J


class floris_wcent_wdiam(Component):
    """ Calculates the center and diameter of each turbine wake at each other turbine """

    def __init__(self, nTurbines):

        print 'entering wcent_wdiam __init__ - Tapenade'

        super(floris_wcent_wdiam, self).__init__()

        # input arrays
        self.add_param('turbineXw', np.zeros(nTurbines), desc='x coordinates of turbines in wind dir. ref. frame')
        self.add_param('turbineYw', np.zeros(nTurbines), desc='y coordinates of turbines in wind dir. ref. frame')
        self.add_param('yaw', np.zeros(nTurbines), desc='yaw of each turbine')
        self.add_param('rotorDiameter', np.zeros(nTurbines), desc='rotor diameter of each turbine')
        self.add_param('Ct', np.zeros(nTurbines), desc='thrust coefficient of each turbine')

        # output arrays
        self.add_output('wakeCentersYT', np.zeros(nTurbines*nTurbines), desc='wake center y position at each turbine')
        self.add_output('wakeDiametersT', np.zeros(3*nTurbines*nTurbines), desc='wake diameter of each zone of each wake at each turbine')

        # FLORIS parameters
        self.add_param('floris_params:pP', 1.88)
        self.add_param('floris_params:ke', 0.065)
        self.add_param('floris_params:keCorrArray', 0.0)
        self.add_param('floris_params:keCorrCT', 0.0)
        self.add_param('floris_params:Region2CT', 4.0*(1.0/3.0)*(1.0-(1.0/3.0)))
        self.add_param('floris_params:kd', 0.15)
        self.add_param('floris_params:me', np.array([-0.5, 0.22, 1.0]))
        self.add_param('floris_params:initialWakeDisplacement', -4.5)
        self.add_param('floris_params:initialWakeAngle', 3.0)
        self.add_param('floris_params:baselineCT', 4./3.*(1.-1./3.))
        self.add_param('floris_params:keCorrTI', 0.0)
        self.add_param('floris_params:baselineTI', 0.045)
        self.add_param('floris_params:keCorrHR', 0.0) # neutral, with heating rate 0, is baseline
        self.add_param('floris_params:keCorrHRTI', 0.0)
        self.add_param('floris_params:keSaturation', 0.0)
        self.add_param('floris_params:kdCorrYawDirection', 0.0)
        self.add_param('floris_params:MU', np.array([0.5, 1.0, 10]))
        self.add_param('floris_params:CTcorrected', True, desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('floris_params:CPcorrected', True, desc = 'CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)')
        self.add_param('floris_params:axialIndProvided', False, desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('floris_params:useWakeAngle', True)
        self.add_param('floris_params:bd', -0.01)
        self.add_param('floris_params:useaUbU', False)
        self.add_param('floris_params:aU', 5.0, units='deg')
        self.add_param('floris_params:bU', 1.66)
        self.add_param('floris_params:adjustInitialWakeDiamToYaw', True)
        self.add_param('floris_params:FLORISoriginal', False, desc='override all parameters and use FLORIS as original in first Wind Energy paper')

    def solve_nonlinear(self, params, unknowns, resids):

        print 'entering wcent_wdiam - tapenade'

        rotorDiameter = params['rotorDiameter']
        Ct = params['Ct']

        Region2CT = params['floris_params:Region2CT']

        if params['floris_params:FLORISoriginal']:
            ke = 0.065
            keCorrCT = 0.0
            kd = 0.15
            me = np.array([-0.5, 0.22, 1.0])
            useWakeAngle = False
            initialWakeDisplacement = 4.5
            initialWakeAngle = params['floris_params:initialWakeAngle']
            bd = -0.01
            adjustInitialWakeDiamToYaw = False

        else:
            ke = params['floris_params:ke']
            keCorrCT = params['floris_params:keCorrCT']
            kd = params['floris_params:kd']
            me = params['floris_params:me']
            initialWakeDisplacement = params['floris_params:initialWakeDisplacement']
            useWakeAngle = params['floris_params:useWakeAngle']
            initialWakeAngle = params['floris_params:initialWakeAngle']
            bd = params['floris_params:bd']
            adjustInitialWakeDiamToYaw = params['floris_params:adjustInitialWakeDiamToYaw']


        # x and y positions w.r.t. the wind direction (wind = +x)
        turbineXw = params['turbineXw']
        turbineYw = params['turbineYw']

        # yaw in degrees
        yaw_deg = params['yaw']
        # print yaw_deg, Ct
        wakeCentersYT_vec, wakeDiametersT_vec = _floris.floris_wcent_wdiam(kd, initialWakeDisplacement, \
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me, bd, useWakeAngle, adjustInitialWakeDiamToYaw)

        # Outputs in vector form so they can be used in Jacobian creation
        unknowns['wakeCentersYT'] = wakeCentersYT_vec
        unknowns['wakeDiametersT'] = wakeDiametersT_vec

    def linearize(self, params, unknowns, resids):

        rotorDiameter = params['rotorDiameter']
        Ct = params['Ct']
        Region2CT = params['floris_params:Region2CT']

        if params['floris_params:FLORISoriginal']:

            ke = 0.065
            keCorrCT = 0.0
            kd = 0.15
            me = np.array([-0.5, 0.22, 1.0])
            useWakeAngle = False
            initialWakeDisplacement = 4.5
            initialWakeAngle = params['floris_params:initialWakeAngle']
            bd = -0.01
            adjustInitialWakeDiamToYaw = False

        else:
            # rename inputs and outputs
            ke = params['floris_params:ke']
            keCorrCT = params['floris_params:keCorrCT']
            kd = params['floris_params:kd']
            me = params['floris_params:me']
            initialWakeDisplacement = params['floris_params:initialWakeDisplacement']
            useWakeAngle = params['floris_params:useWakeAngle']
            initialWakeAngle = params['floris_params:initialWakeAngle']
            bd = params['floris_params:bd']
            adjustInitialWakeDiamToYaw = params['floris_params:adjustInitialWakeDiamToYaw']

        # x and y positions w.r.t. the wind direction (wind = +x)
        turbineXw = params['turbineXw']
        turbineYw = params['turbineYw']

        # turbine yaw w.r.t. wind direction
        yaw_deg = params['yaw']

        # number of turbines
        nTurbines = np.size(turbineXw)

        # number of directions being differentiated in the Jacobian
        nbdirs = 3*nTurbines*nTurbines

        # input arrays to direct differentiation
        wakeDiametersT_vecb = np.zeros((nbdirs, 3*nTurbines*nTurbines))
        wakeCentersYT_vecb = np.eye(nbdirs, nTurbines*nTurbines)

        # initialize linearize dict
        J = {}

        # function call to extract gradients of wakeCentersYT w.r.t. all design vars
        yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb, _, _ = _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement,
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me, bd, useWakeAngle, adjustInitialWakeDiamToYaw, wakeCentersYT_vecb, \
                              wakeDiametersT_vecb)

        # construct Jacobian of wakeCentersYT
        J['wakeCentersYT', 'yaw'] = yawb
        J['wakeCentersYT', 'Ct'] = Ctb
        J['wakeCentersYT', 'turbineXw'] = turbineXwb
        J['wakeCentersYT', 'turbineYw'] = turbineYwb
        J['wakeCentersYT', 'rotorDiameter'] = rotorDiameterb

        # input arrays to direct differentiation
        wakeCentersYT_vecb[:, :] = 0.0
        wakeDiametersT_vecb = np.eye(nbdirs, nbdirs)

        # function call to extract gradients of wakeDiametersT w.r.t. all design vars
        yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb, _, _ = _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement,
							  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                              rotorDiameter, me, bd, useWakeAngle, adjustInitialWakeDiamToYaw, wakeCentersYT_vecb, \
                              wakeDiametersT_vecb)

        # construct Jacobian of wakeDiametersT
        J['wakeDiametersT', 'yaw'] = yawb
        J['wakeDiametersT', 'Ct'] = Ctb
        J['wakeDiametersT', 'turbineXw'] = turbineXwb
        J['wakeDiametersT', 'turbineYw'] = turbineYwb
        J['wakeDiametersT', 'rotorDiameter'] = rotorDiameterb

        return J


class floris_overlap(Component):
    """ Calculates the overlap between each turbine rotor and the existing turbine wakes """

    def __init__(self, nTurbines):

        print 'entering overlap __init__ - Tapenade'

        super(floris_overlap, self).__init__()

        # input arrays
        self.add_param('turbineXw', np.zeros(nTurbines), units='m',
                       desc='X positions of turbines wrt the wind direction')
        self.add_param('turbineYw', np.zeros(nTurbines), units='m',
                       desc='Y positions of turbines wrt the wind direction')
        self.add_param('rotorDiameter', np.zeros(nTurbines), units='m',
                       desc='diameters of all turbine rotors')
        self.add_param('wakeCentersYT', np.zeros(nTurbines*nTurbines), units='m',
                       desc='Y positions of all wakes at each turbine')
        self.add_param('wakeDiametersT', np.zeros(3*nTurbines*nTurbines), units='m',
                       desc='diameters of all turbines wake zones')

        # output arrays
        self.add_output('wakeOverlapTRel', np.zeros(3*nTurbines*nTurbines),
                        desc='relative wake zone overlap to rotor area')

        # etc
        self.nTurbines = nTurbines

    def solve_nonlinear(self, params, unknowns, resids):

        print 'entering overlap - Tapenade'

        # call to fortran code to obtain relative wake overlap values
        wakeOverlapTRel_vec = _floris.floris_overlap(params['turbineXw'], params['turbineYw'], params['rotorDiameter'], \
                                                     params['wakeDiametersT'], params['wakeCentersYT'])

        # pass results to self in the form of a vector for use in Jacobian creation
        unknowns['wakeOverlapTRel'] = wakeOverlapTRel_vec

    def linearize(self, params, unknowns, resids):

        # number of turbines
        nTurbines = self.nTurbines

        # number of directions being differentiated
        nbdirs = 3*nTurbines*nTurbines

        # input array to direct differentiation
        wakeOverlapTRel_vecb = np.eye(nbdirs, 3*nTurbines*nTurbines)

        # function call to fortran to obtain gradients
        turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb \
            = _floris.floris_overlap_bv(params['turbineXw'], params['turbineYw'], params['rotorDiameter'],
                                        params['wakeDiametersT'], params['wakeCentersYT'],
                                        unknowns['wakeOverlapTRel'], wakeOverlapTRel_vecb)

        # construct Jacobian of floris_overlap

        J = {}

        J['wakeOverlapTRel', 'turbineYw'] = turbineYwb
        J['wakeOverlapTRel', 'rotorDiameter'] = rotorDiameterb
        J['wakeOverlapTRel', 'wakeDiametersT'] = wakeDiametersT_vecb
        J['wakeOverlapTRel', 'wakeCentersYT'] = wakeCentersYT_vecb

        return J


class floris_power(Component):
    """ Calculates the turbine power and effective wind speed for each turbine """

    def __init__(self, nTurbines):

        print 'entering power __init__ - Tapenade'

        super(floris_power, self).__init__()

        self.nTurbines = nTurbines

        # inputs
        self.add_param('wind_speed', 8.0, units='m/s', desc='free stream wind velocity')
        self.add_param('air_density', 1.1716, units='kg/(m*m*m)', desc='air density in free stream')
        self.add_param('rotorDiameter', np.zeros(nTurbines), units='m', desc='rotor diameters of all turbine')
        self.add_param('axialInduction', np.zeros(nTurbines), desc='axial induction of all turbines')
        self.add_param('Ct', np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_param('Cp', np.zeros(nTurbines), desc='power coefficient for all turbines')
        self.add_param('generator_efficiency', np.zeros(nTurbines), desc='generator efficiency of all turbines')
        self.add_param('turbineXw', np.zeros(nTurbines), units='m',
                       desc='X positions of turbines in the wind direction reference frame')
        self.add_param('yaw', np.zeros(nTurbines), units='deg',
                       desc='yaw angle of turbines wrt the wind direction')
        self.add_param('wakeCentersYT',  np.zeros(nTurbines*nTurbines), units='m',
                       desc='centers of the wakes at each turbine')
        self.add_param('wakeDiametersT', np.zeros(3*nTurbines*nTurbines), units='m',
                       desc='diameters of each of the wake zones for each of the wakes at each turbine')
        self.add_param('wakeOverlapTRel', np.zeros(3*nTurbines*nTurbines), units='m',
                       desc='ratios of wake overlap area per zone to rotor area')

        # outputs
        self.add_output('velocitiesTurbines', np.zeros(nTurbines), units='m/s',
                       desc='effective hub velocity for each turbine')
        self.add_output('wt_power', np.zeros(nTurbines), units='kW', desc='power output of each turbine')
        # output
        self.add_output('power', 0.0, units='kW', desc='total power output of the wind farm')


        # connect floris_params
        self.add_param('floris_params:pP', 1.88)
        self.add_param('floris_params:ke', 0.065)
        self.add_param('floris_params:keCorrArray', 0.0)
        self.add_param('floris_params:keCorrCT', 0.0)
        self.add_param('floris_params:Region2CT', 4.0*(1.0/3.0)*(1.0-(1.0/3.0)))
        self.add_param('floris_params:kd', 0.15)
        self.add_param('floris_params:me', np.array([-0.5, 0.22, 1.0]))
        self.add_param('floris_params:initialWakeDisplacement', -4.5)
        self.add_param('floris_params:initialWakeAngle', 3.0)
        self.add_param('floris_params:baselineCT', 4./3.*(1.-1./3.))
        self.add_param('floris_params:keCorrTI', 0.0)
        self.add_param('floris_params:baselineTI', 0.045)
        self.add_param('floris_params:keCorrHR', 0.0) # neutral, with heating rate 0, is baseline
        self.add_param('floris_params:keCorrHRTI', 0.0)
        self.add_param('floris_params:keSaturation', 0.0)
        self.add_param('floris_params:kdCorrYawDirection', 0.0)
        self.add_param('floris_params:MU', np.array([0.5, 1.0, 10]))
        self.add_param('floris_params:CTcorrected', True,
                       desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('floris_params:CPcorrected', True,
                       desc='CP factor already corrected by CCBlade calculation '
                            '(assumed with approximately factor cos(yaw)^3)')
        self.add_param('floris_params:axialIndProvided', False,
                       desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('floris_params:useWakeAngle', True)
        self.add_param('floris_params:bd', -0.01)
        self.add_param('floris_params:useaUbU', False)
        self.add_param('floris_params:aU', 5.0, units='deg')
        self.add_param('floris_params:bU', 1.66)
        self.add_param('floris_params:adjustInitialWakeDiamToYaw', True)
        self.add_param('floris_params:FLORISoriginal', False,
                       desc='override all parameters and use FLORIS as original in first Wind Energy paper')

    def solve_nonlinear(self, params, unknowns, resids):
        print 'entering power - tapenade'

        # reassign input variables
        wakeOverlapTRel_v = params['wakeOverlapTRel']
        Region2CT = params['floris_params:Region2CT']
        Ct = params['Ct']
        Vinf = params['wind_speed']
        turbineXw = params['turbineXw']
        axialInduction = params['axialInduction']
        rotorDiameter = params['rotorDiameter']
        rho = params['air_density']
        generator_efficiency = params['generator_efficiency']
        Cp = params['Cp']
        yaw = params['yaw']

        # set floris parameter values
        if params['floris_params:FLORISoriginal']:

            ke = 0.065
            keCorrCT = 0.0
            keCorrArray = 0.0
            MU = np.array([0.5, 1.0, 5.5])
            useaUbU = True
            aU = 5.0
            bU = 1.66

        else:
            ke = params['floris_params:ke']
            keCorrCT = params['floris_params:keCorrCT']
            keCorrArray = params['floris_params:keCorrArray']
            MU = params['floris_params:MU']
            useaUbU = params['floris_params:useaUbU']
            aU = params['floris_params:aU']
            bU = params['floris_params:bU']


        axialIndProvided = params['floris_params:axialIndProvided']

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

        # pass outputs to self
        unknowns['velocitiesTurbines'] = velocitiesTurbines
        unknowns['wt_power'] = wt_power
        unknowns['power'] = power

    def linearize(self, params, unknowns, resids):

        # number of turbines
        nTurbines = self.nTurbines

        # number of directions to differentiate
        nbdirs = nTurbines

        # reassign input variables
        wakeOverlapTRel_v = params['wakeOverlapTRel']
        Region2CT = params['floris_params:Region2CT']
        Ct = params['Ct']
        Vinf = params['wind_speed']
        turbineXw = params['turbineXw']
        axialInduction = params['axialInduction']
        rotorDiameter = params['rotorDiameter']
        rho = params['air_density']
        generator_efficiency = params['generator_efficiency']
        Cp = params['Cp']
        yaw = params['yaw']

        # set floris parameter values
        if params['floris_params:FLORISoriginal']:

            ke = 0.065
            keCorrCT = 0.0
            keCorrArray = 0.0
            MU = np.array([0.5, 1.0, 5.5])
            useaUbU = True
            aU = 5.0
            bU = 1.66

        else:
            ke = params['floris_params:ke']
            keCorrCT = params['floris_params:keCorrCT']
            keCorrArray = params['floris_params:keCorrArray']
            MU = params['floris_params:MU']
            useaUbU = params['floris_params:useaUbU']
            aU = params['floris_params:aU']
            bU = params['floris_params:bU']

        axialIndProvided = params['floris_params:axialIndProvided']

        # see execute(self) for explanation
        p_near0 = self.p_near0

        # create jacobian dict
        J = {}

        # input arrays to direct differentiation
        velocitiesTurbinesb = np.eye(nbdirs, nTurbines)
        wt_powerb = np.zeros((nbdirs, nTurbines))
        powerb = np.zeros(nbdirs)

        # call to fortran to obtain gradients of velocitiesTurbines
        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct, axialInduction, axialIndProvided, useaUbU, keCorrCT, \
                                      Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter, MU, \
                                      rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        # collect values of the jacobian of velocitiesTurbines
        J['velocitiesTurbines', 'wakeOverlapTRel'] = wakeOverlapTRel_vb
        J['velocitiesTurbines', 'Ct'] = Ctb
        J['velocitiesTurbines', 'Cp'] = Cpb
        J['velocitiesTurbines', 'axialInduction'] = axialInductionb
        J['velocitiesTurbines', 'turbineXw'] = turbineXwb
        J['velocitiesTurbines', 'yaw'] = yawb
        J['velocitiesTurbines', 'rotorDiameter'] = rotorDiameterb

        # input arrays to direct differentiation
        velocitiesTurbinesb[:, :] = 0.0
        wt_powerb = np.eye(nbdirs, nTurbines)

        # call to fortran to obtain gradients wt_power
        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct, axialInduction, axialIndProvided, useaUbU, keCorrCT, \
                                      Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter, MU, \
                                      rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        # collect values of the jacobian of wt_power
        J['wt_power', 'wakeOverlapTRel'] = wakeOverlapTRel_vb
        J['wt_power', 'Ct'] = Ctb
        J['wt_power', 'Cp'] = Cpb
        J['wt_power', 'axialInduction'] = axialInductionb
        J['wt_power', 'turbineXw'] = turbineXwb
        J['wt_power', 'yaw'] = yawb
        J['wt_power', 'rotorDiameter'] = rotorDiameterb

        # input arrays to direct differentiation
        wt_powerb[:, :] = 0.0
        powerb[0] = 1.0

        # call to fortran to obtain gradients of power
        wakeOverlapTRel_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, Ct, axialInduction, axialIndProvided, useaUbU, keCorrCT, \
                                      Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter, MU, \
                                      rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb, wt_powerb, powerb)

        # collect values of the jacobian of wt_power
        J['power', 'wakeOverlapTRel'] = wakeOverlapTRel_vb
        J['power', 'Ct'] = Ctb
        J['power', 'Cp'] = Cpb
        J['power', 'axialInduction'] = axialInductionb
        J['power', 'turbineXw'] = turbineXwb
        J['power', 'yaw'] = yawb
        J['power', 'rotorDiameter'] = rotorDiameterb

        return J


if __name__ == "__main__":
    top = Problem()

    root = top.root = Group()

    root.add('p1', IndepVarComp('x', np.array([3.0])))
    root.add('p2', IndepVarComp('y', np.array([2.0])))
    root.add('p3', IndepVarComp('z', np.array([1.0])))
    root.add('p', adjustCtCp_yaw(nTurbines=np.array([1])))

    root.connect('p1.x', 'p.Ct_in')
    root.connect('p2.y', 'p.Cp_in')
    root.connect('p3.z', 'p.yaw')

    top.setup()
    top.run()

    print(root.p.unknowns['Ct_out'])
    print(root.p.unknowns['Cp_out'])

    top = Problem()

    root = top.root = Group()

    root.add('p1', IndepVarComp('x', np.array([10.0])))
    root.add('p2', IndepVarComp('y', np.array([10.0])))
    root.add('p3', IndepVarComp('z', 90))
    root.add('p', floris_windframe(nTurbines=np.array([1]), resolution=0))

    root.connect('p1.x', 'p.turbineX')
    root.connect('p2.y', 'p.turbineY')
    root.connect('p3.z', 'p.wind_direction')

    top.setup()
    top.run()

    print(root.p.unknowns['turbineXw'])
    print(root.p.unknowns['turbineYw'])

    top = Problem()

    root = top.root = Group()

    root.add('p1', IndepVarComp('x', np.array([10.0])))
    root.add('p2', IndepVarComp('y', np.array([10.0])))
    root.add('p', floris_wcent_wdiam(nTurbines=np.array([1])))

    root.connect('p1.x', 'p.turbineXw')
    root.connect('p2.y', 'p.turbineYw')
    root.connect('p1.x', 'p.yaw')
    root.connect('p1.x', 'p.rotorDiameter')
    root.connect('p1.x', 'p.Ct')

    top.setup()
    top.run()

    print(root.p.unknowns['wakeDiametersT'])
    print(root.p.unknowns['wakeCentersYT'])

    top = Problem()

    root = top.root = Group()

    root.add('p1', IndepVarComp('x', np.array([10.0])))
    root.add('p2', IndepVarComp('y', np.array([10.0, 10.0, 10.0])))
    root.add('p', floris_overlap(nTurbines=np.array([1])))

    root.connect('p1.x', 'p.turbineXw')
    root.connect('p1.x', 'p.turbineYw')
    root.connect('p1.x', 'p.rotorDiameter')
    root.connect('p1.x', 'p.wakeCentersYT')
    root.connect('p2.y', 'p.wakeDiametersT')

    top.setup()
    top.run()

    print(root.p.unknowns['wakeOverlapTRel'])

    top = Problem()

    root = top.root = Group()

    root.add('p1', IndepVarComp('x', np.array([10.0])))
    root.add('p2', IndepVarComp('y', np.array([10.0])))
    root.add('p', floris_power(nTurbines=np.array([1])))

    root.connect('p1.x', 'p.turbineXw')
    root.connect('p1.x', 'p.rotorDiameter')
    root.connect('p1.x', 'p.Ct')
    root.connect('p2.y', 'p.Cp')

    top.setup()
    top.run()

    print(root.p.unknowns['power'])
    print(root.p.unknowns['wt_power'])
    print(root.p.unknowns['velocitiesTurbines'])



# class floris_verbosity(Component):
#
#     def __init__(self, nTurbines, verbose):
#         super(floris_verbosity, self).__init__()
#         # variables for verbosity
#         self.add_param('Ct', np.zeros(nTurbines))
#         self.add_param('Cp', np.zeros(nTurbines), desc='power coefficient for all turbines')
#         self.add_param('axialInduction', np.zeros(nTurbines), desc='axial induction of all turbines')
#         self.add_param('yaw', np.zeros(nTurbines), desc='yaw of each turbine')
#
#         self.verbose = verbose
#
#     def solve_nonlinear(self, params, unknowns, resids):
#          # variables to satisfy verbosity
#         axialInd = params['axialInduction']
#         Cp = params['Cp']
#         Ct = params['Ct']
#         yaw = params['yaw * np.pi / 180']
#         windDirection = params['windDirection']
#         Vinf = params['Vinf']
#         verbose = self.verbose
#
#         if verbose:
#             np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
#             print "wind direction %s deg" % (windDirection * 180.0 / np.pi)
#             print "free-stream wind speed %s" % Vinf
#             print "axial induction turbines %s" % axialInd
#             print "C_P turbines %s" % Cp
#             print "C_T turbines %s" % Ct
#             print "yaw turbines %s" % yaw
#
#
#         # optional print statements from power
#         if verbose:
#             print "wind speed at turbines %s [m/s]" % velocitiesTurbines
#             print "rotor area %d" % (np.pi*rotorDiameter[0]*rotorDiameter[0]/4.0)
#             print "rho %s" % rho
#             print "generator_efficiency %s" % generator_efficiency
#             print "powers turbines %s [kW]" % wt_power
#
#     def linearize(self):
#         J = {}
#         return J
