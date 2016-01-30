from openmdao.api import Component, Group, Problem, IndepVarComp

import numpy as np
from scipy import interp


class WindFrame(Component):
    """ Calculates the locations of each turbine in the wind direction reference frame """

    def __init__(self, nTurbines, resolution):

        # print 'entering windframe __init__ - analytic'

        super(WindFrame, self).__init__()

        # set finite difference options (fd used for testing only)
        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-5
        self.fd_options['step_type'] = 'relative'

        self.nTurbines = nTurbines

        # flow property variables
        self.add_param('wind_speed', val=8.0, units='m/s', desc='free stream wind velocity')
        self.add_param('wind_direction', val=270.0, units='deg',
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

    def solve_nonlinear(self, params, unknowns, resids):

        windDirectionDeg = params['wind_direction']

        # get turbine positions and velocity sampling positions
        turbineX = params['turbineX']
        turbineY = params['turbineY']

        # print "in windframe", turbineX, turbineY

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

        unknowns['turbineXw'] = turbineXw
        unknowns['turbineYw'] = turbineYw

    def linearize(self, params, unknowns, resids):

        # print 'entering windframe - provideJ'

        nTurbines = self.nTurbines

        windDirectionDeg = params['wind_direction']

        windDirectionDeg = 270. - windDirectionDeg
        if windDirectionDeg < 0.:
            windDirectionDeg += 360.

        windDirectionRad = np.pi*windDirectionDeg/180.0             # inflow wind direction in radians

        dturbineXw_dturbineX = np.eye(nTurbines, nTurbines)*np.cos(-windDirectionRad)
        dturbineXw_dturbineY = np.eye(nTurbines, nTurbines)*(-np.sin(-windDirectionRad))
        dturbineYw_dturbineX = np.eye(nTurbines, nTurbines)*np.sin(-windDirectionRad)
        dturbineYw_dturbineY = np.eye(nTurbines, nTurbines)*np.cos(-windDirectionRad)

        # print dturbineXw_dturbineY.shape
        J = {}

        J[('turbineXw', 'turbineX')] = dturbineXw_dturbineX
        J[('turbineXw', 'turbineY')] = dturbineXw_dturbineY
        J[('turbineYw', 'turbineX')] = dturbineYw_dturbineX
        J[('turbineYw', 'turbineY')] = dturbineYw_dturbineY


        # print 'end windframe jacobian'

        return J


class AdjustCtCpYaw(Component):
    """ Adjust Cp and Ct to yaw if they are not already adjusted """

    def __init__(self, nTurbines):

        # print 'entering adjustCtCp __init__ - analytic'
        super(AdjustCtCpYaw, self).__init__()

        # set finite difference options (fd used for testing only)
        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-5
        self.fd_options['step_type'] = 'relative'

        # Explicitly size input arrays
        self.add_param('Ct_in', val=np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_param('Cp_in', val=np.zeros(nTurbines), desc='power coefficient for all turbines')
        self.add_param('yaw', val=np.zeros(nTurbines), units='deg', desc='yaw of each turbine')

        # Explicitly size output arrays
        self.add_output('Ct_out', val=np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_output('Cp_out', val=np.zeros(nTurbines), desc='power coefficient for all turbines')

        # parameters since var trees are not supports
        self.add_param('params:pP', 1.88, pass_by_obj=True)
        self.add_param('params:CTcorrected', False,
                       desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)', pass_by_obj=True)
        self.add_param('params:CPcorrected', False,
                       desc='CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)', pass_by_obj=True)
        self.add_param('floris_params:FLORISoriginal', True,
                       desc='override all parameters and use FLORIS as original in first Wind Energy paper', pass_by_obj=True)

    def solve_nonlinear(self, params, unknowns, resids):

        # print 'entering adjustCtCP - analytic'

        # collect inputs
        Ct = params['Ct_in']
        Cp = params['Cp_in']
        yaw = params['yaw'] * np.pi / 180.

        # determine floris_parameter values
        if params['floris_params:FLORISoriginal']:
            pP = 1.88
        else:
            pP = params['params:pP']

        CTcorrected = params['params:CTcorrected']
        CPcorrected = params['params:CPcorrected']

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

        # print 'entering CtCp linearize'
        # collect inputs
        Ct = params['Ct_in']
        Cp = params['Cp_in']
        nTurbines = np.size(Ct)
        yaw = params['yaw'] * np.pi / 180.

        # determine floris_parameter values
        if params['floris_params:FLORISoriginal']:
            pP = 1.88
        else:
            pP = params['params:pP']

        CTcorrected = params['params:CTcorrected']
        CPcorrected = params['params:CPcorrected']

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


class WindFarmAEP(Component):
    """ Estimate the AEP based on power production for each direction and weighted by wind direction frequency  """

    def __init__(self, nDirections):

        super(WindFarmAEP, self).__init__()

        # set finite difference options (fd used for testing only)
        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-5
        self.fd_options['step_type'] = 'relative'

        # define inputs
        self.add_param('power_directions', np.zeros(nDirections), units='kW',
                       desc='vector containing the power production at each wind direction ccw from north')
        self.add_param('windrose_frequencies', np.zeros(nDirections),
                       desc='vector containing the weighted frequency of wind at each direction ccw from east using '
                            'direction too')

        # define output
        self.add_output('AEP', val=0.0, units='kW', desc='total annual energy output of wind farm')

    def solve_nonlinear(self, params, unknowns, resids):

        # # print 'in AEP'

        # locally name input values
        power_directions = params['power_directions']
        windrose_frequencies = params['windrose_frequencies']

        # number of hours in a year
        hours = 8760.0

        # calculate approximate AEP
        AEP = sum(power_directions*windrose_frequencies)*hours

        # promote AEP result to class attribute
        unknowns['AEP'] = AEP

        print 'In AEP, AEP %s' % unknowns['AEP']

    def linearize(self, params, unknowns, resids):

        # # print 'entering AEP - provideJ'

        # assign params to local variables
        windrose_frequencies = params['windrose_frequencies']
        power_directions = params['power_directions']
        ndirs = np.size(windrose_frequencies)

        # number of hours in a year
        hours = 8760.0

        # calculate the derivative of outputs w.r.t. each wind direction
        dAEP_dpower = np.ones(ndirs)*windrose_frequencies*hours
        dAEP_dwindrose_frequencies = np.ones(ndirs)*power_directions*hours

        # print 'dAEP = ', dAEP_dpower
        J = {}

        J['AEP', 'power_directions'] = np.array([dAEP_dpower])
        J['AEP', 'windrose_frequencies'] = np.array([dAEP_dwindrose_frequencies])

        return J


class SpacingComp(Component):
    """
    Calculates inter-turbine spacing for all turbine pairs
    """

    def __init__(self, nTurbines):

        # print 'entering dist_const __init__

        super(SpacingComp, self).__init__()

        # set finite difference options (fd used for testing only)
        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-5
        self.fd_options['step_type'] = 'relative'

        # Explicitly size input arrays
        self.add_param('turbineX', val=np.zeros(nTurbines),
                       desc='x coordinates of turbines in wind dir. ref. frame')
        self.add_param('turbineY', val=np.zeros(nTurbines),
                       desc='y coordinates of turbines in wind dir. ref. frame')

        # Explicitly size output array
        self.add_output('separation', val=np.zeros((nTurbines-1.)*nTurbines/2.),
                        desc='spacing of all turbines in the wind farm')

    def solve_nonlinear(self, params, unknowns, resids):
        # print 'in dist const'

        turbineX = params['turbineX']
        turbineY = params['turbineY']
        nTurbines = turbineX.size
        separation = np.zeros((nTurbines-1.)*nTurbines/2.)

        k = 0
        for i in range(0, nTurbines):
            for j in range(i+1, nTurbines):
                separation[k] = np.sqrt((turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2)
                k += 1
        unknowns['separation'] = separation

    def linearize(self, params, unknowns, resids):
        # print 'entering dist const - linearize'
        turbineX = params['turbineX']
        turbineY = params['turbineY']
        # print turbineX
        # print turbineY
        nTurbines = turbineX.size
        dS = np.zeros(((nTurbines-1.)*nTurbines/2., 2*nTurbines))
        k = 0
        # print 'in dist_const, turbineX = ', turbineX
        # print 'in dist_const, turbineY = ', turbineY

        for i in range(0, nTurbines):
            for j in range(i+1, nTurbines):

                if turbineX[i] != turbineX[j] or turbineY[i] != turbineY[j]:
                    # separation wrt Xj
                    dS[k, j] = (turbineX[j]-turbineX[i])*((turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2)**(-0.5)
                    # separation wrt Xi
                    dS[k, i] = (turbineX[i]-turbineX[j])*((turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2)**(-0.5)
                    dS[k, j+nTurbines] = (turbineY[j]-turbineY[i])*((turbineX[j]-turbineX[i])**2 +
                                                                    (turbineY[j]-turbineY[i])**2)**(-0.5)
                    dS[k, i+nTurbines] = (turbineY[i]-turbineY[j])*((turbineX[j]-turbineX[i])**2 +
                                                                    (turbineY[j]-turbineY[i])**2)**(-0.5)
                else:
                    # separation wrt Xj
                    dS[k, j] = 1.
                    # separation wrt Xi
                    dS[k, i] = 1.
                    dS[k, j+nTurbines] = 1.
                    dS[k, i+nTurbines] = 1.

                k += 1

        J = {}

        J['separation', 'turbineX'] = dS[:, :nTurbines]
        J['separation', 'turbineY'] = dS[:, nTurbines:]
        # print J
        return J


class MUX(Component):
    """ Connect input elements into a single array  """

    def __init__(self, nElements):

        super(MUX, self).__init__()

        # set finite difference options (fd used for testing only)
        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-5
        self.fd_options['step_type'] = 'relative'

        # define inputs
        for i in range(0, nElements):
            self.add_param('input%i' % i, val=0.0, desc='scalar input')

        # define output
        self.add_output('Array', np.zeros(nElements), desc='ndArray of all the scalar inputs')

        self.nElements = nElements

    def solve_nonlinear(self, params, unknowns, resids):

        # print 'in MUX'

        # assign input values to the output array
        for i in range(0, self.nElements):
            exec("unknowns['Array'][%i] = params['input%i']" % (i, i))

        # print unknowns['Array']

    def linearize(self, params, unknowns, resids):

        dArray_dInput = np.zeros(self.nElements)

        J = {}

        for i in range(0, self.nElements):
            dArray_dInput[i] = 1.0
            J['Array', 'input%i' % i] = np.array(dArray_dInput)
            dArray_dInput[i] = 0.0

        return J


class DeMUX(Component):
    """ split a given array into separate elements """

    def __init__(self, nElements):

        super(DeMUX, self).__init__()

        # set finite difference options (fd used for testing only)
        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-5
        self.fd_options['step_type'] = 'relative'

        # define input
        self.add_param('Array', np.zeros(nElements), desc='ndArray of scalars')

        # define outputs
        for i in range(0, nElements):
            self.add_output('output%i' % i, val=0.0, desc='scalar output')
        # print 'demux elements: ', nElements
        self.nElements = nElements

    def solve_nonlinear(self, params, unknowns, resids):

        # print 'in MUX'

        # assign input values to the output array
        for i in range(0, self.nElements):
            exec("unknowns['output%i'] = params['Array'][%i]" % (i, i))

        # print unknowns['Array']

    def linearize(self, params, unknowns, resids):

        doutput_dArray = np.eye(self.nElements)

        J = {}

        for i in range(0, self.nElements):
            # print doutput_dArray
            J['output%i' % i, 'Array'] = np.reshape(doutput_dArray[i, :], (1, self.nElements))

        return J


class DeMUX2D(Component):
    """ split a given 2D array into individual rows """
    # TODO FINISH THIS COMPONENT
    def __init__(self, ArrayShape):

        super(DeMUX2D, self).__init__()

        # define input
        self.add_param('Array', np.zeros(ArrayShape), desc='ndArray of scalars')

        # define outputs
        for i in range(0, ArrayShape[0]):
            self.add_output('output%i' % i, np.zeros(ArrayShape[1]), desc='scalar output')
        # print 'demux rows: ', ArrayShape[0]
        self.ArrayShape = ArrayShape

    def solve_nonlinear(self, params, unknowns, resids):

        # print 'in MUX'

        # assign input values to the output array
        for i in range(0, self.ArrayShape[0]):
            exec("unknowns['output%i[:]'] = params['Array'][%i, :]" % (i, i))

        # print unknowns['Array']

    def linearize(self, params, unknowns, resids):

        doutput_dArray = np.eye(self.ArrayShape)

        J = {}

        for i in range(0, self.ArrayShape[0]):
            # print doutput_dArray
            J['output%i' % i, 'Array'] = np.reshape(doutput_dArray[i, :], (1, 2))

        return J


## ---- if you know wind speed to power and thrust, you can use these tools ----------------
class CPCT_Interpolate_Gradients(Component):

    def __init__(self, nTurbines, direction_id=0, datasize=0):

        super(CPCT_Interpolate_Gradients, self).__init__()

        # set finite difference options (fd used for testing only)
        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-5
        self.fd_options['step_type'] = 'relative'

        self.nTurbines = nTurbines
        self.direction_id = direction_id
        self.datasize = datasize

        # add inputs and outputs
        self.add_param('yaw', np.zeros(nTurbines), desc='yaw error', units='deg')
        self.add_param('velocitiesTurbines%i' % direction_id, np.zeros(nTurbines), units='m/s', desc='hub height wind speed') # Uhub
        self.add_output('Cp_out', np.zeros(nTurbines))
        self.add_output('Ct_out', np.zeros(nTurbines))

        # add variable trees
        self.add_param('params:pP', 3.0, pass_by_obj=True)
        self.add_param('params:windSpeedToCPCT:wind_speed', np.zeros(datasize), units='m/s',
                       desc='range of wind speeds', pass_by_obj=True)
        self.add_param('params:windSpeedToCPCT:CP', np.zeros(datasize), iotype='out',
                       desc='power coefficients', pass_by_obj=True)
        self.add_param('params:windSpeedToCPCT:CT', np.zeros(datasize), iotype='out',
                       desc='thrust coefficients', pass_by_obj=True)

    def solve_nonlinear(self, params, unknowns, resids):
        direction_id = self.direction_id
        pP = self.params['params:pP']
        wind_speed_ax = np.cos(self.params['yaw']*np.pi/180.0)**(pP/3.0)*self.params['velocitiesTurbines%i' % direction_id]
        # use interpolation on precalculated CP-CT curve
        wind_speed_ax = np.maximum(wind_speed_ax, self.params['params:windSpeedToCPCT:wind_speed'][0])
        wind_speed_ax = np.minimum(wind_speed_ax, self.params['params:windSpeedToCPCT:wind_speed'][-1])
        self.unknowns['Cp_out'] = interp(wind_speed_ax, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CP'])
        self.unknowns['Ct_out'] = interp(wind_speed_ax, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CT'])

        # normalize on incoming wind speed to correct coefficients for yaw
        self.unknowns['Cp_out'] = self.unknowns['Cp_out'] * np.cos(self.params['yaw']*np.pi/180.0)**pP
        self.unknowns['Ct_out'] = self.unknowns['Ct_out'] * np.cos(self.params['yaw']*np.pi/180.0)**2
        # print 'in CPCT interp, wind_speed_hub = ', self.params['velocitiesTurbines%i' % direction_id]

    def linearize(self, params, unknowns, resids):  # standard central differencing
        # set step size for finite differencing
        h = 1e-6
        direction_id = self.direction_id

        # calculate upper and lower function values
        wind_speed_ax_high_yaw = np.cos((self.params['yaw']+h)*np.pi/180.0)**(self.params['params:pP']/3.0)*self.params['velocitiesTurbines%i' % direction_id]
        wind_speed_ax_low_yaw = np.cos((self.params['yaw']-h)*np.pi/180.0)**(self.params['params:pP']/3.0)*self.params['velocitiesTurbines%i' % direction_id]
        wind_speed_ax_high_wind = np.cos(self.params['yaw']*np.pi/180.0)**(self.params['params:pP']/3.0)*(self.params['velocitiesTurbines%i' % direction_id]+h)
        wind_speed_ax_low_wind = np.cos(self.params['yaw']*np.pi/180.0)**(self.params['params:pP']/3.0)*(self.params['velocitiesTurbines%i' % direction_id]-h)

        # use interpolation on precalculated CP-CT curve
        wind_speed_ax_high_yaw = np.maximum(wind_speed_ax_high_yaw, self.params['params:windSpeedToCPCT:wind_speed'][0])
        wind_speed_ax_low_yaw = np.maximum(wind_speed_ax_low_yaw, self.params['params:windSpeedToCPCT:wind_speed'][0])
        wind_speed_ax_high_wind = np.maximum(wind_speed_ax_high_wind, self.params['params:windSpeedToCPCT:wind_speed'][0])
        wind_speed_ax_low_wind = np.maximum(wind_speed_ax_low_wind, self.params['params:windSpeedToCPCT:wind_speed'][0])

        wind_speed_ax_high_yaw = np.minimum(wind_speed_ax_high_yaw, self.params['params:windSpeedToCPCT:wind_speed'][-1])
        wind_speed_ax_low_yaw = np.minimum(wind_speed_ax_low_yaw, self.params['params:windSpeedToCPCT:wind_speed'][-1])
        wind_speed_ax_high_wind = np.minimum(wind_speed_ax_high_wind, self.params['params:windSpeedToCPCT:wind_speed'][-1])
        wind_speed_ax_low_wind = np.minimum(wind_speed_ax_low_wind, self.params['params:windSpeedToCPCT:wind_speed'][-1])

        CP_high_yaw = interp(wind_speed_ax_high_yaw, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CP'])
        CP_low_yaw = interp(wind_speed_ax_low_yaw, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CP'])
        CP_high_wind = interp(wind_speed_ax_high_wind, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CP'])
        CP_low_wind = interp(wind_speed_ax_low_wind, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CP'])

        CT_high_yaw = interp(wind_speed_ax_high_yaw, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CT'])
        CT_low_yaw = interp(wind_speed_ax_low_yaw, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CT'])
        CT_high_wind = interp(wind_speed_ax_high_wind, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CT'])
        CT_low_wind = interp(wind_speed_ax_low_wind, self.params['params:windSpeedToCPCT:wind_speed'], self.params['params:windSpeedToCPCT:CT'])

        # normalize on incoming wind speed to correct coefficients for yaw
        CP_high_yaw = CP_high_yaw * np.cos((self.params['yaw']+h)*np.pi/180.0)**self.params['params:pP']
        CP_low_yaw = CP_low_yaw * np.cos((self.params['yaw']-h)*np.pi/180.0)**self.params['params:pP']
        CP_high_wind = CP_high_wind * np.cos((self.params['yaw'])*np.pi/180.0)**self.params['params:pP']
        CP_low_wind = CP_low_wind * np.cos((self.params['yaw'])*np.pi/180.0)**self.params['params:pP']

        CT_high_yaw = CT_high_yaw * np.cos((self.params['yaw']+h)*np.pi/180.0)**2
        CT_low_yaw = CT_low_yaw * np.cos((self.params['yaw']-h)*np.pi/180.0)**2
        CT_high_wind = CT_high_wind * np.cos((self.params['yaw'])*np.pi/180.0)**2
        CT_low_wind = CT_low_wind * np.cos((self.params['yaw'])*np.pi/180.0)**2

        # compute derivative via central differencing and arrange in sub-matrices of the Jacobian
        dCP_dyaw = np.eye(self.nTurbines)*(CP_high_yaw-CP_low_yaw)/(2.0*h)
        dCP_dwind = np.eye(self.nTurbines)*(CP_high_wind-CP_low_wind)/(2.0*h)
        dCT_dyaw = np.eye(self.nTurbines)*(CT_high_yaw-CT_low_yaw)/(2.0*h)
        dCT_dwind = np.eye(self.nTurbines)*(CT_high_wind-CT_low_wind)/(2.0*h)

        # compile Jacobian dict from sub-matrices
        J = {}
        J['Cp_out', 'yaw'] = dCP_dyaw
        J['Cp_out', 'velocitiesTurbines%i' % direction_id] = dCP_dwind
        J['Ct_out', 'yaw'] = dCT_dyaw
        J['Ct_out', 'velocitiesTurbines%i' % direction_id] = dCT_dwind

        return J


if __name__ == "__main__":

    # top = Problem()
    #
    # root = top.root = Group()
    #
    # root.add('p1', IndepVarComp('x', np.array([1.0, 1.0])))
    # root.add('p2', IndepVarComp('y', np.array([0.75, 0.25])))
    # root.add('p', WindFarmAEP(nDirections=2))
    #
    # root.connect('p1.x', 'p.power_directions')
    # root.connect('p2.y', 'p.windrose_frequencies')
    #
    # top.setup()
    # top.run()
    #
    # # should return 8760.0
    # print(root.p.unknowns['AEP'])
    # top.check_partial_derivatives()

    # top = Problem()
    #
    # root = top.root = Group()
    #
    # root.add('p1', IndepVarComp('x', 1.0))
    # root.add('p2', IndepVarComp('y', 2.0))
    # root.add('p', MUX(nElements=2))
    #
    # root.connect('p1.x', 'p.input0')
    # root.connect('p2.y', 'p.input1')
    #
    # top.setup()
    # top.run()
    #
    # # should return 8760.0
    # print(root.p.unknowns['Array'])
    # top.check_partial_derivatives()

    # top = Problem()
    #
    # root = top.root = Group()
    #
    # root.add('p1', IndepVarComp('x', np.zeros(2)))
    # root.add('p', DeMUX(nElements=2))
    #
    # root.connect('p1.x', 'p.Array')
    #
    # top.setup()
    # top.run()
    #
    # # should return 8760.0
    # print(root.p.unknowns['output0'])
    # print(root.p.unknowns['output1'])
    # top.check_partial_derivatives()

    top = Problem()

    root = top.root = Group()

    root.add('p1', IndepVarComp('x', np.array([0, 3])))
    root.add('p2', IndepVarComp('y', np.array([1, 0])))
    root.add('p', SpacingComp(nTurbines=2))

    root.connect('p1.x', 'p.turbineX')
    root.connect('p2.y', 'p.turbineY')

    top.setup()
    top.run()

    # print(root.p.unknowns['output0'])
    # print(root.p.unknowns['output1'])
    top.check_partial_derivatives()

