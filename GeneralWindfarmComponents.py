from openmdao.api import Component, Group, Problem, IndepVarComp

import numpy as np


class WindFrame(Component):
    """ Calculates the locations of each turbine in the wind direction reference frame """

    def __init__(self, nTurbines, resolution):

        # print 'entering windframe __init__ - analytic'

        super(WindFrame, self).__init__()

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

        # print turbineX, turbineY

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

        # Explicitly size input arrays
        self.add_param('Ct_in', val=np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_param('Cp_in', val=np.zeros(nTurbines), desc='power coefficient for all turbines')
        self.add_param('yaw', val=np.zeros(nTurbines), desc='yaw of each turbine')

        # Explicitly size output arrays
        self.add_output('Ct_out', val=np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_output('Cp_out', val=np.zeros(nTurbines), desc='power coefficient for all turbines')

        # parameters since var trees are not supports
        self.add_param('params:pP', 1.88)
        self.add_param('params:CTcorrected', True,
                       desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('params:CPcorrected', True,
                       desc='CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)')
        self.add_param('floris_params:FLORISoriginal', False,
                       desc='override all parameters and use FLORIS as original in first Wind Energy paper')

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

        # print 'In AEP, AEP %s' % unknowns['AEP']

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

        # Explicitly size input arrays
        self.add_param('turbineX', np.zeros(nTurbines),
                       desc='x coordinates of turbines in wind dir. ref. frame')
        self.add_param('turbineY', np.zeros(nTurbines),
                       desc='y coordinates of turbines in wind dir. ref. frame')

        # Explicitly size output array
        self.add_output('separation', np.zeros((nTurbines-1.)*nTurbines/2.),
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
        nTurbines = turbineX.size
        dS = np.zeros(((nTurbines-1.)*nTurbines/2., 2*nTurbines))
        k = 0
        # print 'in dist_const, turbineX = ', turbineX
        # print 'in dist_const, turbineY = ', turbineY

        for i in range(0, nTurbines):
            for j in range(i+1, nTurbines):
                dS[k, j] = (turbineX[j]-turbineX[i])*((turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2)**(-0.5)
                dS[k, i] = (turbineX[i]-turbineX[j])*((turbineX[j]-turbineX[i])**2+(turbineY[j]-turbineY[i])**2)**(-0.5)
                dS[k, j+nTurbines] = (turbineY[j]-turbineY[i])*((turbineX[j]-turbineX[i])**2 +
                                                                (turbineY[j]-turbineY[i])**2)**(-0.5)
                dS[k, i+nTurbines] = (turbineY[i]-turbineY[j])*((turbineX[j]-turbineX[i])**2 +
                                                                (turbineY[j]-turbineY[i])**2)**(-0.5)
                k += 1

        J = {}

        J['separation', 'turbineX'] = dS[:, 0:nTurbines]
        J['separation', 'turbineY'] = dS[:, nTurbines:nTurbines*nTurbines]
        print J
        return J


class MUX(Component):
    """ Estimate the AEP based on power production for each direction and weighted by wind direction frequency  """

    def __init__(self, nElements):

        super(MUX, self).__init__()

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

    top = Problem()

    root = top.root = Group()

    root.add('p1', IndepVarComp('x', np.zeros(2)))
    root.add('p', DeMUX(nElements=2))

    root.connect('p1.x', 'p.Array')

    top.setup()
    top.run()

    # should return 8760.0
    print(root.p.unknowns['output0'])
    print(root.p.unknowns['output1'])
    top.check_partial_derivatives()


