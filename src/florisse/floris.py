import numpy as np

from openmdao.api import Group, Component, Problem, IndepVarComp, ParamComp, ParallelGroup, NLGaussSeidel, LinearGaussSeidel, ScipyGMRES

from GeneralWindFarmComponents import WindFrame, AdjustCtCpYaw, MUX, WindFarmAEP, DeMUX, CPCT_Interpolate_Gradients
from Parameters import FLORISParameters
import _floris
import _florisDiscontinuous


# Components of FLORIS - for full model use FLORIS(Group)
class floris_wcent_wdiam(Component):

    def __init__(self, nTurbines, direction_id=0, differentiable=True):

        super(floris_wcent_wdiam, self).__init__()
        
        self.direction_id = direction_id
        self.differentiable = differentiable

        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-5
        self.fd_options['step_type'] = 'relative'

        if not differentiable:
            self.fd_options['force_fd'] = True
            self.fd_options['form'] = 'forward'

        # print 'entering wcent_wdiam __init__ - Tapenade'

        # input arrays
        self.add_param('turbineXw', np.zeros(nTurbines), units='m', desc='x coordinates of turbines in wind dir. ref. frame')
        self.add_param('turbineYw', np.zeros(nTurbines), units='m', desc='y coordinates of turbines in wind dir. ref. frame')
        self.add_param('yaw%i' % direction_id, np.zeros(nTurbines), units='deg', desc='yaw of each turbine')
        self.add_param('rotorDiameter', np.zeros(nTurbines), units='m', desc='rotor diameter of each turbine')
        self.add_param('Ct', np.zeros(nTurbines), desc='thrust coefficient of each turbine')

        # output arrays
        self.add_output('wakeCentersYT', np.zeros(nTurbines*nTurbines), units='m', desc='wake center y position at each turbine')
        self.add_output('wakeDiametersT', np.zeros(3*nTurbines*nTurbines), units='m', desc='wake diameter of each zone of each wake at each turbine')

        # FLORIS parameters
        self.add_param('floris_params:pP', 1.88, pass_by_obj=True)
        self.add_param('floris_params:ke', 0.065, pass_by_obj=True)
        self.add_param('floris_params:keCorrArray', 0.0, pass_by_obj=True)
        self.add_param('floris_params:keCorrCT', 0.0, pass_by_obj=True)
        self.add_param('floris_params:Region2CT', 4.0*(1.0/3.0)*(1.0-(1.0/3.0)), pass_by_obj=True)
        self.add_param('floris_params:kd', 0.15, pass_by_obj=True)
        self.add_param('floris_params:me', np.array([-0.5, 0.22, 1.0]), pass_by_obj=True)
        self.add_param('floris_params:initialWakeDisplacement', -4.5, pass_by_obj=True)
        self.add_param('floris_params:initialWakeAngle', 3.0, pass_by_obj=True)
        self.add_param('floris_params:baselineCT', 4./3.*(1.-1./3.), pass_by_obj=True)
        self.add_param('floris_params:keCorrTI', 0.0, pass_by_obj=True)
        self.add_param('floris_params:baselineTI', 0.045, pass_by_obj=True)
        self.add_param('floris_params:keCorrHR', 0.0, pass_by_obj=True) # neutral, with heating rate 0, is baseline
        self.add_param('floris_params:keCorrHRTI', 0.0, pass_by_obj=True)
        self.add_param('floris_params:keSaturation', 0.0)
        self.add_param('floris_params:kdCorrYawDirection', 0.0, pass_by_obj=True)
        self.add_param('floris_params:MU', np.array([0.5, 1.0, 10]), pass_by_obj=True)
        self.add_param('floris_params:CTcorrected', False, desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)', pass_by_obj=True)
        self.add_param('floris_params:CPcorrected', False, desc = 'CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)', pass_by_obj=True)
        self.add_param('floris_params:axialIndProvided', False, desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)', pass_by_obj=True)
        self.add_param('floris_params:useWakeAngle', True, pass_by_obj=True)
        self.add_param('floris_params:bd', -0.01, pass_by_obj=True)
        self.add_param('floris_params:useaUbU', False, pass_by_obj=True)
        self.add_param('floris_params:aU', 5.0, units='deg', pass_by_obj=True)
        self.add_param('floris_params:bU', 1.66, pass_by_obj=True)
        self.add_param('floris_params:adjustInitialWakeDiamToYaw', True, pass_by_obj=True)
        self.add_param('floris_params:FLORISoriginal', True, desc='override all parameters and use FLORIS as original in first Wind Energy paper', pass_by_obj=True)

    def solve_nonlinear(self, params, unknowns, resids):

        # print 'entering wcent_wdiam - tapenade'

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
        yaw_deg = params['yaw%i' % self.direction_id]
        # print yaw_deg, Ct
        if self.differentiable:
            wakeCentersYT_vec, wakeDiametersT_vec = _floris.floris_wcent_wdiam(kd, initialWakeDisplacement, \
                                  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                                  rotorDiameter, me, bd, useWakeAngle, adjustInitialWakeDiamToYaw)
        else:
            wakeCentersYT_vec, wakeDiametersT_vec = _florisDiscontinuous.floris_wcent_wdiam(kd, initialWakeDisplacement, \
                                  initialWakeAngle, ke, keCorrCT, Region2CT, yaw_deg, Ct, turbineXw, turbineYw, \
                                  rotorDiameter, me, bd, useWakeAngle, adjustInitialWakeDiamToYaw)

        # Outputs in vector form so they can be used in Jacobian creation
        unknowns['wakeCentersYT'] = wakeCentersYT_vec
        unknowns['wakeDiametersT'] = wakeDiametersT_vec

        # print 'yaw: ', yaw_deg

    def linearize(self, params, unknowns, resids):

        # # # print 'entering wcen wdiam linearize'
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
        yaw_deg = params['yaw%i' % self.direction_id]

        # number of turbines
        nTurbines = np.size(turbineXw)

        # number of directions being differentiated in the Jacobian
        nbdirs = nTurbines*nTurbines

        # input arrays to direct differentiation
        wakeCentersYT_vecb = np.eye(nbdirs, nTurbines*nTurbines)
        wakeDiametersT_vecb = np.zeros((nbdirs, 3*nTurbines*nTurbines))


        # initialize linearize dict
        J = {}

        # function call to extract gradients of wakeCentersYT w.r.t. all design vars
        yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb = \
            _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement, initialWakeAngle, ke, keCorrCT, Region2CT,
                                          yaw_deg, Ct, turbineXw, turbineYw, rotorDiameter, me, bd, useWakeAngle,
                                          adjustInitialWakeDiamToYaw, wakeCentersYT_vecb, wakeDiametersT_vecb)

        # print 'here', turbineXwb, yawb.shape, Ctb.shape, rotorDiameterb.shape

        # construct Jacobian of wakeCentersYT
        J['wakeCentersYT', 'yaw%i' % self.direction_id] = yawb
        J['wakeCentersYT', 'Ct'] = Ctb
        J['wakeCentersYT', 'turbineXw'] = turbineXwb
        J['wakeCentersYT', 'turbineYw'] = turbineYwb
        J['wakeCentersYT', 'rotorDiameter'] = rotorDiameterb

        # number of directions being differentiated in the Jacobian
        nbdirs = 3*nTurbines*nTurbines

        # input arrays to direct differentiation
        # wakeCentersYT_vecb[:, :] = 0.0
        wakeCentersYT_vecb = np.zeros((nbdirs, nTurbines*nTurbines))
        wakeDiametersT_vecb = np.eye(nbdirs, nbdirs)


        # function call to extract gradients of wakeDiametersT w.r.t. all design vars
        yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb = \
            _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement, initialWakeAngle, ke, keCorrCT, Region2CT,
                                          yaw_deg, Ct, turbineXw, turbineYw, rotorDiameter, me, bd, useWakeAngle,
                                          adjustInitialWakeDiamToYaw, wakeCentersYT_vecb, wakeDiametersT_vecb)

        # construct Jacobian of wakeDiametersT
        J['wakeDiametersT', 'yaw%i' % self.direction_id] = yawb
        J['wakeDiametersT', 'Ct'] = Ctb
        J['wakeDiametersT', 'turbineXw'] = turbineXwb
        J['wakeDiametersT', 'turbineYw'] = turbineYwb
        J['wakeDiametersT', 'rotorDiameter'] = rotorDiameterb

        return J


class floris_overlap(Component):
    """ Calculates the overlap between each turbine rotor and the existing turbine wakes """

    def __init__(self, nTurbines, differentiable=True):

        super(floris_overlap, self).__init__()

        self.differentiable = differentiable

        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-6
        self.fd_options['step_type'] = 'relative'

        if not differentiable:
            self.fd_options['force_fd'] = True
            self.fd_options['form'] = 'forward'

        # print 'entering overlap __init__ - Tapenade'

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

        if differentiable:
            self.add_output('cosFac', np.ones(3*nTurbines*nTurbines),
                            desc='cosine factor similar to Jensen 1983')

        # floris parameters
        self.add_param('floris_params:cos_spread', val=3.0, pass_by_obj=True,
                       desc='spread of cosine smoothing factor (percent of sum of wake and rotor radii)')

        # etc
        self.nTurbines = nTurbines

    def solve_nonlinear(self, params, unknowns, resids):

        # print 'entering overlap - Tapenade'

        # call to fortran code to obtain relative wake overlap values
        # print params['turbineXw'], params['turbineYw'], params['rotorDiameter']
        if self.differentiable:
            wakeOverlapTRel_vec, cosFac_vec = _floris.floris_overlap(params['turbineXw'], params['turbineYw'],
                                                                     params['rotorDiameter'], params['wakeDiametersT'],
                                                                     params['wakeCentersYT'], params['floris_params:cos_spread'])
             # pass results to self in the form of a vector for use in Jacobian creation
            unknowns['wakeOverlapTRel'] = wakeOverlapTRel_vec
            unknowns['cosFac'] = cosFac_vec
        else:
            wakeOverlapTRel_vec = _florisDiscontinuous.floris_overlap(params['turbineXw'], params['turbineYw'],
                                                         params['rotorDiameter'], params['wakeDiametersT'],
                                                         params['wakeCentersYT'])
            # pass results to self in the form of a vector for use in Jacobian creation
            unknowns['wakeOverlapTRel'] = wakeOverlapTRel_vec

    def linearize(self, params, unknowns, resids):
        # print 'entering overlap linearize'
        # number of turbines
        nTurbines = self.nTurbines

        # number of directions being differentiated
        nbdirs = 3*nTurbines*nTurbines

        # input array to direct differentiation
        wakeOverlapTRel_vecb = np.eye(nbdirs, 3*nTurbines*nTurbines)
        cosFac_vecb = np.zeros((nbdirs, 3*nTurbines*nTurbines))
        # print params['rotorDiameter'], cosFac_vecb.shape, wakeOverlapTRel_vecb.shape
        # function call to fortran to obtain gradients
        # print params['turbineXw'], params['turbineYw'], params['rotorDiameter']
        turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb, _, _ \
            = _floris.floris_overlap_bv(params['turbineXw'], params['turbineYw'], params['rotorDiameter'],
                                        params['wakeDiametersT'], params['wakeCentersYT'],
                                        params['floris_params:cos_spread'], wakeOverlapTRel_vecb, cosFac_vecb)

        J = {}
        # construct Jacobian of floris_overlap
        J['wakeOverlapTRel', 'turbineYw'] = turbineYwb
        J['wakeOverlapTRel', 'rotorDiameter'] = rotorDiameterb
        J['wakeOverlapTRel', 'wakeDiametersT'] = wakeDiametersT_vecb
        J['wakeOverlapTRel', 'wakeCentersYT'] = wakeCentersYT_vecb

         # input array to direct differentiation
        wakeOverlapTRel_vecb = np.zeros((nbdirs, 3*nTurbines*nTurbines))
        cosFac_vecb = np.eye(nbdirs, 3*nTurbines*nTurbines)

        # function call to fortran to obtain gradients
        turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb, _, _ \
            = _floris.floris_overlap_bv(params['turbineXw'], params['turbineYw'], params['rotorDiameter'],
                                        params['wakeDiametersT'], params['wakeCentersYT'],
                                        params['floris_params:cos_spread'], wakeOverlapTRel_vecb, cosFac_vecb)

        # construct Jacobian of floris_overlap
        J['cosFac', 'turbineYw'] = turbineYwb
        J['cosFac', 'rotorDiameter'] = rotorDiameterb
        J['cosFac', 'wakeDiametersT'] = wakeDiametersT_vecb
        J['cosFac', 'wakeCentersYT'] = wakeCentersYT_vecb

        return J


class floris_power(Component):
    """ Calculates the turbine power and effective wind speed for each turbine """

    def __init__(self, nTurbines, direction_id=0, differentiable=True):

        super(floris_power, self).__init__()

        self.differentiable = differentiable

        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-6
        self.fd_options['step_type'] = 'relative'

        if not differentiable:
            self.fd_options['force_fd'] = True
            self.fd_options['form'] = 'forward'

        # print 'entering power __init__ - Tapenade'
        self.nTurbines = nTurbines
        self.direction_id = direction_id

        # inputs
        self.add_param('wind_speed', 8.0, units='m/s', desc='free stream wind velocity')
        self.add_param('air_density', 1.1716, units='kg/(m*m*m)', desc='air density in free stream')
        self.add_param('rotorDiameter', np.zeros(nTurbines)+126.4, units='m', desc='rotor diameters of all turbine')
        self.add_param('axialInduction', np.zeros(nTurbines)+1./3., desc='axial induction of all turbines')
        self.add_param('Ct', np.zeros(nTurbines)+4.0*(1./3.)*(1.0-(1./3.)), desc='Thrust coefficient for all turbines')
        self.add_param('Cp', np.zeros(nTurbines)+0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2), desc='power coefficient for all turbines')
        self.add_param('generator_efficiency', np.zeros(nTurbines)+0.944, desc='generator efficiency of all turbines')
        self.add_param('turbineXw', np.zeros(nTurbines), units='m',
                       desc='X positions of turbines in the wind direction reference frame')
        self.add_param('yaw%i' % direction_id, np.zeros(nTurbines), units='deg',
                       desc='yaw angle of turbines wrt the wind direction')
        self.add_param('wakeCentersYT',  np.zeros(nTurbines*nTurbines), units='m',
                       desc='centers of the wakes at each turbine')
        self.add_param('wakeDiametersT', np.zeros(3*nTurbines*nTurbines), units='m',
                       desc='diameters of each of the wake zones for each of the wakes at each turbine')
        self.add_param('wakeOverlapTRel', np.zeros(3*nTurbines*nTurbines),
                       desc='ratios of wake overlap area per zone to rotor area')
        if differentiable:
            self.add_param('cosFac', np.zeros(3*nTurbines*nTurbines),
                           desc='cosine factor similar to Jensen 1983')

        # outputs
        # self.add_output('velocitiesTurbines%i' % direction_id, np.zeros(nTurbines), units='m/s',
        #                desc='effective hub velocity for each turbine')
        self.add_output('wt_power%i' % direction_id, np.zeros(nTurbines), units='kW', desc='power output of each turbine')
        # output
        self.add_output('power%i' % direction_id, 0.0, units='kW', desc='total power output of the wind farm')

        # add state for solver
        self.add_output('velocitiesTurbines%i' % direction_id, np.zeros(nTurbines), units='m/s',
                       desc='effective hub velocity for each turbine')

        # connect floris_params
        self.add_param('floris_params:pP', 1.88, pass_by_obj=True)
        self.add_param('floris_params:ke', 0.065, pass_by_obj=True)
        self.add_param('floris_params:keCorrArray', 0.0, pass_by_obj=True)
        self.add_param('floris_params:keCorrCT', 0.0, pass_by_obj=True)
        self.add_param('floris_params:Region2CT', 4.0*(1.0/3.0)*(1.0-(1.0/3.0)), pass_by_obj=True)
        self.add_param('floris_params:kd', 0.15, pass_by_obj=True)
        self.add_param('floris_params:me', np.array([-0.5, 0.22, 1.0]), pass_by_obj=True)
        self.add_param('floris_params:initialWakeDisplacement', -4.5, pass_by_obj=True)
        self.add_param('floris_params:initialWakeAngle', 3.0, pass_by_obj=True)
        self.add_param('floris_params:baselineCT', 4./3.*(1.-1./3.), pass_by_obj=True)
        self.add_param('floris_params:keCorrTI', 0.0, pass_by_obj=True)
        self.add_param('floris_params:baselineTI', 0.045, pass_by_obj=True)
        self.add_param('floris_params:keCorrHR', 0.0, pass_by_obj=True) # neutral, with heating rate 0, is baseline
        self.add_param('floris_params:keCorrHRTI', 0.0, pass_by_obj=True)
        self.add_param('floris_params:keSaturation', 0.0, pass_by_obj=True)
        self.add_param('floris_params:kdCorrYawDirection', 0.0, pass_by_obj=True)
        self.add_param('floris_params:MU', np.array([0.5, 1.0, 10]), pass_by_obj=True)
        self.add_param('floris_params:CTcorrected', False,
                       desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)', pass_by_obj=True)
        self.add_param('floris_params:CPcorrected', False,
                       desc='CP factor already corrected by CCBlade calculation '
                            '(assumed with approximately factor cos(yaw)^3)', pass_by_obj=True)
        self.add_param('floris_params:axialIndProvided', False,
                       desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)', pass_by_obj=True)
        self.add_param('floris_params:useWakeAngle', True, pass_by_obj=True)
        self.add_param('floris_params:bd', -0.01, pass_by_obj=True)
        self.add_param('floris_params:useaUbU', False, pass_by_obj=True)
        self.add_param('floris_params:aU', 5.0, units='deg', pass_by_obj=True)
        self.add_param('floris_params:bU', 1.66, pass_by_obj=True)
        self.add_param('floris_params:adjustInitialWakeDiamToYaw', True, pass_by_obj=True)
        self.add_param('floris_params:FLORISoriginal', True,
                       desc='override all parameters and use FLORIS as original in first Wind Energy paper', pass_by_obj=True)

    def solve_nonlinear(self, params, unknowns, resids):
        # print 'entering power - tapenade'

        # reassign input variables
        wakeOverlapTRel_v = params['wakeOverlapTRel']
        if self.differentiable:
            cosFac_v = params['cosFac']
        Region2CT = params['floris_params:Region2CT']
        Ct = params['Ct']
        Vinf = params['wind_speed']
        turbineXw = params['turbineXw']
        axialInduction = params['axialInduction']
        rotorDiameter = params['rotorDiameter']
        rho = params['air_density']
        generator_efficiency = params['generator_efficiency']
        Cp = params['Cp']
        yaw = params['yaw%i' % self.direction_id]

        # print 'wake OL in power', wakeOverlapTRel_v

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

        # print 'ke, ME, aI, bU: ', ke, MU, aU, bU

        axialIndProvided = params['floris_params:axialIndProvided']

        # how far in front of turbines to use overlap power calculations (in rotor diameters). This must match the
        # value used in floris_wcent_wdiam (hardcoded in fortran as 1)
        if self.differentiable:
            # TODO hard code this parameter in the fortran code and remove the specifier from all functions of this component
            p_near0 = 1.0

            # pass p_near0 to self for use in gradient calculations
            self.p_near0 = p_near0

            # call to fortran code to obtain output values
            velocitiesTurbines, wt_power, power = _floris.floris_power(wakeOverlapTRel_v, cosFac_v, Ct, axialInduction,
                                                                axialIndProvided, useaUbU, keCorrCT, Region2CT, ke,
                                                                Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter,
                                                                MU, rho, aU, bU, Cp, generator_efficiency)
        else:
            velocitiesTurbines, wt_power, power = _florisDiscontinuous.floris_power(wakeOverlapTRel_v, Ct, axialInduction, \
                                                                axialIndProvided, useaUbU, keCorrCT, Region2CT, ke, \
                                                                Vinf, keCorrArray, turbineXw, yaw, rotorDiameter, MU, \
                                                                rho, aU, bU, Cp, generator_efficiency)

        # pass outputs to self
        unknowns['velocitiesTurbines%i' % self.direction_id] = velocitiesTurbines
        unknowns['wt_power%i' % self.direction_id] = wt_power
        unknowns['power%i' % self.direction_id] = power

        # print 'velocitiesTurbines: ', velocitiesTurbines
        # print 'wt_power: ', wt_power
        # print 'power: ', power

    def apply_nonlinear(self, params, unknowns, resids):
        # print 'entering power - tapenade'

        # reassign input variables
        wakeOverlapTRel_v = params['wakeOverlapTRel']
        if self.differentiable:
            cosFac_v = params['cosFac']
        Region2CT = params['floris_params:Region2CT']
        Ct = params['Ct']
        Vinf = params['wind_speed']
        turbineXw = params['turbineXw']
        axialInduction = params['axialInduction']
        rotorDiameter = params['rotorDiameter']
        rho = params['air_density']
        generator_efficiency = params['generator_efficiency']
        Cp = params['Cp']
        yaw = params['yaw%i' % self.direction_id]

        # print 'wake OL in power', wakeOverlapTRel_v

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

        # print 'ke, ME, aI, bU: ', ke, MU, aU, bU

        axialIndProvided = params['floris_params:axialIndProvided']

        # how far in front of turbines to use overlap power calculations (in rotor diameters). This must match the
        # value used in floris_wcent_wdiam (hardcoded in fortran as 1)
        # TODO hard code this parameter in the fortran code and remove the specifier from all functions of this component
        p_near0 = 1.0

        # pass p_near0 to self for use in gradient calculations
        self.p_near0 = p_near0

        if self.differentiable:
            # TODO hard code this parameter in the fortran code and remove the specifier from all functions of this component
            p_near0 = 1.0

            # pass p_near0 to self for use in gradient calculations
            self.p_near0 = p_near0

            # call to fortran code to obtain output values
            velocitiesTurbines, wt_power, power = _floris.floris_power(wakeOverlapTRel_v, cosFac_v, Ct, axialInduction,
                                                                axialIndProvided, useaUbU, keCorrCT, Region2CT, ke,
                                                                Vinf, keCorrArray, turbineXw, yaw, p_near0, rotorDiameter,
                                                                MU, rho, aU, bU, Cp, generator_efficiency)
        else:
            velocitiesTurbines, wt_power, power = _florisDiscontinuous.floris_power(wakeOverlapTRel_v, Ct, axialInduction, \
                                                                axialIndProvided, useaUbU, keCorrCT, Region2CT, ke, \
                                                                Vinf, keCorrArray, turbineXw, yaw, rotorDiameter, MU, \
                                                                rho, aU, bU, Cp, generator_efficiency)
            print 'here'

        # pass outputs to self
        resids['velocitiesTurbines%i' % self.direction_id] = velocitiesTurbines - unknowns['velocitiesTurbines%i' %
                                                                                           self.direction_id]
        resids['wt_power%i' % self.direction_id] = wt_power - unknowns['wt_power%i' % self.direction_id]
        resids['power%i' % self.direction_id] = power - unknowns['power%i' % self.direction_id]
        # print resids['velocitiesTurbines%i' % self.direction_id], unknowns['velocitiesTurbines%i' % self.direction_id]

        # print 'velocitiesTurbines: ', velocitiesTurbines
        # print 'wt_power: ', wt_power
        # print 'power: ', power

    def linearize(self, params, unknowns, resids):
        # print 'entering power linearize'
        # number of turbines
        nTurbines = self.nTurbines
        direction_id = self.direction_id

        # number of directions to differentiate
        nbdirs = nTurbines

        # reassign input variables
        wakeOverlapTRel_v = params['wakeOverlapTRel']
        cosFac = params['cosFac']
        Region2CT = params['floris_params:Region2CT']
        Ct = params['Ct']
        Vinf = params['wind_speed']
        turbineXw = params['turbineXw']
        axialInduction = params['axialInduction']
        rotorDiameter = params['rotorDiameter']
        rho = params['air_density']
        generator_efficiency = params['generator_efficiency']
        Cp = params['Cp']
        yaw = params['yaw%i' % self.direction_id]

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
        # p_near0 = self.p_near0
        p_near0 = 1.0

        # create jacobian dict
        J = {}

        # input arrays to direct differentiation
        velocitiesTurbinesb = np.eye(nbdirs, nTurbines)
        # velocitiesTurbinesb = np.zeros((nbdirs, nTurbines))
        # velocitiesTurbinesb[:, 0] = 1.0
        wt_powerb = np.zeros((nbdirs, nTurbines))
        powerb = np.zeros(nbdirs)

        # call to fortran to obtain gradients of velocitiesTurbines
        wakeOverlapTRel_vb, cosFac_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, cosFac, Ct, axialInduction, axialIndProvided, useaUbU,
                                      keCorrCT, Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0,
                                      rotorDiameter, MU, rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb,
                                      wt_powerb, powerb)
        # print 'dOL', wakeOverlapTRel_vb.shape
        # print 'Ct', Ctb
        # print 'AI', axialInductionb
        # print 'Xw', turbineXwb
        # print 'Y', yawb
        # print 'D', rotorDiameterb
        # print 'Cp', Cpb

        # collect values of the jacobian of velocitiesTurbines
        J['velocitiesTurbines%i' % direction_id, 'wakeOverlapTRel'] = wakeOverlapTRel_vb
        J['velocitiesTurbines%i' % direction_id, 'cosFac'] = cosFac_vb
        J['velocitiesTurbines%i' % direction_id, 'Ct'] = Ctb
        J['velocitiesTurbines%i' % direction_id, 'Cp'] = Cpb
        J['velocitiesTurbines%i' % direction_id, 'axialInduction'] = axialInductionb
        J['velocitiesTurbines%i' % direction_id, 'turbineXw'] = turbineXwb
        J['velocitiesTurbines%i' % direction_id, 'yaw%i' % direction_id] = yawb
        J['velocitiesTurbines%i' % direction_id, 'rotorDiameter'] = rotorDiameterb

        # input arrays to direct differentiation
        velocitiesTurbinesb[:, :] = 0.0
        wt_powerb = np.eye(nbdirs, nTurbines)

        # call to fortran to obtain gradients wt_power
        wakeOverlapTRel_vb, cosFac_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, cosFac, Ct, axialInduction, axialIndProvided, useaUbU,
                                      keCorrCT, Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0,
                                      rotorDiameter, MU, rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb,
                                      wt_powerb, powerb)

        # print wakeOverlapTRel_vb.shape, wt_powerb.shape
        # collect values of the jacobian of wt_power
        J['wt_power%i' % direction_id, 'wakeOverlapTRel'] = wakeOverlapTRel_vb
        J['wt_power%i' % direction_id, 'cosFac'] = cosFac_vb
        J['wt_power%i' % direction_id, 'Ct'] = Ctb
        J['wt_power%i' % direction_id, 'Cp'] = Cpb
        J['wt_power%i' % direction_id, 'axialInduction'] = axialInductionb
        J['wt_power%i' % direction_id, 'turbineXw'] = turbineXwb
        J['wt_power%i' % direction_id, 'yaw%i' % direction_id] = yawb
        J['wt_power%i' % direction_id, 'rotorDiameter'] = rotorDiameterb

        # input arrays to direct differentiation
        wt_powerb[:, :] = 0.0
        powerb[0] = 1.0

        # call to fortran to obtain gradients of power
        wakeOverlapTRel_vb, cosFac_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb, _, _, _ \
            = _floris.floris_power_bv(wakeOverlapTRel_v, cosFac, Ct, axialInduction, axialIndProvided, useaUbU,
                                      keCorrCT, Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0,
                                      rotorDiameter, MU, rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb,
                                      wt_powerb, powerb)

        #print 'OL', wakeOverlapTRel_vb[0, :]#, Ctb[0, :].shape, Cpb[0, :].shape, axialInductionb[0, :], turbineXwb[0, :].shape, yawb[0, :].shape, rotorDiameterb[0, :].shape

        # # print np.array(yawb[:1, :])
        # collect values of the jacobian of wt_power
        # print 'vb', wakeOverlapTRel_vb
        J['power%i' % direction_id, 'wakeOverlapTRel'] = np.array(wakeOverlapTRel_vb[:1, :])
        J['power%i' % direction_id, 'cosFac'] = np.array(cosFac_vb[:1, :])
        J['power%i' % direction_id, 'Ct'] = np.array(Ctb[:1, :])
        J['power%i' % direction_id, 'Cp'] = np.array(Cpb[:1, :])
        J['power%i' % direction_id, 'axialInduction'] = np.array(axialInductionb[:1, :])
        J['power%i' % direction_id, 'turbineXw'] = np.array(turbineXwb[:1, :])
        J['power%i' % direction_id, 'yaw%i' % direction_id] = np.array(yawb[:1, :])
        J['power%i' % direction_id, 'rotorDiameter'] = np.array(rotorDiameterb[:1, :])

        # print 'leaving power linearize'
        # print J
        return J


# Groups using FLORIS
class FLORIS(Group):
    """ Group containing all necessary components of the floris model """

    def __init__(self, nTurbines, resolution, direction_id=0, differentiable=True):
        super(FLORIS, self).__init__()
        self.add('f_1', WindFrame(nTurbines, resolution, differentiable=differentiable), promotes=['*'])
        self.add('f_2', floris_wcent_wdiam(nTurbines, differentiable=differentiable), promotes=['*'])
        self.add('f_3', floris_overlap(nTurbines, differentiable=differentiable), promotes=['*'])
        self.add('f_4', floris_power(nTurbines, direction_id=direction_id, differentiable=differentiable),
                 promotes=['*'])


class DirectionGroupFLORIS(Group):
    """
    Group containing all necessary components for wind plant calculations
    in a single direction
    """

    def __init__(self, nTurbines, resolution=0, direction_id=0, use_rotor_components=False, datasize=0,
                 differentiable=True):
        super(DirectionGroupFLORIS, self).__init__()
        epsilon = 1e-6

        # self.add('p0', IndepVarComp('wind_direction', val=0.0, units='deg'), promotes=['*'])

        # self.add('fp', FLORISParameters(), promotes=['*'])
        if use_rotor_components:
            self.nl_solver = NLGaussSeidel()
            self.nl_solver.options['atol'] = epsilon
            self.nl_solver.options['rtol'] = epsilon;
            # self.ln_solver = LinearGaussSeidel()
            self.ln_solver = ScipyGMRES()
            # self.ln_solver.options['atol'] = epsilon
            # self.ln_solver.options['rtol'] = epsilon
            # self.ln_solver.setup('direction_group%i' % direction_id)
            # self.ln_solver.solve('velocitiesTurbines%i' % direction_id)
            self.add('CtCp', CPCT_Interpolate_Gradients(nTurbines, direction_id=direction_id, datasize=0),
                     promotes=['params:*', 'yaw%i' % direction_id,
                               'velocitiesTurbines%i' % direction_id])
        else:
            self.add('CtCp', AdjustCtCpYaw(nTurbines, direction_id, differentiable),
                     promotes=['Ct_in', 'Cp_in', 'params:*', 'floris_params:*', 'yaw%i' % direction_id])

        self.add('myFloris', FLORIS(nTurbines, resolution, direction_id, differentiable),
                 promotes=['floris_params:*', 'wind_speed', 'wind_direction', 'air_density', 'axialInduction',
                           'generator_efficiency', 'turbineX', 'turbineY', 'rotorDiameter', 'yaw%i' % direction_id,
                           'velocitiesTurbines%i' % direction_id, 'wt_power%i' % direction_id, 'power%i' % direction_id, 'wakeCentersYT', 'wakeDiametersT', 'wakeOverlapTRel'])
        if use_rotor_components:
            self.connect('CtCp.Ct_out', 'myFloris.Ct')
            self.connect('CtCp.Cp_out', 'myFloris.Cp')
        else:
            self.connect('floris_params:CTcorrected', 'params:CTcorrected')
            self.connect('floris_params:CPcorrected', 'params:CPcorrected')
            self.connect('CtCp.Ct_out', 'myFloris.Ct')
            self.connect('CtCp.Cp_out', 'myFloris.Cp')


class AEPGroupFLORIS(Group):
    """
    Group containing all necessary components for wind plant AEP calculations using the FLORIS model
    """

    def __init__(self, nTurbines, resolution=0, nDirections=1, use_rotor_components=False, datasize=0,
                 differentiable=True):

        super(AEPGroupFLORIS, self).__init__()

        # providing default unit types for general MUX/DeMUX components
        power_units = 'kW'
        direction_units = 'deg'
        wind_speed_units = 'm/s'

        # add necessary inputs for group
        self.add('p0', IndepVarComp('windDirections', np.zeros(nDirections), units=direction_units), promotes=['*'])
        self.add('p1', IndepVarComp('windSpeeds', np.zeros(nDirections), units=wind_speed_units), promotes=['*'])
        self.add('p2', IndepVarComp('windrose_frequencies', np.ones(nDirections)), promotes=['*'])
        self.add('p3', IndepVarComp('turbineX', np.zeros(nTurbines), units='m'), promotes=['*'])
        self.add('p4', IndepVarComp('turbineY', np.zeros(nTurbines), units='m'), promotes=['*'])

        # add vars to be seen by MPI and gradient calculations
        self.add('p5', IndepVarComp('rotorDiameter', np.zeros(nTurbines), units='m'), promotes=['*'])
        self.add('p6', IndepVarComp('axialInduction', np.zeros(nTurbines)), promotes=['*'])
        self.add('p7', IndepVarComp('generator_efficiency', np.zeros(nTurbines)), promotes=['*'])
        self.add('p8', IndepVarComp('air_density', val=1.1716, units='kg/(m*m*m)'), promotes=['*'])

        if not use_rotor_components:
            self.add('p9', IndepVarComp('Ct_in', np.zeros(nTurbines)), promotes=['*'])
            self.add('p10', IndepVarComp('Cp_in', np.zeros(nTurbines)), promotes=['*'])

        # add components and groups
        self.add('windDirectionsDeMUX', DeMUX(nDirections, units=direction_units))
        self.add('windSpeedsDeMUX', DeMUX(nDirections, units=wind_speed_units))

        pg = self.add('all_directions', ParallelGroup(), promotes=['*'])
        if use_rotor_components:
            for direction_id in range(0, nDirections):
                pg.add('direction_group%i' % direction_id,
                       DirectionGroupFLORIS(nTurbines=nTurbines, resolution=resolution, direction_id=direction_id,
                                            use_rotor_components=use_rotor_components, datasize=datasize,
                                            differentiable=differentiable),
                       promotes=['params:*', 'floris_params:*', 'air_density',
                                 'axialInduction', 'generator_efficiency', 'turbineX', 'turbineY',
                                 'yaw%i' % direction_id, 'rotorDiameter', 'velocitiesTurbines%i' % direction_id,
                                 'wt_power%i' % direction_id, 'power%i' % direction_id])#, 'wakeCentersYT', 'wakeDiametersT'])
        else:
            for direction_id in range(0, nDirections):
                pg.add('direction_group%i' % direction_id,
                       DirectionGroupFLORIS(nTurbines=nTurbines, resolution=resolution, direction_id=direction_id,
                                            use_rotor_components=use_rotor_components, datasize=datasize,
                                            differentiable=differentiable),
                       promotes=['Ct_in', 'Cp_in', 'params:*', 'floris_params:*', 'air_density', 'axialInduction',
                                 'generator_efficiency', 'turbineX', 'turbineY', 'yaw%i' % direction_id, 'rotorDiameter',
                                 'velocitiesTurbines%i' % direction_id, 'wt_power%i' % direction_id,
                                 'power%i' % direction_id])#, 'wakeCentersYT', 'wakeDiametersT'])

        self.add('powerMUX', MUX(nDirections, units=power_units))
        self.add('AEPcomp', WindFarmAEP(nDirections), promotes=['*'])

        # connect components
        self.connect('windDirections', 'windDirectionsDeMUX.Array')
        self.connect('windSpeeds', 'windSpeedsDeMUX.Array')
        for direction_id in range(0, nDirections):
            self.add('y%i' % direction_id, IndepVarComp('yaw%i' % direction_id, np.zeros(nTurbines), units='deg'), promotes=['*'])
            self.connect('windDirectionsDeMUX.output%i' % direction_id, 'direction_group%i.wind_direction' % direction_id)
            self.connect('windSpeedsDeMUX.output%i' % direction_id, 'direction_group%i.wind_speed' % direction_id)
            self.connect('power%i' % direction_id, 'powerMUX.input%i' % direction_id)
        self.connect('powerMUX.Array', 'power_directions')


# Testing code for development only
if __name__ == "__main__":
    # top = Problem()
    #
    # root = top.root = Group()
    #
    # root.add('p1', IndepVarComp('x', np.array([3.0])))
    # root.add('p2', IndepVarComp('y', np.array([2.0])))
    # root.add('p3', IndepVarComp('z', np.array([10.0])))
    # root.add('p', AdjustCtCpYaw(nTurbines=np.array([1])))
    #
    # root.connect('p1.x', 'p.Ct_in')
    # root.connect('p2.y', 'p.Cp_in')
    # root.connect('p3.z', 'p.yaw')
    #
    # top.setup()
    # top.check_partial_derivatives()
    # top.run()
    #
    # print(root.p.unknowns['Ct_out'])
    # print(root.p.unknowns['Cp_out'])
    #
    # top = Problem()
    #
    # root = top.root = Group()
    #
    # root.add('p1', IndepVarComp('x', np.array([10.0])))
    # root.add('p2', IndepVarComp('y', np.array([10.0])))
    # root.add('p3', IndepVarComp('z', 90.))
    # root.add('p', WindFrame(nTurbines=np.array([1]), resolution=0))
    #
    # root.connect('p1.x', 'p.turbineX')
    # root.connect('p2.y', 'p.turbineY')
    # root.connect('p3.z', 'p.wind_direction')
    #
    # top.setup()
    # top.check_partial_derivatives()
    # top.run()
    #
    # print(root.p.unknowns['turbineXw'])
    # print(root.p.unknowns['turbineYw'])
    #
    # top = Problem()
    #
    # root = top.root = Group()
    #
    # root.add('p1', IndepVarComp('x', np.array([10.0])))
    # root.add('p2', IndepVarComp('y', np.array([10.0])))
    # root.add('p', floris_wcent_wdiam(nTurbines=np.array([1])))
    #
    # root.connect('p1.x', 'p.turbineXw')
    # root.connect('p2.y', 'p.turbineYw')
    # root.connect('p1.x', 'p.yaw')
    # root.connect('p1.x', 'p.rotorDiameter')
    # root.connect('p1.x', 'p.Ct')
    #
    # top.setup()
    # top.check_partial_derivatives()
    # #top.run()
    #
    # print(root.p.unknowns['wakeDiametersT'])
    # print(root.p.unknowns['wakeCentersYT'])
    #
    # top = Problem()
    #
    # root = top.root = Group()
    #
    # root.add('p1', IndepVarComp('x', np.array([10.0])))
    # root.add('p2', IndepVarComp('y', np.array([10.0, 10.0, 10.0])))
    # root.add('p', floris_overlap(nTurbines=np.array([1])))
    #
    # root.connect('p1.x', 'p.turbineXw')
    # root.connect('p1.x', 'p.turbineYw')
    # root.connect('p1.x', 'p.rotorDiameter')
    # root.connect('p1.x', 'p.wakeCentersYT')
    # root.connect('p2.y', 'p.wakeDiametersT')
    #
    # top.setup()
    # top.check_partial_derivatives()
    # #top.run()
    #
    # print(root.p.unknowns['wakeOverlapTRel'])

    top = Problem()

    root = top.root = Group()
    wakeOL = np.array([ 0.94882764,  0.,          0.,          0.,          0.005853,    0.,          0.,
    0. ,         0.00603356,  0.,          0.,          0.,          0.,
    0.94876119,  0. ,         0. ,         0.,          0.00585258 , 0. ,         0.,
    0.  ,        0.00603362,  0. ,         0. ,         0. ,         0.,
    0.94882764 , 0.    ,      0.  ,        0. ,         0.005853   , 0. ,         0.,
    0.        ,  0.00603356 , 0.  ,        0.   ,       0. ,         0.,
    0.94837338,  0.  ,        0.   ,       0.  ,        0.00585014 , 0. ,         0.,
    0.    ,      0.00603391])

    root.add('p1', IndepVarComp('x', np.array([10.0, 10.0, 20, 20])))
    root.add('p2', IndepVarComp('y', wakeOL), promotes=['*'])
    root.add('p', floris_power(nTurbines=np.array([4])))

    root.connect('y', 'p.wakeOverlapTRel')

    top.setup()

    top.run()
    top.check_partial_derivatives()
    print(root.p.unknowns['power'])
    print(root.p.unknowns['wt_power'])
    print(root.p.unknowns['velocitiesTurbines'])
    d_wakeOL = np.zeros([wakeOL.size])
    step = 200
    for i in range(0, 2):
        top.run()
        shifthigh = np.zeros_like(wakeOL)
        shifthigh[i] = step
        shiftlow = np.zeros_like(wakeOL)
        shiftlow[i] = -step
        print shifthigh, shiftlow
        top['y'] = wakeOL+shifthigh
        print 'topy = ', top['y']
        top.setup()
        top.run()
        high = top.root.p.unknowns['power']
        top.root.p.params['wakeOverlapTRel'] = wakeOL+shiftlow
        print 'topy = ', top['y'], wakeOL+shiftlow
        top.setup()
        top.run()
        low = top.root.p.unknowns['power']
        print high, low
        d_wakeOL[i] = (high-low)/(2*step)
    print 'd_wakeOL: ', d_wakeOL



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
