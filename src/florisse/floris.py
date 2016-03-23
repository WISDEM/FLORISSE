import numpy as np

from openmdao.api import Group, Component, Problem, IndepVarComp, ParamComp, ParallelGroup
from openmdao.api import NLGaussSeidel, ScipyGMRES

from GeneralWindFarmComponents import WindFrame, AdjustCtCpYaw, MUX, WindFarmAEP, DeMUX, \
    CPCT_Interpolate_Gradients_Smooth, WindDirectionPower, add_gen_params_IdepVarComps, \
    CPCT_Interpolate_Gradients

import _floris
import _florisDiscontinuous
# import _florisHubSmooth as _floris


def add_floris_parameters(openmdao_comp, use_rotor_components=True):
    # altering the values in this function will have no effect during optimization. To change defaults permanently,
    # alter the values in add_floris_IndepVarComps().

    ###################   wake deflection   ##################

    ### parameters
    # original model
    openmdao_comp.add_param('floris_params:kd', 0.15 if not use_rotor_components else 0.17,
                            desc='model parameter that defines the sensitivity of the wake deflection to yaw')
    openmdao_comp.add_param('floris_params:initialWakeDisplacement', -4.5, pass_by_obj=True,
                            desc='defines the wake at the rotor to be slightly offset from the rotor. This is'
                                 'necessary for tuning purposes')
    openmdao_comp.add_param('floris_params:bd', -0.01, pass_by_obj=True,
                            desc='defines rate of wake displacement if initialWakeAngle is not used')
    # added
    openmdao_comp.add_param('floris_params:initialWakeAngle', 0.5*3.0, pass_by_obj=True,
                            desc='sets how angled the wake flow should be at the rotor')

    ### flags
    openmdao_comp.add_param('floris_params:useWakeAngle', False if not use_rotor_components else True, pass_by_obj=True,
                            desc='define whether an initial angle or initial offset should be used for wake center. '
                                 'if True, then bd will be ignored and initialWakeAngle will'
                                 'be used. The reverse is also true')


    ###################   wake expansion   ##################

    ### parameters
    # original model
    openmdao_comp.add_param('floris_params:ke', 0.065 if not use_rotor_components else 0.05, pass_by_obj=True,
                            desc='parameter defining overall wake expansion')
    openmdao_comp.add_param('floris_params:me', np.array([-0.5, 0.22, 1.0]) if not use_rotor_components else np.array([-0.5, 0.3, 1.0]),
                            pass_by_obj=True,
                            desc='parameters defining relative zone expansion. Mixing zone (me[2]) must always be 1.0')

    ### flags
    openmdao_comp.add_param('floris_params:adjustInitialWakeDiamToYaw', False if not use_rotor_components else True,
                            pass_by_obj=True,
                            desc='if True then initial wake diameter will be set to rotorDiameter*cos(yaw)')


    ###################   wake velocity   ##################

    ### parameters
    # original model
    openmdao_comp.add_param('floris_params:MU', np.array([0.5, 1.0, 5.5]), pass_by_obj=True,
                            desc='velocity deficit decay rates for each zone. Middle zone must always be 1.0')
    openmdao_comp.add_param('floris_params:aU', 5.0 if not use_rotor_components else 12.0, units='deg', pass_by_obj=True,
                            desc='zone decay adjustment parameter independent of yaw')
    openmdao_comp.add_param('floris_params:bU', 1.66 if not use_rotor_components else 1.3, pass_by_obj=True,
                            desc='zone decay adjustment parameter dependent yaw')
    # added
    openmdao_comp.add_param('floris_params:cos_spread', 1E12 if not use_rotor_components else 2.0, pass_by_obj=True,
                            desc='spread of cosine smoothing factor (multiple of sum of wake and rotor radii)')
    openmdao_comp.add_param('floris_params:keCorrArray', 0.0, pass_by_obj=True,
                            desc='multiplies the ke value by 1+keCorrArray*(sum of rotors relative overlap with '
                                 'inner two zones for including array affects')
    openmdao_comp.add_param('floris_params:keCorrCT', 0.0, pass_by_obj=True,
                            desc='adjust ke by adding a precentage of the difference of CT and ideal CT as defined in'
                                 'Region2CT')
    openmdao_comp.add_param('floris_params:Region2CT', 4.0*(1.0/3.0)*(1.0-(1.0/3.0)), pass_by_obj=True,
                            desc='defines ideal CT value for use in adjusting ke to yaw adjust CT if keCorrCT>0.0')

    # flags
    openmdao_comp.add_param('floris_params:axialIndProvided', True if not use_rotor_components else False, pass_by_obj=True,
                            desc='if axial induction is not provided, then it will be calculated based on CT')
    openmdao_comp.add_param('floris_params:useaUbU', True, pass_by_obj=True,
                            desc='if True then zone velocity decay rates (MU) will be adjusted based on yaw')


    ###################   other   ##################
    openmdao_comp.add_param('floris_params:FLORISoriginal', False, pass_by_obj=True,
                            desc='override all parameters and use FLORIS as original in first Wind Energy paper')


    ####### apply to things now external to the Floris model
    # openmdao_comp.add_param('floris_params:CTcorrected', True, pass_by_obj=True,
    #                desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
    #
    # openmdao_comp.add_param('floris_params:CPcorrected', True, pass_by_obj=True,
    #                desc = 'CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)')
    # openmdao_comp.add_param('floris_params:pP', 1.88, pass_by_obj=True)

    ####### unused
    # openmdao_comp.add_param('floris_params:baselineCT', 4./3.*(1.-1./3.), pass_by_obj=True)
    # openmdao_comp.add_param('floris_params:keCorrTI', 0.0, pass_by_obj=True)
    # openmdao_comp.add_param('floris_params:baselineTI', 0.045, pass_by_obj=True)
    # openmdao_comp.add_param('floris_params:keCorrHR', 0.0, pass_by_obj=True) # neutral, with heating rate 0, is baseline
    # openmdao_comp.add_param('floris_params:keCorrHRTI', 0.0, pass_by_obj=True)
    # openmdao_comp.add_param('floris_params:keSaturation', 0.0, pass_by_obj=True)
    # openmdao_comp.add_param('floris_params:kdCorrYawDirection', 0.0, pass_by_obj=True)


def add_floris_params_IndepVarComps(openmdao_object, use_rotor_components=True):
    # permanently alter defaults here
    ###################   wake deflection   ##################

    ### parameters
    # original model
    openmdao_object.add('fp00', IndepVarComp('floris_params:kd', 0.15 if not use_rotor_components else 0.17,
                                             desc='model parameter that defines the sensitivity of the wake deflection '
                                                  'to yaw'),
                        promotes=['*'])
    openmdao_object.add('fp01', IndepVarComp('floris_params:initialWakeDisplacement', -4.5, pass_by_obj=True,
                                             desc='defines the wake at the rotor to be slightly offset from the rotor. '
                                                  'This is necessary for tuning purposes'),
                        promotes=['*'])
    openmdao_object.add('fp02', IndepVarComp('floris_params:bd', -0.01, pass_by_obj=True,
                                             desc='defines rate of wake displacement if initialWakeAngle is not used'),
                        promotes=['*'])
    # added
    openmdao_object.add('fp03', IndepVarComp('floris_params:initialWakeAngle', 1.5, pass_by_obj=True,
                                             desc='sets how angled the wake flow should be at the rotor'),
                        promotes=['*'])

    ### flags
    openmdao_object.add('fp04', IndepVarComp('floris_params:useWakeAngle', False if not use_rotor_components else True,
                                             pass_by_obj=True,
                                             desc='define whether an initial angle or initial offset should be used for'
                                                  'wake center. If True, then bd will be ignored and initialWakeAngle '
                                                  'will be used. The reverse is also true'),
                        promotes=['*'])


    ###################   wake expansion   ##################

    ### parameters
    # original model
    openmdao_object.add('fp05', IndepVarComp('floris_params:ke', 0.065 if not use_rotor_components else 0.05, pass_by_obj=True,
                                             desc='parameter defining overall wake expansion'),
                        promotes=['*'])
    openmdao_object.add('fp06', IndepVarComp('floris_params:me', np.array([-0.5, 0.22, 1.0]) if not use_rotor_components else np.array([-0.5, 0.3, 1.0]),
                                             pass_by_obj=True,
                                             desc='parameters defining relative zone expansion. Mixing zone (me[2]) '
                                                  'must always be 1.0'),
                        promotes=['*'])

    ### flags
    openmdao_object.add('fp07', IndepVarComp('floris_params:adjustInitialWakeDiamToYaw',
                                             False if not use_rotor_components else True, pass_by_obj=True,
                                             desc='if True then initial wake diameter will be set to '
                                                  'rotorDiameter*cos(yaw)'),
                        promotes=['*'])


    ###################   wake velocity   ##################

    ### parameters
    # original model
    openmdao_object.add('fp08', IndepVarComp('floris_params:MU', np.array([0.5, 1.0, 5.5]), pass_by_obj=True,
                                             desc='velocity deficit decay rates for each zone. Middle zone must always '
                                                  'be 1.0'),
                        promotes=['*'])
    openmdao_object.add('fp09', IndepVarComp('floris_params:aU', 5.0 if not use_rotor_components else 12.0, units='deg', pass_by_obj=True,
                                             desc='zone decay adjustment parameter independent of yaw'),
                        promotes=['*'])
    openmdao_object.add('fp10', IndepVarComp('floris_params:bU', 1.66 if not use_rotor_components else 1.3, pass_by_obj=True,
                                             desc='zone decay adjustment parameter dependent yaw'),
                        promotes=['*'])
    # added
    openmdao_object.add('fp11', IndepVarComp('floris_params:cos_spread', 2.0, pass_by_obj=True,
                                             desc='spread of cosine smoothing factor (multiple of sum of wake and rotor '
                                                  'radii)'),
                        promotes=['*'])
    openmdao_object.add('fp12', IndepVarComp('floris_params:keCorrArray', 0.0, pass_by_obj=True,
                                             desc='multiplies the ke value by 1+keCorrArray*(sum of rotors relative '
                                                  'overlap with inner two zones for including array affects'),
                        promotes=['*'])
    openmdao_object.add('fp13', IndepVarComp('floris_params:keCorrCT', 0.0, pass_by_obj=True,
                                             desc='adjust ke by adding a precentage of the difference of CT and ideal '
                                                  'CT as defined in Region2CT'),
                        promotes=['*'])
    openmdao_object.add('fp14', IndepVarComp('floris_params:Region2CT', 4.0*(1.0/3.0)*(1.0-(1.0/3.0)), pass_by_obj=True,
                                             desc='defines ideal CT value for use in adjusting ke to yaw adjust CT if '
                                                  'keCorrCT>0.0'),
                        promotes=['*'])

    # flags
    openmdao_object.add('fp15', IndepVarComp('floris_params:axialIndProvided', True if not use_rotor_components else False,
                                             pass_by_obj=True,
                                             desc='if axial induction is not provided, then it will be calculated based '
                                                  'on CT'),
                        promotes=['*'])
    openmdao_object.add('fp16', IndepVarComp('floris_params:useaUbU', True, pass_by_obj=True,
                                             desc='if True then zone velocity decay rates (MU) will be adjusted based '
                                                  'on yaw'),
                        promotes=['*'])


    ###################   other   ##################
    openmdao_object.add('fp17', IndepVarComp('floris_params:FLORISoriginal', False, pass_by_obj=True,
                                             desc='override all parameters and use FLORIS as original in Gebraad et al.'
                                                  '2014, Wind plant power optimization through yaw control using a '
                                                  'parametric model for wake effect-a CFD simulation study'),
                        promotes=['*'])


# Components of FLORIS - for full model use FLORIS(Group)
class floris_wcent_wdiam(Component):

    def __init__(self, nTurbines, direction_id=0, differentiable=True, splineshift=0.0, use_rotor_components=True):

        super(floris_wcent_wdiam, self).__init__()
        
        self.direction_id = direction_id
        self.differentiable = differentiable
        self.splineshift = splineshift

        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-3
        self.fd_options['step_type'] = 'relative'

        if not differentiable:
            self.fd_options['force_fd'] = True
            self.fd_options['form'] = 'forward'

        # print 'entering wcent_wdiam __init__ - Tapenade'

        # input arrays
        self.add_param('turbineXw', np.zeros(nTurbines), units='m', desc='x coordinates of turbines in wind dir. ref. frame')
        self.add_param('turbineYw', np.zeros(nTurbines), units='m', desc='y coordinates of turbines in wind dir. ref. frame')
        self.add_param('yaw%i' % direction_id, np.zeros(nTurbines), units='deg', desc='yaw of each turbine')
        self.add_param('rotorDiameter', np.zeros(nTurbines) + 126.4, units='m', desc='rotor diameter of each turbine')
        self.add_param('Ct', np.zeros(nTurbines)+4.0*(1./3.)*(1.0-(1./3.)), desc='thrust coefficient of each turbine')

        # output arrays
        self.add_output('wakeCentersYT', np.zeros(nTurbines*nTurbines), units='m', desc='wake center y position at each turbine')
        self.add_output('wakeDiametersT', np.zeros(3*nTurbines*nTurbines), units='m', desc='wake diameter of each zone of each wake at each turbine')

        # FLORIS parameters
        add_floris_parameters(self, use_rotor_components=use_rotor_components)

    def solve_nonlinear(self, params, unknowns, resids):

        # print 'entering wcent_wdiam - tapenade'

        rotorDiameter = params['rotorDiameter']
        Ct = params['Ct']

        Region2CT = params['floris_params:Region2CT']

        # if params['floris_params:FLORISoriginal']:
        #     ke = 0.065
        #     keCorrCT = 0.0
        #     kd = 0.15
        #     me = np.array([-0.5, 0.22, 1.0])
        #     useWakeAngle = False
        #     initialWakeDisplacement = -4.5
        #     initialWakeAngle = params['floris_params:initialWakeAngle']
        #     bd = -0.01
        #     adjustInitialWakeDiamToYaw = False
        #
        # else:
        ke = params['floris_params:ke']
        keCorrCT = params['floris_params:keCorrCT']
        kd = params['floris_params:kd']
        me = params['floris_params:me']
        initialWakeDisplacement = params['floris_params:initialWakeDisplacement']
        useWakeAngle = params['floris_params:useWakeAngle']
        initialWakeAngle = params['floris_params:initialWakeAngle']
        bd = params['floris_params:bd']
        adjustInitialWakeDiamToYaw = params['floris_params:adjustInitialWakeDiamToYaw']
        # print "ke = %f, kd = %f" % (ke, kd)

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
            # print "not differentiable"
            # print "ke = ", ke
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

        # if params['floris_params:FLORISoriginal']:
        #
        #     ke = 0.065
        #     keCorrCT = 0.0
        #     kd = 0.15
        #     me = np.array([-0.5, 0.22, 1.0])
        #     useWakeAngle = False
        #     initialWakeDisplacement = -4.5
        #     initialWakeAngle = params['floris_params:initialWakeAngle']
        #     bd = -0.01
        #     adjustInitialWakeDiamToYaw = False
        #
        # else:
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
        # yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb = \
        #     _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement, initialWakeAngle, ke, keCorrCT, Region2CT,
        #                                   yaw_deg, Ct, turbineXw, turbineYw, rotorDiameter, me, bd, self.splineshift, useWakeAngle,
        #                                   adjustInitialWakeDiamToYaw, wakeCentersYT_vecb, wakeDiametersT_vecb)
        #
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
        # yawb, Ctb, turbineXwb, turbineYwb, rotorDiameterb = \
        #     _floris.floris_wcent_wdiam_bv(kd, initialWakeDisplacement, initialWakeAngle, ke, keCorrCT, Region2CT,
        #                                   yaw_deg, Ct, turbineXw, turbineYw, rotorDiameter, me, bd, self.splineshift, useWakeAngle,
        #                                   adjustInitialWakeDiamToYaw, wakeCentersYT_vecb, wakeDiametersT_vecb)
        #     #
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

    def __init__(self, nTurbines, differentiable=True, splineshift=0.0, use_rotor_components=True):

        super(floris_overlap, self).__init__()

        self.differentiable = differentiable
        self.splineshift = splineshift

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
        self.add_param('rotorDiameter', np.zeros(nTurbines)+126.4, units='m',
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

        # FLORIS parameters
        add_floris_parameters(self, use_rotor_components=use_rotor_components)

        # etc
        self.nTurbines = nTurbines

    def solve_nonlinear(self, params, unknowns, resids):

        # print 'entering overlap - Tapenade'

        # call to fortran code to obtain relative wake overlap values
        # print params['turbineXw'], params['turbineYw'], params['rotorDiameter']
        if self.differentiable:
            wakeOverlapTRel_vec, cosFac_vec = _floris.floris_overlap(params['turbineXw'], params['turbineYw'],
                                                                     params['rotorDiameter'], params['wakeDiametersT'],
                                                                     params['wakeCentersYT'],
                                                                     params['floris_params:cos_spread'])
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
        # turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb\
        #     = _floris.floris_overlap_bv(params['turbineXw'], params['turbineYw'], params['rotorDiameter'],
        #                                 params['wakeDiametersT'], params['wakeCentersYT'],
        #                                 params['floris_params:cos_spread'], self.splineshift, wakeOverlapTRel_vecb, cosFac_vecb)

        turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb \
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
        # turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb\
        #     = _floris.floris_overlap_bv(params['turbineXw'], params['turbineYw'], params['rotorDiameter'],
        #                                 params['wakeDiametersT'], params['wakeCentersYT'],
        #                                 params['floris_params:cos_spread'], self.splineshift, wakeOverlapTRel_vecb, cosFac_vecb)
        #
        turbineYwb, rotorDiameterb, wakeDiametersT_vecb, wakeCentersYT_vecb \
            = _floris.floris_overlap_bv(params['turbineXw'], params['turbineYw'], params['rotorDiameter'],
                                        params['wakeDiametersT'], params['wakeCentersYT'],
                                        params['floris_params:cos_spread'], wakeOverlapTRel_vecb, cosFac_vecb)

        # construct Jacobian of floris_overlap
        J['cosFac', 'turbineYw'] = turbineYwb
        J['cosFac', 'rotorDiameter'] = rotorDiameterb
        J['cosFac', 'wakeDiametersT'] = wakeDiametersT_vecb
        J['cosFac', 'wakeCentersYT'] = wakeCentersYT_vecb

        return J


class floris_velocity(Component):
    """ Calculates the turbine power and effective wind speed for each turbine """

    def __init__(self, nTurbines, direction_id=0, differentiable=True, splineshift=0.0, use_rotor_components=True):

        super(floris_velocity, self).__init__()

        self.differentiable = differentiable
        self.splineshift = splineshift
        self.nTurbines = nTurbines
        self.direction_id = direction_id

        self.fd_options['form'] = 'central'
        self.fd_options['step_size'] = 1.0e-6
        self.fd_options['step_type'] = 'relative'

        if not differentiable:
            self.fd_options['force_fd'] = True
            self.fd_options['form'] = 'forward'

        # print 'entering power __init__ - Tapenade'


        # inputs
        self.add_param('wind_speed', 8.0, units='m/s', desc='free stream wind velocity')
        self.add_param('air_density', 1.1716, units='kg/(m*m*m)', desc='air density in free stream')
        self.add_param('rotorDiameter', np.zeros(nTurbines)+126.4, units='m', desc='rotor diameters of all turbine')
        self.add_param('axialInduction', np.zeros(nTurbines)+1./3., desc='axial induction of all turbines')
        self.add_param('Ct', np.zeros(nTurbines)+4.0*(1./3.)*(1.0-(1./3.)), desc='Thrust coefficient for all turbines')
        # self.add_param('Cp', np.zeros(nTurbines)+0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2), desc='power coefficient for all turbines')
        # self.add_param('generator_efficiency', np.zeros(nTurbines)+0.944, desc='generator efficiency of all turbines')
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
        self.add_output('velocitiesTurbines%i' % direction_id, val=np.zeros(nTurbines), units='m/s',
                       desc='effective hub velocity for each turbine')

        # FLORIS parameters
        add_floris_parameters(self, use_rotor_components=use_rotor_components)

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
        yaw = params['yaw%i' % self.direction_id]

        # print 'wake OL in power', wakeOverlapTRel_v

        # set floris parameter values
        # if params['floris_params:FLORISoriginal']:
        #
        #     ke = 0.065
        #     keCorrCT = 0.0
        #     keCorrArray = 0.0
        #     MU = np.array([0.5, 1.0, 5.5])
        #     useaUbU = True
        #     aU = 5.0
        #     bU = 1.66
        #     axialIndProvided = True
        #
        # else:
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
        # print self.differentiable
        if self.differentiable:
            # TODO hard code this parameter in the fortran code and remove the specifier from all functions of this component
            p_near0 = -1.0
            # p_near0 = 0.0
            # splineshift = self.splineshift

            # pass p_near0 to self for use in gradient calculations
            self.p_near0 = p_near0

            # call to fortran code to obtain output values
            velocitiesTurbines = _floris.floris_velocity(wakeOverlapTRel_v, cosFac_v, Ct, axialInduction,
                                                                       axialIndProvided, useaUbU, keCorrCT, Region2CT,
                                                                       ke, Vinf, keCorrArray, turbineXw, yaw,
                                                                       rotorDiameter, MU, aU, bU)
        else:
            velocitiesTurbines = _florisDiscontinuous.floris_velocity(wakeOverlapTRel_v, Ct, axialInduction, \
                                                                axialIndProvided, useaUbU, keCorrCT, Region2CT, ke, \
                                                                Vinf, keCorrArray, turbineXw, yaw, rotorDiameter, MU, \
                                                                aU, bU)

        # pass outputs to self
        unknowns['velocitiesTurbines%i' % self.direction_id] = velocitiesTurbines
        # unknowns['wt_power%i' % self.direction_id] = wt_power
        # unknowns['power%i' % self.direction_id] = power

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
        # rho = params['air_density']
        # generator_efficiency = params['generator_efficiency']
        # Cp = params['Cp']
        yaw = params['yaw%i' % self.direction_id]

        # print 'wake OL in power', wakeOverlapTRel_v

        # set floris parameter values
        # if params['floris_params:FLORISoriginal']:
        #
        #     ke = 0.065
        #     keCorrCT = 0.0
        #     keCorrArray = 0.0
        #     MU = np.array([0.5, 1.0, 5.5])
        #     useaUbU = True
        #     aU = 5.0
        #     bU = 1.66
        #     axialIndProvided = True
        #
        # else:

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
        # p_near0 = 1.0

        # pass p_near0 to self for use in gradient calculations
        # self.p_near0 = p_near0
        # print self.differentiable
        if self.differentiable:
            # TODO hard code this parameter in the fortran code and remove the specifier from all functions of this component
            # p_near0 = self.p_near0
            # splineshift = self.splineshift

            # pass p_near0 to self for use in gradient calculations
            # self.p_near0 = p_near0

            # call to fortran code to obtain output values
            velocitiesTurbines = _floris.floris_velocity(wakeOverlapTRel_v, cosFac_v, Ct, axialInduction,
                                                                       axialIndProvided, useaUbU, keCorrCT, Region2CT,
                                                                       ke, Vinf, keCorrArray, turbineXw, yaw,
                                                                       rotorDiameter, MU, aU, bU)
        else:
            velocitiesTurbines = _florisDiscontinuous.floris_velocity(wakeOverlapTRel_v, Ct, axialInduction, \
                                                                axialIndProvided, useaUbU, keCorrCT, Region2CT, ke, \
                                                                Vinf, keCorrArray, turbineXw, yaw, rotorDiameter, MU, \
                                                                aU, bU)
            # print 'here'

        # pass outputs to self
        resids['velocitiesTurbines%i' % self.direction_id] = velocitiesTurbines - unknowns['velocitiesTurbines%i' %
                                                                                           self.direction_id]
        # resids['wt_power%i' % self.direction_id] = wt_power - unknowns['wt_power%i' % self.direction_id]
        # resids['power%i' % self.direction_id] = power - unknowns['power%i' % self.direction_id]
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
        # rho = params['air_density']
        # generator_efficiency = params['generator_efficiency']
        # Cp = params['Cp']
        yaw = params['yaw%i' % self.direction_id]

        # set floris parameter values
        # if params['floris_params:FLORISoriginal']:
        #
        #     ke = 0.065
        #     keCorrCT = 0.0
        #     keCorrArray = 0.0
        #     MU = np.array([0.5, 1.0, 5.5])
        #     useaUbU = True
        #     aU = 5.0
        #     bU = 1.66
        #     axialIndProvided = True
        #
        # else:
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
        # p_near0 = -1.0

        # create jacobian dict
        J = {}

        # input arrays to direct differentiation
        velocitiesTurbinesb = np.eye(nbdirs, nTurbines)
        # velocitiesTurbinesb = np.zeros((nbdirs, nTurbines))
        # velocitiesTurbinesb[:, 0] = 1.0
        # wt_powerb = np.zeros((nbdirs, nTurbines))
        # powerb = np.zeros(nbdirs)
        # print np.size(velocitiesTurbinesb), velocitiesTurbinesb
        # call to fortran to obtain gradients of velocitiesTurbines
        # wakeOverlapTRel_vb, cosFac_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb\
        #     = _floris.floris_power_bv(wakeOverlapTRel_v, cosFac, Ct, axialInduction, axialIndProvided, useaUbU,
        #                               keCorrCT, Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, self.splineshift,
        #                               rotorDiameter, MU, rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb,
        #                               wt_powerb, powerb)

          # call to fortran to obtain gradients of velocitiesTurbines
        wakeOverlapTRel_vb, cosFac_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb\
            = _floris.floris_velocity_bv(wakeOverlapTRel_v, cosFac, Ct, axialInduction, axialIndProvided, useaUbU,
                                      keCorrCT, Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw,
                                      rotorDiameter, MU, aU, bU, velocitiesTurbinesb)
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
        # J['velocitiesTurbines%i' % direction_id, 'Cp'] = Cpb
        J['velocitiesTurbines%i' % direction_id, 'axialInduction'] = axialInductionb
        J['velocitiesTurbines%i' % direction_id, 'turbineXw'] = turbineXwb
        J['velocitiesTurbines%i' % direction_id, 'yaw%i' % direction_id] = yawb
        J['velocitiesTurbines%i' % direction_id, 'rotorDiameter'] = rotorDiameterb

        # input arrays to direct differentiation
        # velocitiesTurbinesb[:, :] = 0.0
        # wt_powerb = np.eye(nbdirs, nTurbines)
        #
        # # call to fortran to obtain gradients wt_power
        # wakeOverlapTRel_vb, cosFac_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb\
        #     = _floris.floris_power_bv(wakeOverlapTRel_v, cosFac, Ct, axialInduction, axialIndProvided, useaUbU,
        #                               keCorrCT, Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, self.splineshift,
        #                               rotorDiameter, MU, rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb,
        #                               wt_powerb, powerb)

        # # print wakeOverlapTRel_vb.shape, wt_powerb.shape
        # # collect values of the jacobian of wt_power
        # J['wt_power%i' % direction_id, 'wakeOverlapTRel'] = wakeOverlapTRel_vb
        # J['wt_power%i' % direction_id, 'cosFac'] = cosFac_vb
        # J['wt_power%i' % direction_id, 'Ct'] = Ctb
        # J['wt_power%i' % direction_id, 'Cp'] = Cpb
        # J['wt_power%i' % direction_id, 'axialInduction'] = axialInductionb
        # J['wt_power%i' % direction_id, 'turbineXw'] = turbineXwb
        # J['wt_power%i' % direction_id, 'yaw%i' % direction_id] = yawb
        # J['wt_power%i' % direction_id, 'rotorDiameter'] = rotorDiameterb
        #
        # # input arrays to direct differentiation
        # wt_powerb[:, :] = 0.0
        # powerb[0] = 1.0
        #
        # # call to fortran to obtain gradients of power
        # wakeOverlapTRel_vb, cosFac_vb, Ctb, axialInductionb, turbineXwb, yawb, rotorDiameterb, Cpb\
        #     = _floris.floris_power_bv(wakeOverlapTRel_v, cosFac, Ct, axialInduction, axialIndProvided, useaUbU,
        #                               keCorrCT, Region2CT, ke, Vinf, keCorrArray, turbineXw, yaw, p_near0, self.splineshift,
        #                               rotorDiameter, MU, rho, aU, bU, Cp, generator_efficiency, velocitiesTurbinesb,
        #                               wt_powerb, powerb)
        #
        # #print 'OL', wakeOverlapTRel_vb[0, :]#, Ctb[0, :].shape, Cpb[0, :].shape, axialInductionb[0, :], turbineXwb[0, :].shape, yawb[0, :].shape, rotorDiameterb[0, :].shape
        #
        # # # print np.array(yawb[:1, :])
        # # collect values of the jacobian of wt_power
        # # print 'vb', wakeOverlapTRel_vb
        # J['power%i' % direction_id, 'wakeOverlapTRel'] = np.array(wakeOverlapTRel_vb[:1, :])
        # J['power%i' % direction_id, 'cosFac'] = np.array(cosFac_vb[:1, :])
        # J['power%i' % direction_id, 'Ct'] = np.array(Ctb[:1, :])
        # J['power%i' % direction_id, 'Cp'] = np.array(Cpb[:1, :])
        # J['power%i' % direction_id, 'axialInduction'] = np.array(axialInductionb[:1, :])
        # J['power%i' % direction_id, 'turbineXw'] = np.array(turbineXwb[:1, :])
        # J['power%i' % direction_id, 'yaw%i' % direction_id] = np.array(yawb[:1, :])
        # J['power%i' % direction_id, 'rotorDiameter'] = np.array(rotorDiameterb[:1, :])

        # print 'leaving power linearize'
        # print J
        return J


# Groups using FLORIS
class FLORIS(Group):
    """ Group containing all necessary components of the floris model """

    def __init__(self, nTurbines, resolution, direction_id=0, differentiable=True, optimizingLayout=False,
                 use_rotor_components=True):
        super(FLORIS, self).__init__()
        splineshift = 0.0

        # print differentiable
        # print 'in myFloris direction %i' % direction_id
        if optimizingLayout:
            splineshift = 1.0
        self.add('f_1', WindFrame(nTurbines, resolution, differentiable=differentiable), promotes=['*'])
        self.add('f_2', floris_wcent_wdiam(nTurbines, direction_id=direction_id, differentiable=differentiable,
                                           splineshift=splineshift, use_rotor_components=use_rotor_components),
                 promotes=['*'])
        self.add('f_3', floris_overlap(nTurbines, differentiable=differentiable, splineshift=splineshift,
                                       use_rotor_components=use_rotor_components),
                 promotes=['*'])
        self.add('f_4', floris_velocity(nTurbines, direction_id=direction_id, differentiable=differentiable,
                                     splineshift=splineshift, use_rotor_components=use_rotor_components),
                 promotes=['*'])


class RotorSolveGroup(Group):

    def __init__(self, nTurbines, resolution=0, direction_id=0, datasize=0,
                 differentiable=True, optimizingLayout=False, use_rotor_components=True):

        super(RotorSolveGroup, self).__init__()
        epsilon = 1E-12

        # self.nl_solver = NLGaussSeidel()
        # self.ln_solver = Brent()
        # self.nl_solver.options['atol'] = epsilon
        # self.nl_solver.options['rtol'] = epsilon
        # self.ln_solver = LinearGaussSeidel()
        self.ln_solver = ScipyGMRES()
        # self.ln_solver = DirectSolver()
        # self.ln_solver = DirectSolver()

        # self.ln_solver = PetscKSP()
        self.nl_solver = NLGaussSeidel()
        # self.nl_solver = Newton()
        # self.nl_solver.options['iprint'] = 1
        # self.ln_solver.options['iprint'] = 1
        # self.nl_solver.options['maxiter'] = 200
        # self.ln_solver = PetscKSP()
        # self.ln_solver = Backtracking()
        # self.ln_solver = Brent()
        # self.ln_solver.options['state_var'] = 'velocitiesTurbines%i' % direction_id
        # self.nl_solver.options['lower_bound'] = 0.
        # self.nl_solver.options['upper_bound'] = 100.
        self.ln_solver.options['atol'] = epsilon
        # self.ln_solver.options['rtol'] = epsilon
        # self.ln_solver.setup('direction_group%i' % direction_id)
        # self.ln_solver.solve('velocitiesTurbines%i' % direction_id)

        self.add('CtCp', CPCT_Interpolate_Gradients_Smooth(nTurbines, direction_id=direction_id, datasize=datasize),
                     promotes=['gen_params:*', 'yaw%i' % direction_id,
                               'velocitiesTurbines%i' % direction_id, 'Cp_out'])

        # self.add('CtCp', CPCT_Interpolate_Gradients(nTurbines, direction_id=direction_id, datasize=datasize),
        #              promotes=['gen_params:*', 'yaw%i' % direction_id,
        #                        'velocitiesTurbines%i' % direction_id, 'Cp_out'])

        self.add('floris', FLORIS(nTurbines, resolution=resolution, direction_id=direction_id,
                                  differentiable=differentiable, optimizingLayout=optimizingLayout,
                                  use_rotor_components=use_rotor_components),
                 promotes=['floris_params:*', 'wind_speed', 'wind_direction', 'air_density', 'axialInduction',
                           'turbineX', 'turbineY', 'rotorDiameter', 'yaw%i' % direction_id,
                           'velocitiesTurbines%i' % direction_id, 'wakeCentersYT', 'wakeDiametersT',
                           'wakeOverlapTRel'])
        self.connect('CtCp.Ct_out', 'floris.Ct')
        # self.connect('CtCp.Cp_out', 'Cp')


class DirectionGroupFLORIS(Group):
    """
    Group containing all necessary components for wind plant calculations
    in a single direction
    """

    def __init__(self, nTurbines, resolution=0, direction_id=0, use_rotor_components=True, datasize=0,
                 differentiable=True, optimizingLayout=False, add_IdepVarComps=True, forcefd=False):
        super(DirectionGroupFLORIS, self).__init__()
        epsilon = 1e-6

        self.fd_options['force_fd'] = forcefd

        # self.add('p0', IndepVarComp('wind_direction', val=0.0, units='deg'), promotes=['*'])
        if add_IdepVarComps:
            add_floris_params_IndepVarComps(self, use_rotor_components=use_rotor_components)
            add_gen_params_IdepVarComps(self, datasize=datasize)

        # self.add('fp', FLORISParameters(), promotes=['*'])
        if use_rotor_components:
            self.add('myFloris', RotorSolveGroup(nTurbines, resolution=resolution, direction_id=direction_id,
                                                 datasize=datasize, differentiable=differentiable,
                                                 optimizingLayout=optimizingLayout,
                                                 use_rotor_components=use_rotor_components),
                     promotes=['gen_params:*', 'yaw%i' % direction_id, 'velocitiesTurbines%i' % direction_id,
                               'floris_params:*', 'wind_speed', 'wind_direction', 'air_density', 'axialInduction',
                               'turbineX', 'turbineY', 'rotorDiameter', 'wakeCentersYT', 'wakeDiametersT',
                               'wakeOverlapTRel'])
        else:
            self.add('CtCp', AdjustCtCpYaw(nTurbines, direction_id, differentiable),
                     promotes=['Ct_in', 'Cp_in', 'gen_params:*', 'yaw%i' % direction_id])

            self.add('myFloris', FLORIS(nTurbines, resolution=resolution, direction_id=direction_id,
                                        differentiable=differentiable, optimizingLayout=optimizingLayout,
                                        use_rotor_components=use_rotor_components),
                     promotes=['floris_params:*', 'wind_speed', 'wind_direction', 'air_density', 'axialInduction',
                               'turbineX', 'turbineY', 'rotorDiameter', 'yaw%i' % direction_id,
                               'velocitiesTurbines%i' % direction_id, 'wakeCentersYT', 'wakeDiametersT',
                               'wakeOverlapTRel'])

        self.add('powerComp', WindDirectionPower(nTurbines=nTurbines, direction_id=direction_id, differentiable=True,
                                                 use_rotor_components=use_rotor_components),
                 promotes=['air_density', 'generator_efficiency', 'rotorDiameter',
                           'velocitiesTurbines%i' % direction_id,
                           'wt_power%i' % direction_id, 'power%i' % direction_id])

        if use_rotor_components:
            self.connect('myFloris.Cp_out', 'powerComp.Cp')
        else:
            self.connect('CtCp.Ct_out', 'myFloris.Ct')
            self.connect('CtCp.Cp_out', 'powerComp.Cp')


class AEPGroupFLORIS(Group):
    """
    Group containing all necessary components for wind plant AEP calculations using the FLORIS model
    """

    def __init__(self, nTurbines, resolution=0, nDirections=1, use_rotor_components=True, datasize=0,
                 differentiable=True, optimizingLayout=False, forcefd=False):

        super(AEPGroupFLORIS, self).__init__()

        # providing default unit types for general MUX/DeMUX components
        power_units = 'kW'
        direction_units = 'deg'
        wind_speed_units = 'm/s'

        # add necessary inputs for group
        self.add('dv0', IndepVarComp('windDirections', np.zeros(nDirections), units=direction_units), promotes=['*'])
        self.add('dv1', IndepVarComp('windSpeeds', np.zeros(nDirections), units=wind_speed_units), promotes=['*'])
        self.add('dv2', IndepVarComp('windrose_frequencies', np.ones(nDirections)), promotes=['*'])
        self.add('dv3', IndepVarComp('turbineX', np.zeros(nTurbines), units='m'), promotes=['*'])
        self.add('dv4', IndepVarComp('turbineY', np.zeros(nTurbines), units='m'), promotes=['*'])

        # add vars to be seen by MPI and gradient calculations
        self.add('dv5', IndepVarComp('rotorDiameter', np.zeros(nTurbines), units='m'), promotes=['*'])
        self.add('dv6', IndepVarComp('axialInduction', np.zeros(nTurbines)), promotes=['*'])
        self.add('dv7', IndepVarComp('generator_efficiency', np.zeros(nTurbines)), promotes=['*'])
        self.add('dv8', IndepVarComp('air_density', val=1.1716, units='kg/(m*m*m)'), promotes=['*'])

        # add variable tree IndepVarComps
        add_floris_params_IndepVarComps(self, use_rotor_components=use_rotor_components)
        add_gen_params_IdepVarComps(self, datasize=datasize)

        if not use_rotor_components:
            self.add('dv9', IndepVarComp('Ct_in', np.zeros(nTurbines)), promotes=['*'])
            self.add('dv10', IndepVarComp('Cp_in', np.zeros(nTurbines)), promotes=['*'])

        # add components and groups
        self.add('windDirectionsDeMUX', DeMUX(nDirections, units=direction_units))
        self.add('windSpeedsDeMUX', DeMUX(nDirections, units=wind_speed_units))

        pg = self.add('all_directions', ParallelGroup(), promotes=['*'])
        if use_rotor_components:
            for direction_id in np.arange(0, nDirections):
                # print 'assigning direction group %i' % direction_id
                pg.add('direction_group%i' % direction_id,
                       DirectionGroupFLORIS(nTurbines=nTurbines, resolution=resolution, direction_id=direction_id,
                                            use_rotor_components=use_rotor_components, datasize=datasize,
                                            differentiable=differentiable, optimizingLayout=optimizingLayout,
                                            add_IdepVarComps=False, forcefd=forcefd),
                       promotes=['gen_params:*', 'floris_params:*', 'air_density',
                                 'axialInduction', 'generator_efficiency', 'turbineX', 'turbineY',
                                 'yaw%i' % direction_id, 'rotorDiameter', 'velocitiesTurbines%i' % direction_id,
                                 'wt_power%i' % direction_id, 'power%i' % direction_id])#, 'wakeCentersYT', 'wakeDiametersT'])
        else:
            for direction_id in np.arange(0, nDirections):
                # print 'assigning direction group %i' % direction_id
                pg.add('direction_group%i' % direction_id,
                       DirectionGroupFLORIS(nTurbines=nTurbines, resolution=resolution, direction_id=direction_id,
                                            use_rotor_components=use_rotor_components, datasize=datasize,
                                            differentiable=differentiable, add_IdepVarComps=False, forcefd=forcefd),
                       promotes=['Ct_in', 'Cp_in', 'gen_params:*', 'floris_params:*', 'air_density', 'axialInduction',
                                 'generator_efficiency', 'turbineX', 'turbineY', 'yaw%i' % direction_id, 'rotorDiameter',
                                 'velocitiesTurbines%i' % direction_id, 'wt_power%i' % direction_id,
                                 'power%i' % direction_id])#, 'wakeCentersYT', 'wakeDiametersT'])

        self.add('powerMUX', MUX(nDirections, units=power_units))
        self.add('AEPcomp', WindFarmAEP(nDirections), promotes=['*'])

        # connect components
        self.connect('windDirections', 'windDirectionsDeMUX.Array')
        self.connect('windSpeeds', 'windSpeedsDeMUX.Array')
        for direction_id in np.arange(0, nDirections):
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
    root.add('p', floris_velocity(nTurbines=np.array([4])))

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
