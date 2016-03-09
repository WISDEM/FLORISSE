from openmdao.main.api import Assembly
from openmdao.lib.datatypes.api import Array, Bool, Float, VarTree
from openmdao.lib.drivers.api import FixedPointIterator, SLSQPdriver
from Parameters import FLORISParameters

import numpy as np

# ###########    imports for smooth model with analytic gradients    ##################################################
from Analytic_components import AEP
from Analytic_components import dist_const

# ###########    imports from python model    #########################################
from Circle_components import floris_windframe
from Circle_components import floris_wcent_wdiam
from Circle_components import floris_overlap
from Circle_components import floris_power

# ###########    imports for rotor modeling    ########################################################################
from rotor_components import *


class floris_assembly_opt_AEP(Assembly):
    """ Defines the connections between each Component used in the FLORIS model """

    # general input variables
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')
    # Flow property variables
    air_density = Float(iotype='in', units='kg/(m*m*m)', desc='air density in free stream')

    # output
    AEP = Float(iotype='out', units='kW', desc='total windfarm AEP')

    def __init__(self, nTurbines, nDirections, optimize_position=False, nSamples=0, optimize_yaw=False, datasize=0, nSpeeds=False, maxiter=100):

        super(floris_assembly_opt_AEP, self).__init__()

        if nSpeeds == False:
            nSpeeds = nDirections

        self.nTurbines = nTurbines
        self.nSamples = nSamples
        self.nDirections = nDirections
        self.optimize_yaw = optimize_yaw
        self.optimize_position = optimize_position
        self.datasize = datasize
        self.nSpeeds = nSpeeds
        self.maxiter = maxiter

        # wt_layout input variables
        self.add('rotorDiameter', Array(np.zeros(nTurbines), dtype='float', iotype='in', units='m',
                                        desc='rotor diameters of all turbine'))
        self.add('axialInduction', Array(np.zeros(nTurbines), iotype='in', dtype='float',
                                         desc='axial induction of all turbines'))
        self.add('hubHeight', Array(np.zeros(nTurbines), dtype='float', iotype='in', units='m', \
                desc='hub heights of all turbines'))

        # turbine properties for ccblade and pre-calculated controller
        self.add('curve_CP', Array(np.zeros(datasize), iotype='in', desc='pre-calculated CPCT'))
        self.add('curve_CT', Array(np.zeros(datasize), iotype='in', desc='pre-calculated CPCT'))
        self.add('curve_wind_speed', Array(np.zeros(datasize), iotype='in', desc='pre-calculated CPCT'))
        self.add('initVelocitiesTurbines', Array(np.zeros(nTurbines), iotype='in', units='m/s'))
        self.add('generator_efficiency', Array(np.zeros(nTurbines), iotype='in', dtype='float',
                                               desc='generator efficiency of all turbines'))
        self.add('turbineX', Array(np.zeros(nTurbines), iotype='in', dtype='float',
                                   desc='x positions of turbines in original ref. frame'))
        self.add('turbineY', Array(np.zeros(nTurbines), iotype='in', dtype='float',
                                   desc='y positions of turbines in original ref. frame'))
        if optimize_yaw:
            for direction in range(0, nDirections):
                self.add('yaw_%d' % direction, Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                         desc='yaw of each turbine for each direction'))
        else:
            self.add('yaw', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                              desc='yaw of each turbine'))

        # windrose input variables
        self.add('windrose_directions', Array(np.zeros(nDirections), dtype='float', iotype='in',
                                              desc='windrose directions in degrees ccw from east'))
        self.add('windrose_frequencies', Array(np.ones(nDirections), dtype='float', iotype='in',
					desc='windrose frequencies corresponding to windrose_directions'))
        if nSpeeds == 1:
            self.add('windrose_speeds', Float(iotype='in', units='m/s', 
                                          desc='wind speeds for each direction given in windrose_directions'))
        else:
            self.add('windrose_speeds', Array(np.zeros(nDirections), dtype='float', iotype='in', units='m/s', 
                                          desc='wind speeds for each direction given in windrose_directions'))

        # Explicitly size output arrays

        # variables added to test individual components
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='out', units='m',
                                    desc='X positions of turbines in the wind direction reference frame'))
        self.add('turbineYw', Array(np.zeros(nTurbines), iotype='out', units='m',
                                    desc='Y positions of turbines in the wind direction reference frame'))
        self.add('wakeCentersYT', Array(np.zeros(nTurbines), dtype='float', iotype='out', units='m',
                                        desc='centers of the wakes at each turbine'))
        self.add('wakeDiametersT', Array(np.zeros(nTurbines), dtype='float', iotype='out', units='m',
                                         desc='diameters of each of the wake zones for each of the \
                                         wakes at each turbine'))
        self.add('wakeOverlapTRel', Array(np.zeros(nTurbines), dtype='float', iotype='out', units='m',
                                          desc='ratio of overlap area of each zone to rotor area'))

        # standard output
        self.add('velocitiesTurbines_directions', Array(np.zeros([nDirections, nTurbines]), iotype='out', units='m/s',
                                                        dtype='float', desc='effective windspeed at each turbine \
                                                        in each direction ccw from east using direction to'))
        self.add('wt_power_directions', Array(np.zeros([nDirections, nTurbines]), iotype='out', units='kW',
                                              dtype='float', desc='power of each turbine in each direction ccw from \
                                              east using direction to'))
        self.add('power_directions', Array(np.zeros(nDirections), iotype='out', units='kW', desc='total windfarm power \
                                           in each direction ccw from east using direction to'))

        if nSamples>0:
            # flow samples
            self.add('ws_positionX', Array(np.zeros(nSamples), iotype='in', units='m',
                                        desc='X positions of sampling points'))
            self.add('ws_positionY', Array(np.zeros(nSamples), iotype='in', units='m',
                                        desc='Y position of sampling points'))
            self.add('ws_positionZ', Array(np.zeros(nSamples), iotype='in', units='m',
                                        desc='Z position of sampling points'))
            for direction in range(0, nDirections):
                self.add('ws_array_%d' % direction, Array(np.zeros(nSamples), iotype='out', units='m/s', desc='predicted wind speed at sampling points'))

    def configure(self):

        # rename options
        nTurbines = self.nTurbines
        nDirections = self.nDirections
        optimize_position = self.optimize_position
        optimize_yaw = self.optimize_yaw
        datasize = self.datasize
        nSamples = self.nSamples
        nSpeeds = self.nSpeeds
        maxiter = self.maxiter

        # add driver so the workflow is not overwritten later
        if optimize_position or optimize_yaw:
            self.add('driver', SLSQPdriver())

        # add AEP component first so it can be connected to
        F6 = self.add('floris_AEP', AEP(nDirections=nDirections))
        F6.missing_deriv_policy = 'assume_zero'
        self.connect('windrose_frequencies', 'floris_AEP.windrose_frequencies')
        self.connect('floris_AEP.AEP', 'AEP')
        self.connect('floris_AEP.power_directions_out', 'power_directions')

        # set up constraints
        self.add('floris_dist_const', dist_const(nTurbines=nTurbines))
        self.connect('turbineX', 'floris_dist_const.turbineX')
        self.connect('turbineY', 'floris_dist_const.turbineY')

        if nSamples>0:
            samplingNonSampling = ['','Sampling_']
        else:
            samplingNonSampling = ['']

        for i in range(0, nDirections):

            # add fixed point iterator
            self.add('FPIdriver_%d' % i, FixedPointIterator())
            self.add('rotor_CPCT_%d' % i, CPCT_Interpolate(nTurbines=self.nTurbines, datasize=self.datasize))
            CP = 'rotor_CPCT_%d.CP' % i
            CT = 'rotor_CPCT_%d.CT' % i
            CPCT = 'rotor_CPCT_%d' % i

            # add components of floris to assembly
            F2 = self.add('floris_windframe_%d' % i, floris_windframe(nTurbines=nTurbines))
            F2.missing_deriv_policy = 'assume_zero'
            self.add('floris_wcent_wdiam_%d' % i, floris_wcent_wdiam(nTurbines=nTurbines))
            F4 = self.add('floris_overlap_%d' % i, floris_overlap(nTurbines=nTurbines))
            F4.missing_deriv_policy = 'assume_zero'
            self.add('floris_power_%d' % i, floris_power(nTurbines=nTurbines))

            # add visualization components of floris to assembly
            if nSamples>0:
                self.add('Sampling_floris_windframe_%d' % i, floris_windframe(nTurbines=nTurbines, nSamples=nSamples))
                self.add('Sampling_floris_wcent_wdiam_%d' % i, floris_wcent_wdiam(nTurbines=nTurbines, nSamples=nSamples))
                self.add('Sampling_floris_overlap_%d' % i, floris_overlap(nTurbines=nTurbines))
                self.add('Sampling_floris_power_%d' % i, floris_power(nTurbines=nTurbines, nSamples=nSamples))

            # connect inputs to components
            self.connect('curve_CP', 'rotor_CPCT_%d.windSpeedToCPCT.CP' % i)
            self.connect('curve_CT', 'rotor_CPCT_%d.windSpeedToCPCT.CT' % i)
            self.connect('curve_wind_speed', 'rotor_CPCT_%d.windSpeedToCPCT.wind_speed' % i)
            self.connect('parameters.pP', 'rotor_CPCT_%d.pP' % i)

            for ssn in samplingNonSampling:
                self.connect('parameters', ['%sfloris_wcent_wdiam_%d.parameters' % (ssn,i), '%sfloris_power_%d.parameters' % (ssn,i)])
                self.connect('verbose', ['%sfloris_windframe_%d.verbose' % (ssn,i), '%sfloris_wcent_wdiam_%d.verbose' % (ssn,i),
                                     '%sfloris_power_%d.verbose' % (ssn,i)])
                self.connect('turbineX', '%sfloris_windframe_%d.turbineX' % (ssn,i))
                self.connect('turbineY', '%sfloris_windframe_%d.turbineY' % (ssn,i))
                self.connect('rotorDiameter', ['%sfloris_wcent_wdiam_%d.rotorDiameter' % (ssn,i),
                                           '%sfloris_overlap_%d.rotorDiameter' % (ssn,i), '%sfloris_power_%d.rotorDiameter' % (ssn,i)])
                self.connect('axialInduction', '%sfloris_power_%d.axialInduction' % (ssn,i))
                self.connect('generator_efficiency', '%sfloris_power_%d.generator_efficiency' % (ssn,i))

            if nSamples>0:
                # connections needed for visualization
                self.connect('ws_positionX', 'Sampling_floris_windframe_%d.ws_positionX' % i)
                self.connect('ws_positionY', 'Sampling_floris_windframe_%d.ws_positionY' % i)
                self.connect('ws_positionZ', 'Sampling_floris_windframe_%d.ws_positionZ' % i)
                self.connect('hubHeight', 'Sampling_floris_wcent_wdiam_%d.hubHeight' % i)


            if optimize_yaw:
                yawToConnect = 'yaw_%d' % i
            else:
                yawToConnect = 'yaw'

            self.connect(yawToConnect, '%s.yaw' % CPCT)
            for ssn in samplingNonSampling:
                self.connect(yawToConnect, ['%sfloris_wcent_wdiam_%d.yaw' % (ssn,i), '%sfloris_power_%d.yaw' % (ssn,i)])

            for ssn in samplingNonSampling:
                self.connect('air_density', '%sfloris_power_%d.air_density' % (ssn,i))
                self.connect('windrose_directions[%d]' % i, '%sfloris_windframe_%d.wind_direction' % (ssn,i))

            # for satisfying the verbosity in windframe
            for ssn in samplingNonSampling:
                self.connect(CT, '%sfloris_windframe_%d.Ct' % (ssn,i))
                self.connect(CP, '%sfloris_windframe_%d.Cp' % (ssn,i))
                self.connect(yawToConnect, '%sfloris_windframe_%d.yaw' % (ssn,i))
                self.connect('axialInduction', '%sfloris_windframe_%d.axialInduction' % (ssn,i))

            # ############### Connections between components ##################
            # connections from CtCp calculation to other components
            for ssn in samplingNonSampling:
                self.connect(CT, ['%sfloris_wcent_wdiam_%d.Ct' % (ssn,i), '%sfloris_power_%d.Ct' % (ssn,i)])
                self.connect(CP, '%sfloris_power_%d.Cp' % (ssn,i))

                # connections from floris_windframe to floris_wcent_wdiam
                self.connect('%sfloris_windframe_%d.turbineXw' % (ssn,i), '%sfloris_wcent_wdiam_%d.turbineXw' % (ssn,i))
                self.connect('%sfloris_windframe_%d.turbineYw' % (ssn,i), '%sfloris_wcent_wdiam_%d.turbineYw' % (ssn,i))

                # connections from floris_wcent_wdiam to floris_overlap
                self.connect('%sfloris_wcent_wdiam_%d.wakeCentersYT' % (ssn,i), '%sfloris_overlap_%d.wakeCentersYT' % (ssn,i))
                self.connect('%sfloris_wcent_wdiam_%d.wakeDiametersT' % (ssn,i), '%sfloris_overlap_%d.wakeDiametersT' % (ssn,i))

                # connections from floris_windframe to floris_overlap
                self.connect('%sfloris_windframe_%d.turbineXw' % (ssn,i), '%sfloris_overlap_%d.turbineXw' % (ssn,i))
                self.connect('%sfloris_windframe_%d.turbineYw' % (ssn,i), '%sfloris_overlap_%d.turbineYw' % (ssn,i))

                # connections from floris_windframe to floris_power
                self.connect('%sfloris_windframe_%d.turbineXw' % (ssn,i), '%sfloris_power_%d.turbineXw' % (ssn,i))

                # connections from floris_overlap to floris_power
                self.connect('%sfloris_overlap_%d.wakeOverlapTRel' % (ssn,i), '%sfloris_power_%d.wakeOverlapTRel' % (ssn,i))

            # additional connections needed for visualization
            if nSamples>0:
                self.connect('Sampling_floris_windframe_%d.wsw_position' % i, ['Sampling_floris_wcent_wdiam_%d.wsw_position' % i, 'Sampling_floris_power_%d.wsw_position' % i])
                self.connect('Sampling_floris_wcent_wdiam_%d.wakeCentersY' % i, 'Sampling_floris_power_%d.wakeCentersY' % i)
                self.connect('Sampling_floris_wcent_wdiam_%d.wakeCentersZ' % i, 'Sampling_floris_power_%d.wakeCentersZ' % i)
                self.connect('Sampling_floris_wcent_wdiam_%d.wakeDiameters' % i, 'Sampling_floris_power_%d.wakeDiameters' % i)
                self.connect('Sampling_floris_power_%d.ws_array' % i, 'ws_array_%d' % i)

            # connections from floris_power to floris_AEP
            self.connect('floris_power_%d.power' % i, 'floris_AEP.power_directions[%d]' % i)
            # #################################################################

            # add to workflow
            exec("self.FPIdriver_%d.workflow.add(['rotor_CPCT_%d', 'floris_windframe_%d', \
                 'floris_wcent_wdiam_%d', 'floris_overlap_%d', 'floris_power_%d'])" % (i, i, i, i, i, i))
            exec("self.FPIdriver_%d.add_parameter('rotor_CPCT_%d.wind_speed_hub', low=0., high=100.)" % (i, i))
            exec("self.FPIdriver_%d.add_constraint('rotor_CPCT_%d.wind_speed_hub = \
                  floris_power_%d.velocitiesTurbines')" % (i, i, i))
            self.driver.workflow.add('FPIdriver_%d' % i)
            if nSamples>0:
                self.driver.workflow.add(['Sampling_floris_windframe_%d' % i,
                                          'Sampling_floris_wcent_wdiam_%d' % i, 'Sampling_floris_overlap_%d' % i, 'Sampling_floris_power_%d' % i])

        if nSpeeds>1:
            for i in range(0, nSpeeds):
                for ssn in samplingNonSampling:
                    self.connect('windrose_speeds[%d]' % i, '%sfloris_power_%d.wind_speed' % (ssn,i))
                    self.connect('windrose_speeds[%d]' % i, '%sfloris_windframe_%d.wind_speed' % (ssn,i))
        else:
            for i in range(0, nDirections):
                for ssn in samplingNonSampling:
                    self.connect('windrose_speeds', '%sfloris_power_%d.wind_speed' % (ssn,i))
                    self.connect('windrose_speeds', '%sfloris_windframe_%d.wind_speed' % (ssn,i))

        # add AEP calculations to workflow
        self.driver.workflow.add(['floris_AEP', 'floris_dist_const'])
        if optimize_position or optimize_yaw:
            # set up driver
            self.driver.iprint = 3
            self.driver.accuracy = 1.0e-12
            self.driver.maxiter = maxiter
            self.driver.add_objective('-floris_AEP.AEP')
            if optimize_position:
                self.driver.add_parameter('turbineX', low=7*126.4, high=np.sqrt(self.nTurbines)*7*126.4)
                self.driver.add_parameter('turbineY', low=7*126.4, high=np.sqrt(self.nTurbines)*7*126.4)
                self.driver.add_constraint('floris_dist_const.separation > 2*rotorDiameter[0]')
            if optimize_yaw:
                for direction in range(0, self.nDirections):
                    self.driver.add_parameter('yaw_%d' % direction, low=-30., high=30., scaler=1.)
