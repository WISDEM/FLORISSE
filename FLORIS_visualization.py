from openmdao.main.api import Assembly
from openmdao.lib.datatypes.api import Array, Bool, Float, VarTree
from openmdao.lib.drivers.api import SLSQPdriver
from Parameters import FLORISParameters

# ###########    imports for discontinuous (original) model    ##########
from Original_components import floris_windframe
from Original_components import floris_wcent_wdiam
from Original_components import floris_overlap
from Original_components import floris_power

# ###########    imports for smooth model w/o gradients    ##############
# from RegSmooth_components import floris_windframe
# from RegSmooth_components import floris_wcent_wdiam
# from RegSmooth_components import floris_overlap
# from RegSmooth_components import floris_power

# ###########    imports for smooth model with analytic gradients    ####
# from Analytic_components import floris_windframe

# ###########    imports for AlgoPy components    #######################
# from AlgoPy_components import floris_windframe #-- not available
# from AlgoPy_components import floris_wcent_wdiam #-- not available
# from AlgoPy_components import floris_overlap #-- not available
# from AlgoPy_components import floris_power #-- unfinished

# ###########    imports for Tapenade components    #####################
# from Tapenade_components import floris_windframe -- not available
# from Tapenade_components import floris_wcent_wdiam
# from Tapenade_components import floris_overlap
# from Tapenade_components import floris_power # velocitiesTurbines gradients are incorrect
#

# ###########    imports for Fortran components (not gradient)    #####################
# from Fortran_components import floris_windframe -- not available
# from Fortran_components import floris_wcent_wdiam
# from Fortran_components import floris_overlap
# from Fortran_components import floris_power
#figure_1b

class floris_assembly_opt(Assembly):
    """ Defines the connections between each Component used in the FLORIS model """
    #
    # original input variables in Pieter's OpenMDAO stand-alone version of FLORIS
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    # input variables added so I don't have to use WISDEM while developing gradients
    # position = Array(iotype='in', desc='position of turbines in original ref. frame')
    # turbineX = Array(iotype='in', desc='x positions of turbines in original ref. frame')
    # turbineY = Array(iotype='in', desc='y positions of turbines in original ref. frame')
    ws_position = Array(iotype='in', units='m', desc='position where you want measurements in ref. frame')
    rotorDiameter = Array(dtype='float', iotype='in', units='m', desc='rotor diameters of all turbine')
    rotorArea = Array(iotype='in', dtype='float', units='m*m', desc='rotor area of all turbines')
    axialInduction = Array(iotype='in', dtype='float', desc='axial induction of all turbines')
    Ct = Array(iotype='in', desc='Thrust coefficient for all turbines')
    Cp = Array(iotype='in', dtype='float', desc='power coefficient for all turbines')
    generator_efficiency = Array(iotype='in', dtype='float', desc='generator efficiency of all turbines')
    yaw = Array(iotype='in', desc='yaw of each turbine')

    # Flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    air_density = Float(iotype='in', units='kg/(m*m*m)', desc='air density in free stream')
    wind_direction = Float(iotype='in', units='deg', desc='overall wind direction for wind farm')

    # original output variables in Pieter's OpenMDAO stand-alone version of FLORIS
    velocitiesTurbines = Array(iotype='out', units='m/s')

    # output variables added so I don't have to use WISDEM while developing gradients
    wt_power = Array(iotype='out', units='kW')
    ws_array = Array(iotype='out', units='m/s', desc='wind speed at measurement locations')

    # variables added to test individual components
    turbineXw = Array(iotype='out', units='m', desc='X positions of turbines in the wind direction reference frame')
    turbineYw = Array(iotype='out', units='m', desc='Y positions of turbines in the wind direction reference frame')
    wakeCentersYT = Array(dtype='float', iotype='out', units='m', desc='centers of the wakes at each turbine')
    wakeDiametersT = Array(dtype='float', iotype='out', units='m', desc='diameters of each of the wake zones for each of the wakes at each turbine')
    wakeOverlapTRel = Array(dtype='float', iotype='out', units='m', desc='ratio of overlap area of each zone to rotor area')

    # testing
    # p_near0 = Float(iotyp='out', desc='upwind location of diameter spline in rotor diameters')

    # final output
    power = Float(iotype='out', units='kW', desc='total windfarm power')


    def configure(self):

        # add components to floris assembly
        self.add('floris_windframe', floris_windframe())
        self.add('floris_wcent_wdiam', floris_wcent_wdiam())
        self.add('floris_overlap', floris_overlap())
        self.add('floris_power', floris_power())

        # add driver to floris assembly
        self.driver.workflow.add(['floris_windframe', 'floris_wcent_wdiam', 'floris_overlap', 'floris_power'])

        # added for gradient testing
        self.add('driver', SLSQPdriver())
        self.driver.iprint = 3
        self.driver.accuracy = 1.0e-6
        self.driver.maxiter = 100
        self.driver.add_objective('-floris_power.power')
        self.driver.add_parameter('floris_windframe.turbineX', low=0., high=3000.)
        self.driver.add_parameter('floris_windframe.turbineY', low=0., high=3000.)
        self.driver.add_parameter('yaw', low=-20., high=20.)

        # connect inputs to components
        self.connect('parameters', ['floris_wcent_wdiam.parameters', 'floris_power.parameters'])
        self.connect('verbose', ['floris_windframe.verbose', 'floris_wcent_wdiam.verbose', 'floris_power.verbose'])
        # self.connect('position', 'floris_windframe.position')
        # self.connect('turbineX', 'floris_windframe.turbineX')
        # self.connect('turbineY', 'floris_windframe.turbineY')
        self.connect('ws_position', 'floris_windframe.ws_position')
        self.connect('rotorDiameter', ['floris_wcent_wdiam.rotorDiameter', 'floris_overlap.rotorDiameter', 'floris_power.rotorDiameter'])
        self.connect('rotorArea', ['floris_overlap.rotorArea', 'floris_power.rotorArea'])
        self.connect('axialInduction', 'floris_power.axialInduction')
        self.connect('Ct', ['floris_wcent_wdiam.Ct', 'floris_power.Ct'])
        self.connect('Cp', 'floris_power.Cp')
        self.connect('generator_efficiency', 'floris_power.generator_efficiency')
        self.connect('yaw', ['floris_wcent_wdiam.yaw', 'floris_power.yaw'])
        self.connect('wind_speed', 'floris_power.wind_speed')
        self.connect('air_density', 'floris_power.air_density')
        self.connect('wind_direction', 'floris_windframe.wind_direction')

        # for satisfying the verbosity in windframe
        self.connect('Cp', 'floris_windframe.Cp')
        self.connect('Ct', 'floris_windframe.Ct')
        # self.connect('yaw', 'floris_windframe.yaw')
        self.connect('wind_speed', 'floris_windframe.wind_speed')
        self.connect('axialInduction', 'floris_windframe.axialInduction')

        # ############### Connections between components ##################
        # connections from floris_windframe to floris_wcent_wdiam
        self.connect("floris_windframe.turbineXw", "floris_wcent_wdiam.turbineXw")
        self.connect("floris_windframe.turbineYw", "floris_wcent_wdiam.turbineYw")
        self.connect("floris_windframe.wsw_position", "floris_wcent_wdiam.wsw_position")

        # connections from floris_wcent_wdiam to floris_overlap
        self.connect("floris_wcent_wdiam.wakeCentersYT", "floris_overlap.wakeCentersYT")
        self.connect("floris_wcent_wdiam.wakeDiametersT", "floris_overlap.wakeDiametersT")

        # connections from floris_windframe to floris_overlap
        self.connect("floris_windframe.turbineXw", "floris_overlap.turbineXw")
        self.connect("floris_windframe.turbineYw", "floris_overlap.turbineYw")

        # connections from floris_windframe to floris_power
        self.connect('floris_windframe.turbineXw', 'floris_power.turbineXw')
        self.connect('floris_windframe.wsw_position', 'floris_power.wsw_position')


        # test
        # self.connect('floris_wcent_wdiam.p_near0', 'floris_overlap.p_near0')

        # connections from floris_wcent_wdiam to floris_power
        self.connect("floris_wcent_wdiam.wakeCentersY", "floris_power.wakeCentersY")
        self.connect("floris_wcent_wdiam.wakeDiameters", "floris_power.wakeDiameters")

        # connections from floris_overlap to floris_power
        self.connect("floris_overlap.wakeOverlapTRel", "floris_power.wakeOverlapTRel")
        # #################################################################

        # output connections
        self.connect("floris_power.velocitiesTurbines", "velocitiesTurbines")
        self.connect("floris_power.wt_power", "wt_power")
        self.connect("floris_power.power", "power")
        self.connect("floris_power.ws_array", "ws_array")

        # outputs for testing only
        self.connect("floris_windframe.turbineXw", "turbineXw")
        self.connect("floris_windframe.turbineYw", "turbineYw")
        self.connect("floris_wcent_wdiam.wakeCentersYT", "wakeCentersYT")
        self.connect("floris_wcent_wdiam.wakeDiametersT", "wakeDiametersT")
        self.connect("floris_overlap.wakeOverlapTRel", "wakeOverlapTRel")
        # self.connect("floris_wcent_wdiam.p_near0", "p_near0")


class floris_assembly(Assembly):
    """ Defines the connections between each Component used in the FLORIS model """
    #
    # original input variables in Pieter's OpenMDAO stand-alone version of FLORIS
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    # input variables added so I don't have to use WISDEM while developing gradients
    # position = Array(iotype='in', desc='position of turbines in original ref. frame')
    turbineX = Array(iotype='in', desc='x positions of turbines in original ref. frame')
    turbineY = Array(iotype='in', desc='y positions of turbines in original ref. frame')
    ws_position = Array(iotype='in', units='m', desc='position where you want measurements in ref. frame')
    rotorDiameter = Array(dtype='float', iotype='in', units='m', desc='rotor diameters of all turbine')
    rotorArea = Array(iotype='in', dtype='float', units='m*m', desc='rotor area of all turbines')
    axialInduction = Array(iotype='in', dtype='float', desc='axial induction of all turbines')
    Ct = Array(iotype='in', desc='Thrust coefficient for all turbines')
    Cp = Array(iotype='in', dtype='float', desc='power coefficient for all turbines')
    generator_efficiency = Array(iotype='in', dtype='float', desc='generator efficiency of all turbines')
    yaw = Array(iotype='in', desc='yaw of each turbine')

    # Flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    air_density = Float(iotype='in', units='kg/(m*m*m)', desc='air density in free stream')
    wind_direction = Float(iotype='in', units='deg', desc='overall wind direction for wind farm')

    # original output variables in Pieter's OpenMDAO stand-alone version of FLORIS
    velocitiesTurbines = Array(iotype='out', units='m/s')

    # output variables added so I don't have to use WISDEM while developing gradients
    wt_power = Array(iotype='out', units='kW')
    ws_array = Array(iotype='out', units='m/s', desc='wind speed at measurement locations')

    # variables added to test individual components
    turbineXw = Array(iotype='out', units='m', desc='X positions of turbines in the wind direction reference frame')
    turbineYw = Array(iotype='out', units='m', desc='Y positions of turbines in the wind direction reference frame')
    wakeCentersYT = Array(dtype='float', iotype='out', units='m', desc='centers of the wakes at each turbine')
    wakeDiametersT = Array(dtype='float', iotype='out', units='m', desc='diameters of each of the wake zones for each of the wakes at each turbine')
    wakeOverlapTRel = Array(dtype='float', iotype='out', units='m', desc='ratio of overlap area of each zone to rotor area')

    # testing
    # p_near0 = Float(iotyp='out', desc='upwind location of diameter spline in rotor diameters')

    # final output
    power = Float(iotype='out', units='kW', desc='total windfarm power')


    def configure(self):

        # add components to floris assembly
        self.add('floris_windframe', floris_windframe())
        self.add('floris_wcent_wdiam', floris_wcent_wdiam())
        self.add('floris_overlap', floris_overlap())
        self.add('floris_power', floris_power())

        # add driver to floris assembly
        self.driver.workflow.add(['floris_windframe', 'floris_wcent_wdiam', 'floris_overlap', 'floris_power'])

        # added for gradient testing
        # self.add('driver', SLSQPdriver())
        # self.driver.iprint = 0
        # self.driver.add_objective('-floris_power.power')
        # self.driver.add_parameter('floris_windframe.turbineX', low=0., high=3000.)
        # self.driver.add_parameter('floris_windframe.turbineY', low=0., high=3000.)

        # connect inputs to components
        self.connect('parameters', ['floris_wcent_wdiam.parameters', 'floris_power.parameters'])
        self.connect('verbose', ['floris_windframe.verbose', 'floris_wcent_wdiam.verbose', 'floris_power.verbose'])
        # self.connect('position', 'floris_windframe.position')
        self.connect('turbineX', 'floris_windframe.turbineX')
        self.connect('turbineY', 'floris_windframe.turbineY')
        self.connect('ws_position', 'floris_windframe.ws_position')
        self.connect('rotorDiameter', ['floris_wcent_wdiam.rotorDiameter', 'floris_overlap.rotorDiameter', 'floris_power.rotorDiameter'])
        self.connect('rotorArea', ['floris_overlap.rotorArea', 'floris_power.rotorArea'])
        self.connect('axialInduction', 'floris_power.axialInduction')
        self.connect('Ct', ['floris_wcent_wdiam.Ct', 'floris_power.Ct'])
        self.connect('Cp', 'floris_power.Cp')
        self.connect('generator_efficiency', 'floris_power.generator_efficiency')
        self.connect('yaw', ['floris_wcent_wdiam.yaw', 'floris_power.yaw'])
        self.connect('wind_speed', 'floris_power.wind_speed')
        self.connect('air_density', 'floris_power.air_density')
        self.connect('wind_direction', 'floris_windframe.wind_direction')


        # for satisfying the verbosity in windframe
        self.connect('Cp', 'floris_windframe.Cp')
        self.connect('Ct', 'floris_windframe.Ct')
        # self.connect('yaw', 'floris_windframe.yaw')
        self.connect('wind_speed', 'floris_windframe.wind_speed')
        self.connect('axialInduction', 'floris_windframe.axialInduction')

        # ############### Connections between components ##################
        # connections from floris_windframe to floris_wcent_wdiam
        self.connect("floris_windframe.turbineXw", "floris_wcent_wdiam.turbineXw")
        self.connect("floris_windframe.turbineYw", "floris_wcent_wdiam.turbineYw")
        self.connect("floris_windframe.wsw_position", "floris_wcent_wdiam.wsw_position")

        # connections from floris_wcent_wdiam to floris_overlap
        self.connect("floris_wcent_wdiam.wakeCentersYT", "floris_overlap.wakeCentersYT")
        self.connect("floris_wcent_wdiam.wakeDiametersT", "floris_overlap.wakeDiametersT")

        # connections from floris_windframe to floris_overlap
        self.connect("floris_windframe.turbineXw", "floris_overlap.turbineXw")
        self.connect("floris_windframe.turbineYw", "floris_overlap.turbineYw")

        # connections from floris_windframe to floris_power
        self.connect('floris_windframe.turbineXw', 'floris_power.turbineXw')
        self.connect('floris_windframe.wsw_position', 'floris_power.wsw_position')


        # test
        # self.connect('floris_wcent_wdiam.p_near0', 'floris_overlap.p_near0')

        # connections from floris_wcent_wdiam to floris_power
        self.connect("floris_wcent_wdiam.wakeCentersY", "floris_power.wakeCentersY")
        self.connect("floris_wcent_wdiam.wakeDiameters", "floris_power.wakeDiameters")

        # connections from floris_overlap to floris_power
        self.connect("floris_overlap.wakeOverlapTRel", "floris_power.wakeOverlapTRel")
        # #################################################################

        # output connections
        self.connect("floris_power.velocitiesTurbines", "velocitiesTurbines")
        self.connect("floris_power.wt_power", "wt_power")
        self.connect("floris_power.power", "power")
        self.connect("floris_power.ws_array", "ws_array")

        # outputs for testing only
        self.connect("floris_windframe.turbineXw", "turbineXw")
        self.connect("floris_windframe.turbineYw", "turbineYw")
        self.connect("floris_wcent_wdiam.wakeCentersYT", "wakeCentersYT")
        self.connect("floris_wcent_wdiam.wakeDiametersT", "wakeDiametersT")
        self.connect("floris_overlap.wakeOverlapTRel", "wakeOverlapTRel")
        # self.connect("floris_wcent_wdiam.p_near0", "p_near0")


