from openmdao.main.api import Assembly, Component
from openmdao.lib.datatypes.api import Array, Float, Bool, Int, List, Str, VarTree

from rotorse.rotoraerodefaults import CCBlade
from rotorse.rotoraero import Coefficients
from floris import FLORIS
from fusedwind.turbine.turbine_vt import AeroelasticHAWTVT
from fusedwind.plant_flow.comp import GenericWindFarm
from openmdao.lib.drivers.api import BroydenSolver, CaseIteratorDriver, FixedPointIterator
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout
from floris import FLORIS

import os
import numpy as np

class AeroelasticHAWTVT_CCBlade_floris(AeroelasticHAWTVT):
    """ AeroelasticHAWTVT with extra CCblade inputs and control inputs"""

    # inherits from AeroelasticHAWTVT:
    # turbine_name
    # rotor_diameter
    # nblades = 3 --> B
    # hub_height = 90.0 --> hubHt
    # tilt_angle --> tilt
    # cone_angle --> precone
    # hub_radius --> Rhub

    rotor_area = Float(iotype='in', units='m*m', desc='rotor swept area')
    position = Array(iotype='in', units='m', desc='position')
    name = Str(iotype='in', desc='turbine tag in wind plant')
    turbineName = Str(iotype='in', desc='turbine tag in wind plant')
    # double names needed, that needs to be fixed! Related to name being overwritten by signal name (e.g. turbineIn)

    # flow
    wind_speed_hub = Float(iotype='in', units='m/s', desc='hub height wind speed') # Uhub
    air_density = Float(1.225, iotype='in', units='kg/m**3', desc='density of air') # rho
    dynamic_viscosity = Float(1.81206e-5, iotype='in', units='kg/(m*s)', desc='dynamic viscosity of air') # mu
    shear_exponent = Float(0.2, iotype='in', desc='shear exponent') # shearExp

    # CCblade inputs
    tip_radius = Float(iotype='in', units='m', desc='tip radius') # Rtip
    radial_locations = Array(iotype='in', units='m', desc='radial locations where blade is defined (should be increasing and not go all the way to hub or tip)') #r
    chord = Array(iotype='in', units='m', desc='chord length at each section')
    theta = Array(iotype='in', units='deg', desc='twist angle at each section (positive decreases angle of attack)')
    precurve = Array(iotype='in', units='m', desc='precurve at each section')
    precurveTip = Float(0.0, iotype='in', units='m', desc='precurve at tip')
    airfoil_files = List(Str, iotype='in', desc='names of airfoil file')
    nSector = Int(4, iotype='in', desc='number of sectors to divide rotor face into in computing thrust and power')
    tiploss = Bool(True, iotype='in', desc='include Prandtl tip loss model')
    hubloss = Bool(True, iotype='in', desc='include Prandtl hub loss model')
    wakerotation = Bool(True, iotype='in', desc='include effect of wake rotation (i.e., tangential induction factor is nonzero)')
    usecd = Bool(True, iotype='in', desc='use drag coefficient in computing induction factors')

    # control DOFs
    yaw = Float(iotype='in', desc='yaw error', units='deg')
    pitch = Float(iotype='in', desc='blade pitch angle', units='deg')
    rotor_speed = Float(iotype='in', desc='rotor speed', units='rpm') #Omega

    # CCblade outputs
    power = Float(iotype='out', desc='rotor power', units='W')
    thrust = Float(iotype='out', desc='aerodynamic thrust on rotor', units='N')
    torque = Float(iotype='out', desc='aerodynamic torque on rotor', units='N*m')

    # coefficients
    CP = Float(iotype='out', desc='power coefficient')
    CT = Float(iotype='out', desc='thrust coefficient')
    CQ = Float(iotype='out', desc='torque coefficient')
    axial_induction = Float(iotype='out', desc='axial induction')

    # control
    opt_tsr = Float(iotype='in', desc='TSR to track in below-rated conditions')
    rated_generator_speed = Float(iotype='in', units='rpm', desc='rated generator speed')
    gearbox_ratio = Float(iotype='in', desc='gearbox ratio')
    rated_power = Float(iotype='in', units='W', desc='rated electric power')
    generator_efficiency = Float(iotype='in', desc='generator efficiency')

    # FLORIS outputs
    wind_speed_eff = Float(iotype='out', desc='effective wind speed', units='m/s')


class CCBladeCoefficients(Assembly):
    """ Couples components with CCBlade so that it outputs CP and CT coefficients and axial-induction factor
    takes AeroelasticHAWTVT_CCBlade as input
    """

    turbineIn = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='in')
    turbineOut = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='out')

    def configure(self):
        """ Creates a new Assembly containing a Paraboloid component"""

        self.add("CCBlade", CCBlade())
        self.add("Coefficients", Coefficients())
        self.add("CTtoAxialInd", CTtoAxialInd())

        self.driver.workflow.add(['CCBlade', 'Coefficients', 'CTtoAxialInd'])
        self.CCBlade.run_case = 'power'

        self.connect("turbineIn.nblades", "CCBlade.B")
        self.connect("turbineIn.hub_height", "CCBlade.hubHt")
        self.connect("turbineIn.tilt_angle", "CCBlade.tilt")
        self.connect("turbineIn.cone_angle", "CCBlade.precone")
        self.connect("turbineIn.hub_radius", "CCBlade.Rhub")
        self.connect("turbineIn.tip_radius", "CCBlade.Rtip")
        self.connect("turbineIn.wind_speed_hub", "CCBlade.Uhub[0]")
        self.connect("turbineIn.air_density", "CCBlade.rho")
        self.connect("turbineIn.dynamic_viscosity", "CCBlade.mu")
        self.connect("turbineIn.shear_exponent", "CCBlade.shearExp")
        self.connect("turbineIn.radial_locations", "CCBlade.r")
        self.connect("turbineIn.chord", "CCBlade.chord")
        self.connect("turbineIn.theta", "CCBlade.theta")
        self.connect("turbineIn.precurve", "CCBlade.precurve")
        self.connect("turbineIn.precurveTip", "CCBlade.precurveTip")
        self.connect("turbineIn.airfoil_files", "CCBlade.airfoil_files")
        self.connect("turbineIn.nSector", "CCBlade.nSector")
        self.connect("turbineIn.tiploss", "CCBlade.tiploss")
        self.connect("turbineIn.hubloss", "CCBlade.hubloss")
        self.connect("turbineIn.wakerotation", "CCBlade.wakerotation")
        self.connect("turbineIn.usecd", "CCBlade.usecd")
        self.connect("turbineIn.yaw", "CCBlade.yaw")
        self.connect("turbineIn.pitch", "CCBlade.pitch[0]")
        self.connect("turbineIn.rotor_speed", "CCBlade.Omega[0]")

        self.connect("CCBlade.P", "Coefficients.P")
        self.connect("CCBlade.T", "Coefficients.T")
        self.connect("CCBlade.Q", "Coefficients.Q")

        self.connect("CCBlade.P[0]", "turbineOut.power")
        self.connect("CCBlade.T[0]", "turbineOut.thrust")
        self.connect("CCBlade.Q[0]", "turbineOut.torque")

        self.connect("turbineIn.wind_speed_hub", "Coefficients.V[0]")
        self.connect("turbineIn.tip_radius", "Coefficients.R")
        self.connect("turbineIn.air_density", "Coefficients.rho")

        self.connect("Coefficients.CT[0]", "CTtoAxialInd.CT")
        self.connect("Coefficients.CP[0]", "turbineOut.CP")
        self.connect("Coefficients.CT[0]", "turbineOut.CT")
        self.connect("Coefficients.CQ[0]", "turbineOut.CQ")
        self.connect("CTtoAxialInd.axial_induction", "turbineOut.axial_induction")

        # directly feed-through rest of turbine properties In to Out
        attributes = ("rotor_area","position","name","turbine_name","rotor_diameter","nblades", "hub_height",
                      "tilt_angle", "cone_angle", "hub_radius", "tip_radius", "wind_speed_hub","air_density",
                      "dynamic_viscosity", "shear_exponent", "radial_locations", "chord", "theta", "precurve",
                      "precurveTip", "airfoil_files", "nSector", "tiploss", "hubloss", "wakerotation", "usecd",
                      "pitch", "yaw", "rotor_speed", "turbineName")
        for attribute in attributes:
            self.connect("turbineIn.%s" % attribute, "turbineOut.%s" % attribute)


class CTtoAxialInd(Component):
    """Convert thrust coefficient to axial induction factor"""

    CT = Float(iotype='in')
    axial_induction = Float(iotype='out')

    def execute(self):
        self.axial_induction = 0.5*(1-np.sqrt(1-self.CT))


class windSpeedDistributor(Component):
    wind_speed_in = Array(iotype='in', units='m/s', desc='effective wind speed')
    wind_speed_out = Array(iotype='out', units='m/s', desc='effective wind speed')

    def execute(self):
        self.wind_speed_out = self.wind_speed_in


class CCBladePlant(Assembly):

    listOfTurbinesIn = Array(iotype='in')
    listOfTurbinesOut = Array(iotype='out')

    def configure(self):
        self.add('CCBlade', CCBladeCoefficients())
        self.add('driver', CaseIteratorDriver())
        self.driver.workflow.add(['CCBlade'])
        self.driver.add_parameter('CCBlade.turbineIn')
        self.driver.add_response('CCBlade.turbineOut')
        self.connect('listOfTurbinesIn', 'driver.case_inputs.CCBlade.turbineIn')
        self.connect('driver.case_outputs.CCBlade.turbineOut', 'listOfTurbinesOut')


class AeroelasticHAWTVT_floris(AeroelasticHAWTVT):
    # inherits from AeroelasticHAWTVT:
    # turbine_name
    # rotor_diameter
    # and other parameters that are not used:
    #   nblades
    #   hub_height
    #   tilt_angle
    #   cone_angle
    #   hub_radius

    CP = Float(iotype='in',desc='power coefficient')
    axial_induction = Float(iotype='in',desc='axial induction')
    position = Array(iotype='in', units='m', desc='position')
    rotor_area = Float(iotype='in', units='m*m', desc='rotor swept area')
    name = Str(iotype='in', desc='turbine tag in wind plant')
    yaw = Float(iotype='in', desc='yaw error', units='deg')
    power = Float(units='W', desc='power')
    wind_speed_eff = Float(iotype='out', units='m/s', desc='effective wind speed')


class florisPlant(Component):

    listOfTurbines = Array(iotype='in')
    #yaw = Array(iotype='in')
    wind_speed = Float(iotype='in')
    wind_direction = Float(iotype='in')
    air_density = Float(iotype='in')

    velocitiesTurbines = Array(iotype='out')

    def execute(self):

        self.florisWindPlant = FLORIS()

        for turbineI in range(0, len(self.listOfTurbines)):
            wt = AeroelasticHAWTVT_floris()

            wt.name = self.listOfTurbines[turbineI].turbineName
            attributes = ('turbine_name', 'rotor_diameter', 'CP', 'axial_induction', 'position', 'rotor_area', 'yaw')
            for attribute in attributes:
                setattr(wt, attribute, getattr(self.listOfTurbines[turbineI], attribute))
            #wt.yaw = yaw[turbineI]

            self.florisWindPlant.wt_layout.add_wt(wt)

        self.florisWindPlant.wind_speed = self.wind_speed
        self.florisWindPlant.wind_direction = self.wind_direction
        self.florisWindPlant.air_density = self.air_density

        self.florisWindPlant.run()

        self.velocitiesTurbines = self.florisWindPlant.velocitiesTurbines


class ControllerPlant(Assembly):

    listOfTurbinesIn = Array(iotype='in')
    listOfTurbinesOut = Array(iotype='out')
    verbose = Bool(False, iotype='in', desc='verbosity of controller, False is no output')

    def configure(self):
        self.add('PowerSpeedController', PowerSpeedController())
        self.add('driver', CaseIteratorDriver())
        self.driver.workflow.add(['PowerSpeedController'])
        self.driver.add_parameter('PowerSpeedController.turbineIn')
        self.driver.add_response('PowerSpeedController.turbineOut')
        self.connect('verbose', 'PowerSpeedController.verbose')
        self.connect('listOfTurbinesIn', 'driver.case_inputs.PowerSpeedController.turbineIn')
        self.connect('driver.case_outputs.PowerSpeedController.turbineOut', 'listOfTurbinesOut')


class PowerSpeedController(Component):
    turbineIn = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='in')
    turbineOut = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='out')
    verbose = Bool(False, iotype='in', desc='verbosity of controller, False is no output')

    def execute(self):
        turbine = self.turbineIn

        # set optimal rotor speed and zero pitch initially
        wind_speed_ax = np.cos(turbine.yaw*np.pi/180.0)*turbine.wind_speed_hub
        turbine.rotor_speed = wind_speed_ax*turbine.opt_tsr/turbine.tip_radius * 30.0/np.pi
        turbine.pitch = 0.0
        if self.verbose:
            if hasattr(turbine,'turbineName'):
                TopDelimiter = "_____Controller %s_______________________" % turbine.turbineName
            else:
                TopDelimiter = "_____Controller_______________________________"
            print TopDelimiter
            print "axial component wind speed %s" % wind_speed_ax

        # define constraints
        rated_rotor_speed = turbine.rated_generator_speed / turbine.gearbox_ratio
        rated_mechanical_power = turbine.rated_power/turbine.generator_efficiency
        rated_generator_torque = rated_mechanical_power/(turbine.rated_generator_speed*np.pi/30.0)
        rated_rotor_torque = rated_generator_torque*turbine.gearbox_ratio

        # if optimal rotor speed exceeds rated rotor speed, use rated rotor speed instead
        if turbine.rotor_speed > rated_rotor_speed:
            turbine.rotor_speed = rated_rotor_speed
            region = 2.5
        else:
            region = 2

        if self.verbose:
            print "rotor speed %s" % turbine.rotor_speed

        # calculate generator torque
        CCBlade = CCBladeCoefficients()
        CCBlade.turbineIn = turbine
        CCBlade.turbineIn.wind_speed_hub = wind_speed_ax
        CCBlade.turbineIn.yaw = 0.0
        CCBlade.run()

        # if above-rated torque, calculate pitch_angle that will result in rated torque
        if CCBlade.turbineOut.torque>rated_rotor_torque:
            calcPitch = CalculateAboveRatedPitch()
            calcPitch.turbine = turbine
            calcPitch.rated_rotor_torque = rated_rotor_torque
            calcPitch.turbine.wind_speed_hub = wind_speed_ax
            calcPitch.turbine.yaw = 0.0
            calcPitch.run()
            turbine.pitch = calcPitch.CCBlade.turbineIn.pitch
            region = 3

        if self.verbose:
            print "pitch %s" % turbine.pitch
            print "region %s" % region
            print "-"*len(TopDelimiter)

        self.turbineOut = turbine


class CalculateAboveRatedPitch(Assembly):
    turbine = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='in')
    rated_rotor_torque = Float(iotype='in', units='N*m', desc='rated rotor torque')

    def configure(self):
        self.add("CCBlade", CCBladeCoefficients())
        self.add("CQfromTorque", CQfromTorque())

        self.add('optimizer', BroydenSolver())

        self.driver.workflow.add(['CQfromTorque', 'optimizer'])
        self.connect('rated_rotor_torque', 'CQfromTorque.torque')
        self.connect('turbine.air_density', 'CQfromTorque.air_density')
        self.connect('turbine.rotor_area', 'CQfromTorque.rotor_area')
        self.connect('turbine.wind_speed_hub', 'CQfromTorque.wind_speed')
        self.connect('turbine.tip_radius', 'CQfromTorque.tip_radius')

        # connect all CCBlade turbine properties except for pitch, for which we use optimizer
        attributes = ("nblades", "hub_height", "tilt_angle", "cone_angle", "hub_radius", "tip_radius", "wind_speed_hub",
                      "air_density","dynamic_viscosity","shear_exponent","radial_locations","chord","theta","precurve",
                      "precurveTip","airfoil_files","nSector","tiploss","hubloss","wakerotation","usecd","yaw","rotor_speed")
        for attribute in attributes:
            self.connect("turbine.%s" % attribute, "CCBlade.turbineIn.%s" % attribute)

        self.optimizer.itmax = 20
        self.optimizer.tol = 0.000001
        self.optimizer.workflow.add(['CCBlade'])
        self.optimizer.add_parameter('CCBlade.turbineIn.pitch', low=0.0, high=90.0, start=10.0)
        self.optimizer.add_constraint('CCBlade.turbineOut.CQ = CQfromTorque.CQ')


class CQfromTorque(Component):
    torque = Float(iotype='in', units='N*m', desc='rotor torque')
    air_density = Float(iotype='in', units='kg/m**3', desc='density of air')
    rotor_area = Float(iotype='in', units='m*m', desc='rotor swept area')
    tip_radius = Float(iotype='in', units='m', desc='tip radius')
    wind_speed = Float(iotype='in', units='m/s', desc='effective wind speed')
    CQ = Float(iotype='out', desc='torque coefficient')

    def execute(self):
        self.CQ = self.torque / (0.50 * self.air_density * self.rotor_area * self.tip_radius * (self.wind_speed ** 2))

def constructCoupledFLORIS_CCBlade_control(number_of_turbines, freestream_wind_speed=103, startpoint_wind_speed=10, tolerance=0.00001, max_iteration=100):
    """Returns the assembly of the FLORIS model with a coupled CC-Blade model to calculate axial induction, and a power and speed controller for each turbine"""

    FLORIS_CCBlade_control = Assembly()
    FLORIS_CCBlade_control.add('driver', FixedPointIterator())  # Note that FixedPointIterator is top-level (driver)

    FLORIS_CCBlade_control.add('windSpeedDistributor', windSpeedDistributor())
    FLORIS_CCBlade_control.add('control', ControllerPlant())
    FLORIS_CCBlade_control.add('ccblade', CCBladePlant())
    FLORIS_CCBlade_control.add('floris', florisPlant())

    FLORIS_CCBlade_control.driver.workflow.add(['windSpeedDistributor', 'control', 'ccblade', 'floris'])

    FLORIS_CCBlade_control.connect('control.listOfTurbinesOut','ccblade.listOfTurbinesIn')
    FLORIS_CCBlade_control.connect('ccblade.listOfTurbinesOut','floris.listOfTurbines')

    for turbineI in range(0, number_of_turbines):
        FLORIS_CCBlade_control.connect('windSpeedDistributor.wind_speed_out[%s]' % turbineI, 'control.listOfTurbinesIn[%s].wind_speed_hub' % turbineI)

    FLORIS_CCBlade_control.windSpeedDistributor.wind_speed_in = np.zeros(number_of_turbines)

    # connect FLORIS-predicted wind_speed_eff to windSpeedDistributor through FixedPointIterator
    FLORIS_CCBlade_control.driver.tolerance = tolerance
    FLORIS_CCBlade_control.driver.max_iteration = max_iteration
    for turbineI in range(0, number_of_turbines):
        FLORIS_CCBlade_control.driver.add_parameter('windSpeedDistributor.wind_speed_in[%s]' % turbineI, low = FLORIS_CCBlade_control.driver.tolerance, high = freestream_wind_speed, start = startpoint_wind_speed)
        FLORIS_CCBlade_control.driver.add_constraint("windSpeedDistributor.wind_speed_in[%s] = floris.velocitiesTurbines[%s]" % (turbineI, turbineI))

    return FLORIS_CCBlade_control