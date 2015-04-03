from openmdao.main.api import Assembly, Component
from openmdao.lib.datatypes.api import Array, Float, Bool, Int, List, Str, VarTree

from rotorse.rotoraerodefaults import CCBlade
from rotorse.rotoraero import Coefficients
from fusedwind.turbine.turbine_vt import AeroelasticHAWTVT
from openmdao.lib.drivers.api import BroydenSolver

import numpy as np


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


class AeroelasticHAWTVT_CCBlade(AeroelasticHAWTVT):
    """ AeroelasticHAWTVT with extra CCblade inputs and control inputs"""

    # inherits from AeroelasticHAWTVT:
    # turbine_name
    # rotor_diameter
    # nblades = 3 --> B
    # hub_height = 90.0 --> hubHt
    # tilt_angle --> tilt
    # cone_angle --> precone
    # hub_radius --> Rhub

    # flow
    rotor_area = Float(iotype='in', units='m*m', desc='rotor swept area')

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


class CCBladeAeroelasticHAWTVT(Assembly):
    """ CCBlade that takes AeroelasticHAWTVT_CCBlade as input
    """

    turbine = VarTree(AeroelasticHAWTVT_CCBlade(), iotype='in')
    power = Array(iotype='out', desc='rotor power', units='W')
    thrust = Array(iotype='out', desc='aerodynamic thrust on rotor', units='N')
    torque = Array(iotype='out', desc='aerodynamic torque on rotor', units='N*m')

    def configure(self):

        self.add("CCBlade", CCBlade())

        self.driver.workflow.add(['CCBlade'])

        self.CCBlade.run_case = 'power'

        self.connect("turbine.nblades", "CCBlade.B")
        self.connect("turbine.hub_height", "CCBlade.hubHt")
        self.connect("turbine.tilt_angle", "CCBlade.tilt")
        self.connect("turbine.cone_angle", "CCBlade.precone")
        self.connect("turbine.hub_radius", "CCBlade.Rhub")
        self.connect("turbine.tip_radius", "CCBlade.Rtip")
        self.connect("turbine.wind_speed_hub", "CCBlade.Uhub[0]")
        self.connect("turbine.air_density", "CCBlade.rho")
        self.connect("turbine.dynamic_viscosity", "CCBlade.mu")
        self.connect("turbine.shear_exponent", "CCBlade.shearExp")
        self.connect("turbine.radial_locations", "CCBlade.r")
        self.connect("turbine.chord", "CCBlade.chord")
        self.connect("turbine.theta", "CCBlade.theta")
        self.connect("turbine.precurve", "CCBlade.precurve")
        self.connect("turbine.precurveTip", "CCBlade.precurveTip")
        self.connect("turbine.airfoil_files", "CCBlade.airfoil_files")
        self.connect("turbine.nSector", "CCBlade.nSector")
        self.connect("turbine.tiploss", "CCBlade.tiploss")
        self.connect("turbine.hubloss", "CCBlade.hubloss")
        self.connect("turbine.wakerotation", "CCBlade.wakerotation")
        self.connect("turbine.usecd", "CCBlade.usecd")
        self.connect("turbine.yaw", "CCBlade.yaw")
        self.connect("turbine.pitch", "CCBlade.pitch[0]")
        self.connect("turbine.rotor_speed", "CCBlade.Omega[0]")

        self.connect("CCBlade.P", "power")
        self.connect("CCBlade.T", "thrust")
        self.connect("CCBlade.Q", "torque")

class CCBladeCoefficients(Assembly):
    """ Couples components with CCBlade so that it outputs CP and CT coefficients and axial-induction factor
    takes AeroelasticHAWTVT_CCBlade as input
    """

    turbine = VarTree(AeroelasticHAWTVT_CCBlade(), iotype='in')

    CP = Float(iotype='out', desc='power coefficient')
    CT = Float(iotype='out', desc='thrust coefficient')
    CQ = Float(iotype='out', desc='torque coefficient')
    axial_induction = Float(iotype='out', desc='axial induction')


    def configure(self):
        """ Creates a new Assembly containing a Paraboloid component"""

        self.add("CCBlade", CCBladeAeroelasticHAWTVT())
        self.add("Coefficients", Coefficients())
        self.add("CTtoAxialInd", CTtoAxialInd())

        self.driver.workflow.add(['CCBlade', 'Coefficients', 'CTtoAxialInd'])

        self.connect("turbine", "CCBlade.turbine")
        self.connect("CCBlade.power", "Coefficients.P")
        self.connect("CCBlade.thrust", "Coefficients.T")
        self.connect("CCBlade.torque", "Coefficients.Q")
        self.connect("turbine.wind_speed_hub", "Coefficients.V[0]")
        self.connect("turbine.tip_radius", "Coefficients.R")
        self.connect("turbine.air_density", "Coefficients.rho")
        self.connect("Coefficients.CT[0]", "CTtoAxialInd.CT")
        self.connect("Coefficients.CP[0]", "CP")
        self.connect("Coefficients.CT[0]", "CT")
        self.connect("Coefficients.CQ[0]", "CQ")
        self.connect("CTtoAxialInd.axial_induction", "axial_induction")


class CTtoAxialInd(Component):
    """Convert thrust coefficient to axial induction factor"""

    CT = Float(iotype='in')
    axial_induction = Float(iotype='out')

    def execute(self):
        self.axial_induction = 0.5*(1-np.sqrt(1-self.CT))


class AeroelasticHAWTVT_CCBlade_control(AeroelasticHAWTVT_CCBlade):
    opt_tsr = Float(iotype='in', desc='TSR to track in below-rated conditions')
    rated_generator_speed = Float(iotype='in', units='rpm', desc='rated generator speed')
    gearbox_ratio = Float(iotype='in', desc='gearbox ratio')
    rated_power = Float(iotype='in', units='W', desc='rated electric power')
    generator_efficiency = Float(iotype='in', desc='generator efficiency')


class windSpeedDistributor(Component):
    wind_speed = Float(iotype='in', units='m/s', desc='effective wind speed')
    wind_speed_controller = Float(iotype='out', units='m/s', desc='effective wind speed')
    wind_speed_turbine = Float(iotype='out', units='m/s', desc='effective wind speed')

    def execute(self):
        self.wind_speed_controller = self.wind_speed
        self.wind_speed_turbine = self.wind_speed


class Controller(Component):
    turbine = VarTree(AeroelasticHAWTVT_CCBlade_control(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of controller, False is no output')

    def execute(self):

        # set optimal rotor speed and zero pitch initially
        wind_speed_ax = np.cos(self.turbine.yaw*np.pi/180.0)*self.turbine.wind_speed_hub
        self.turbine.rotor_speed = wind_speed_ax*self.turbine.opt_tsr/self.turbine.tip_radius * 30.0/np.pi
        self.turbine.pitch = 0.0
        if self.verbose:
            if hasattr(self.turbine,'turbineName'):
                TopDelimiter = "_____Controller %s_______________________" % self.turbine.turbineName
            else:
                TopDelimiter = "_____Controller_______________________________"
            print TopDelimiter
            print "axial component wind speed %s" % wind_speed_ax

        # define constraints
        rated_rotor_speed = self.turbine.rated_generator_speed / self.turbine.gearbox_ratio
        rated_mechanical_power = self.turbine.rated_power/self.turbine.generator_efficiency
        rated_generator_torque = rated_mechanical_power/(self.turbine.rated_generator_speed*np.pi/30.0)
        rated_rotor_torque = rated_generator_torque*self.turbine.gearbox_ratio

        # if optimal rotor speed exceeds rated rotor speed, use rated rotor speed instead
        if self.turbine.rotor_speed > rated_rotor_speed:
            self.turbine.rotor_speed = rated_rotor_speed
            region = 2.5
        else:
            region = 2

        if self.verbose:
            print "rotor speed %s" % self.turbine.rotor_speed

        # calculate generator torque
        CCBlade = CCBladeAeroelasticHAWTVT()
        CCBlade.turbine = self.turbine
        CCBlade.turbine.wind_speed_hub = wind_speed_ax
        CCBlade.turbine.yaw = 0.0
        CCBlade.run()

        # if above-rated torque, calculate pitch_angle that will result in rated torque
        if CCBlade.torque>rated_rotor_torque:
            calcPitch = CalculateAboveRatedPitch()
            calcPitch.turbine = self.turbine
            calcPitch.rated_rotor_torque = rated_rotor_torque
            calcPitch.turbine.wind_speed_hub = wind_speed_ax
            calcPitch.turbine.yaw = 0.0
            calcPitch.run()
            self.turbine.pitch = calcPitch.CCBlade.turbine.pitch
            region = 3

        if self.verbose:
            print "pitch %s" % self.turbine.pitch
            print "region %s" % region
            print "-"*len(TopDelimiter)


class CalculateAboveRatedPitch(Assembly):
    turbine = VarTree(AeroelasticHAWTVT_CCBlade_control(), iotype='in')
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
            self.connect("turbine.%s" % attribute, "CCBlade.turbine.%s" % attribute)

        self.optimizer.itmax = 20
        self.optimizer.tol = 0.000001
        self.optimizer.workflow.add(['CCBlade'])
        self.optimizer.add_parameter('CCBlade.turbine.pitch', low=0.0, high=90.0, start=10.0)
        self.optimizer.add_constraint('CCBlade.CQ = CQfromTorque.CQ')


class CQfromTorque(Component):
    torque = Float(iotype='in', units='N*m', desc='rotor torque')
    air_density = Float(iotype='in', units='kg/m**3', desc='density of air')
    rotor_area = Float(iotype='in', units='m*m', desc='rotor swept area')
    tip_radius = Float(iotype='in', units='m', desc='tip radius')
    wind_speed = Float(iotype='in', units='m/s', desc='effective wind speed')
    CQ = Float(iotype='out', desc='torque coefficient')

    def execute(self):
        self.CQ = self.torque / (0.50 * self.air_density * self.rotor_area * self.tip_radius * (self.wind_speed ** 2))