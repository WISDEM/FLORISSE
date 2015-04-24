from openmdao.main.api import Assembly, Component, VariableTree
from openmdao.lib.datatypes.api import Array, Float, Bool, Int, List, Str, VarTree

from rotorse.rotoraerodefaults import CCBlade
from rotorse.rotoraero import Coefficients
from floris import FLORIS
from fusedwind.turbine.turbine_vt import AeroelasticHAWTVT
from fusedwind.plant_flow.comp import GenericWindFarm
from openmdao.lib.drivers.api import BroydenSolver, CaseIteratorDriver, FixedPointIterator
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout
from floris import FLORIS
from scipy import interp
from copy import copy

import os
import numpy as np
import matplotlib.pyplot as plt


class PowerSpeedControllerPreCalculated(VariableTree):
    wind_speeds = Array(iotype='in', units='m/s', desc='range of wind speeds')
    pitch = Array(iotype='out', units='deg', desc='pitch control setting')
    rotor_speed = Array(iotype='out', units='rpm', desc='rotor speed control setting')


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

    # information required for control
    opt_tsr = Float(iotype='in', desc='TSR to track in below-rated conditions')
    rated_generator_speed = Float(iotype='in', units='rpm', desc='rated generator speed')
    cut_in_wind_speed = Float(iotype='in', units='m/s', desc='cut-in wind speed')
    cut_out_wind_speed = Float(iotype='in', units='m/s', desc='cut-out wind speed')
    minimum_generator_speed = Float(iotype='in', units='rpm', desc='minimum generator speed')
    transitional_generator_speed = Float(iotype='in', desc='transitional generator speed between regions 1 and 1.5')

    gearbox_ratio = Float(iotype='in', desc='gearbox ratio')
    rated_power = Float(iotype='in', units='W', desc='rated electric power')
    generator_efficiency = Float(iotype='in', desc='generator efficiency')
    preCalculated = Bool(False, iotype='in', desc='choice on whether or not to use pre-calculated controller')
    PreCalculatedPowerSpeedController = VarTree(PowerSpeedControllerPreCalculated(), iotype='in', desc='pre-calculated control schedule')

    # FLORIS outputs
    wind_speed_eff = Float(iotype='out', desc='effective wind speed', units='m/s')

    def preCalculateController(self, wind_speeds=None, wind_speed_resolution=200, visual=False):

        # calculate default range of wind speeds
        if wind_speeds is None:

            # calculate rated wind speed
            calcRatedWindSpeed = CalculateRatedWindSpeed()
            calcRatedWindSpeed.turbine = self
            calcRatedWindSpeed.run()

            wind_speeds = np.hstack((self.cut_in_wind_speed-0.01, np.linspace(self.cut_in_wind_speed, calcRatedWindSpeed.rated_wind_speed, np.round(wind_speed_resolution/3.0)), np.logspace(np.log10(calcRatedWindSpeed.rated_wind_speed+0.001), np.log10(self.cut_out_wind_speed), np.round(2.0*wind_speed_resolution/3.0)), self.cut_out_wind_speed+0.001))

        precalc = preCalculatePowerAndSpeedController()
        precalc.wind_speeds = wind_speeds
        precalc.turbine = self
        precalc.run()
        self.PreCalculatedPowerSpeedController = precalc.PreCalculatedPowerSpeedController
        self.preCalculated = True

        if visual:
            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax1.plot(self.PreCalculatedPowerSpeedController.wind_speeds, self.PreCalculatedPowerSpeedController.pitch,'.-')
            ax1.set_xlabel('wind speed (m/s)')
            ax1.set_ylabel('pitch (deg)')
            ax2 = fig.add_subplot(212)
            ax2.plot(self.PreCalculatedPowerSpeedController.wind_speeds, self.PreCalculatedPowerSpeedController.rotor_speed,'.-')
            ax2.set_xlabel('wind speed (m/s)')
            ax2.set_ylabel('rotor speed (rpm)')
            plt.show()


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
        if self.CT > 0.96: # Glauert condition
            self.axial_induction = 0.143+np.sqrt(0.0203-0.6427*(0.889-self.CT))
        else:
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
    wind_speed = Float(iotype='in')
    wind_direction = Float(iotype='in')
    air_density = Float(iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of floris, False is no output')

    velocitiesTurbines = Array(iotype='out')
    power = Float(iotype='out')
    wt_power = Array(iotype='out')

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
        self.florisWindPlant.verbose = self.verbose

        self.florisWindPlant.run()

        self.velocitiesTurbines = self.florisWindPlant.velocitiesTurbines
        self.power = self.florisWindPlant.power
        self.wt_power = self.florisWindPlant.wt_power


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
    PowerSpeedControllerPreCalculated = VarTree(PowerSpeedControllerPreCalculated(), iotype='in', desc='pre-calculated pitch and torque control settings for range of wind speeds')

    def execute(self):
        turbine = self.turbineIn

        # calculate axial wind speed
        wind_speed_ax = np.cos(turbine.yaw*np.pi/180.0)*turbine.wind_speed_hub

        if self.verbose:
            if hasattr(turbine,'turbineName'):
                TopDelimiter = "_____Controller %s_______________________" % turbine.turbineName
            else:
                TopDelimiter = "_____Controller_______________________________"
            print TopDelimiter
            print "axial component wind speed %s" % wind_speed_ax

        if self.turbineIn.preCalculated:
            # interpolate pre-calculated pitch and torque policy
            if self.verbose:
                print 'interpolates precalculated controller'
            turbine.pitch = interp(wind_speed_ax,turbine.PreCalculatedPowerSpeedController.wind_speeds, turbine.PreCalculatedPowerSpeedController.pitch)
            turbine.rotor_speed = interp(wind_speed_ax,turbine.PreCalculatedPowerSpeedController.wind_speeds, turbine.PreCalculatedPowerSpeedController.rotor_speed)
        else:
            # control not pre-calculated, find out which region and calculate pitch and rotor speed

            # define constraints
            rated_rotor_speed = turbine.rated_generator_speed / turbine.gearbox_ratio
            transitional_rotor_speed = turbine.transitional_generator_speed / turbine.gearbox_ratio
            rated_mechanical_power = turbine.rated_power/turbine.generator_efficiency
            rated_generator_torque = rated_mechanical_power/(turbine.rated_generator_speed*np.pi/30.0)
            rated_rotor_torque = rated_generator_torque*turbine.gearbox_ratio

            if wind_speed_ax > turbine.cut_out_wind_speed:
                region = 4

                # in region 4, use zero rotor speed and pitch angle for cut-out wind speed
                calcPitch = CalculateAboveRatedPitch()
                calcPitch.turbine = copy(turbine)
                calcPitch.turbine.rotor_speed = rated_rotor_speed
                calcPitch.rated_rotor_torque = rated_rotor_torque
                calcPitch.turbine.wind_speed_hub = turbine.cut_out_wind_speed
                calcPitch.turbine.yaw = 0.0
                calcPitch.run()
                calcPitch.turbine.rotor_speed = rated_rotor_speed

                turbine.pitch = calcPitch.CCBlade.turbineIn.pitch
                turbine.rotor_speed = 0.0

            elif wind_speed_ax < turbine.cut_in_wind_speed:
                region = 1

                # in region 1, use zero rotor speed and pitch angle
                turbine.rotor_speed = 0.0
                turbine.pitch = 0.0

            else:
                # set optimal rotor speed and zero pitch initially (region 2), then check if not actually in
                # region 3 or 1.5
                region = 2
                turbine.rotor_speed = wind_speed_ax*turbine.opt_tsr/turbine.tip_radius * 30.0/np.pi
                turbine.pitch = 0.0

                if turbine.rotor_speed > rated_rotor_speed:
                    region = 2.5

                    # if optimal rotor speed exceeds rated rotor speed, use rated rotor speed instead
                    turbine.rotor_speed = rated_rotor_speed

                    # calculate torque to determine if in region 3
                    CCBlade = CCBladeCoefficients()
                    CCBlade.turbineIn = copy(turbine)
                    CCBlade.turbineIn.wind_speed_hub = wind_speed_ax
                    CCBlade.turbineIn.yaw = 0.0
                    CCBlade.run()

                    if CCBlade.turbineOut.torque>rated_rotor_torque:
                        region = 3

                        # if above-rated, calculate pitch_angle that will result in rated torque
                        calcPitch = CalculateAboveRatedPitch()
                        calcPitch.turbine = copy(turbine)
                        calcPitch.turbine.rotor_speed = rated_rotor_speed
                        calcPitch.rated_rotor_torque = rated_rotor_torque
                        calcPitch.turbine.wind_speed_hub = wind_speed_ax
                        calcPitch.turbine.yaw = 0.0
                        calcPitch.run()
                        turbine.pitch = calcPitch.CCBlade.turbineIn.pitch
                        turbine.rotor_speed = rated_rotor_speed
                elif turbine.rotor_speed < transitional_rotor_speed:

                    # if in region 1.5, use torque balance to calculate rotor speed
                    region = 1.5

                    calcRotorSpeed = calculateRotorSpeed15()
                    calcRotorSpeed.turbine = copy(turbine)
                    calcRotorSpeed.run()
                    turbine.pitch = 0.0
                    turbine.rotor_speed = calcRotorSpeed.CCBlade.turbineIn.rotor_speed
            if self.verbose:
                print "region %s" % region

        if self.verbose:
            print "rotor speed %s" % turbine.rotor_speed
            print "pitch %s" % turbine.pitch

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


class CalculateRatedWindSpeed(Assembly):
    turbine = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='in')
    rated_wind_speed = Float(iotype='out', units='m/s', desc='rated wind speed (with zero yaw)')

    def configure(self):

        # first calculate rated rotor speed from rated generator speed
        self.add('calculateRatedRotorSpeed', calculateRatedRotorSpeed())
        self.add('calculateRatedMechanicalPower', calculateRatedMechanicalPower())
        self.connect("turbine.rated_generator_speed", "calculateRatedRotorSpeed.rated_generator_speed")
        self.connect("turbine.gearbox_ratio", "calculateRatedRotorSpeed.gearbox_ratio")
        self.connect("turbine.rated_power", "calculateRatedMechanicalPower.rated_power")
        self.connect("turbine.generator_efficiency", "calculateRatedMechanicalPower.generator_efficiency")

        self.add("CCBlade", CCBladeCoefficients())
        self.add('optimizer', BroydenSolver())

        self.driver.workflow.add(['calculateRatedRotorSpeed', 'calculateRatedMechanicalPower', 'optimizer'])

        # define all CCBlade turbine properties except for wind_speed_hub, for which we use optimizer
        attributes = ("nblades", "hub_height", "tilt_angle", "cone_angle", "hub_radius", "tip_radius","generator_efficiency",
                      "air_density", "dynamic_viscosity", "shear_exponent", "radial_locations", "chord", "theta", "precurve",
                      "precurveTip", "airfoil_files", "nSector", "tiploss", "hubloss", "wakerotation", "usecd")
        for attribute in attributes:
            self.connect("turbine.%s" % attribute, "CCBlade.turbineIn.%s" % attribute)
        # set rotor speed to rated rotor speed
        self.connect("calculateRatedRotorSpeed.rated_rotor_speed", "CCBlade.turbineIn.rotor_speed")
        # set pitch and yaw to zero
        self.CCBlade.turbineIn.pitch = 0.0
        self.CCBlade.turbineIn.yaw = 0.0

        # let CCBlade calculate wind speed for which turbine will hit rated power
        self.optimizer.itmax = 20
        self.optimizer.tol = 0.000001
        self.optimizer.workflow.add(['CCBlade'])
        self.optimizer.add_parameter('CCBlade.turbineIn.wind_speed_hub', low=0.0, high=40.0, start=10.0)
        self.optimizer.add_constraint('CCBlade.turbineOut.power = calculateRatedMechanicalPower.rated_mechanical_power')
        self.connect('CCBlade.turbineIn.wind_speed_hub', 'rated_wind_speed')


class calculateRatedRotorSpeed(Component):
    rated_generator_speed = Float(iotype='in', units='rpm', desc='rated generator speed')
    gearbox_ratio = Float(iotype='in', desc='gearbox ratio')
    rated_rotor_speed = Float(iotype='out', units='rpm', desc='rated rotor speed')

    def execute(self):
        self.rated_rotor_speed = self.rated_generator_speed / self.gearbox_ratio

class calculateRatedMechanicalPower(Component):
    rated_power = Float(iotype='in', units='W', desc='rated generator power')
    generator_efficiency = Float(iotype='in', desc='generator efficiency')
    rated_mechanical_power = Float(iotype='out', units='W', desc='rated mechanical power')

    def execute(self):
        self.rated_mechanical_power = self.rated_power / self.generator_efficiency


class calculateRotorSpeed15(Assembly):
    # calculates region 1.5 rotor speed for a certain wind speed, by balancing aerodynamic torque with mechanical torque
    turbine = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='in')

    def configure(self):
        self.add('calculateMaximumTorqueRegion15',calculateMaximumTorqueRegion15())
        self.add('optimizer',BroydenSolver())

        self.driver.workflow.add(['calculateMaximumTorqueRegion15', 'optimizer'])

        self.connect('turbine', 'calculateMaximumTorqueRegion15.turbine')

        self.add("CCBlade", CCBladeCoefficients())
        self.add("calculateTorqueRegion15", calculateTorqueRegion15())

        attributes = ("nblades", "hub_height", "tilt_angle", "cone_angle", "hub_radius", "tip_radius",
              "air_density", "dynamic_viscosity", "shear_exponent", "radial_locations", "chord", "theta", "precurve",
              "precurveTip", "airfoil_files", "nSector", "tiploss", "hubloss", "wakerotation", "usecd", "wind_speed_hub")
        for attribute in attributes:
            self.connect("turbine.%s" % attribute, "CCBlade.turbineIn.%s" % attribute)
            self.connect("turbine.%s" % attribute, "calculateTorqueRegion15.turbine.%s" % attribute)

        attributes = ("minimum_generator_speed", "transitional_generator_speed", "gearbox_ratio")
        for attribute in attributes:
            self.connect("turbine.%s" % attribute, "calculateTorqueRegion15.turbine.%s" % attribute)

        self.optimizer.itmax = 20
        self.optimizer.tol = 0.000001
        self.optimizer.workflow.add(['CCBlade','calculateTorqueRegion15'])
        self.optimizer.add_parameter('CCBlade.turbineIn.rotor_speed', low=0.01, high=50.0)
        self.connect('CCBlade.turbineIn.rotor_speed', 'calculateTorqueRegion15.turbine.rotor_speed')
        self.connect('calculateMaximumTorqueRegion15.maximum_generator_torque_region_15','calculateTorqueRegion15.maximum_generator_torque_region_15')
        self.optimizer.add_constraint('CCBlade.turbineOut.torque = calculateTorqueRegion15.mechanical_torque_region_15')


class calculateTorqueRegion15(Component):
    # calculates region 1.5 mechanical torque based on rotor speed
    # (linear transition between zero torque at minimum_generator_speed and torque = maximum_generator_torque_region_15 at
    # transitional_generator_speed)
    turbine = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='in')
    maximum_generator_torque_region_15 = Float(iotype='in', desc='maximum generator torque region 1.5', units='N*m')
    mechanical_torque_region_15 = Float(iotype='out', desc='mechanical torque in region 1.5', units='N*m')

    def execute(self):
        generator_speed = self.turbine.rotor_speed*self.turbine.gearbox_ratio
        generator_torque = ((generator_speed-self.turbine.minimum_generator_speed)/(self.turbine.transitional_generator_speed-self.turbine.minimum_generator_speed))*self.maximum_generator_torque_region_15
        self.mechanical_torque_region_15 = generator_torque*self.turbine.gearbox_ratio


class calculateMaximumTorqueRegion15(Assembly):
    turbine = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='in')
    maximum_generator_torque_region_15 = Float(iotype='out', desc='maximum generator torque region 1.5', units='N*m')

    def execute(self):
        CCBlade = CCBladeCoefficients()
        CCBlade.turbineIn = self.turbine
        CCBlade.turbineIn.wind_speed_hub = (((self.turbine.transitional_generator_speed*np.pi/30.0)/self.turbine.gearbox_ratio)
                                            *self.turbine.tip_radius)/self.turbine.opt_tsr  # set wind speed to transitional wind speed
        CCBlade.turbineIn.rotor_speed = self.turbine.transitional_generator_speed/self.turbine.gearbox_ratio
        CCBlade.turbineIn.pitch = 0.0
        CCBlade.turbineIn.yaw = 0.0
        CCBlade.run()
        self.maximum_generator_torque_region_15 = CCBlade.turbineOut.torque/self.turbine.gearbox_ratio


class preCalculatePowerAndSpeedController(Assembly):
    wind_speeds = Array(iotype='in', units='m/s', desc='effective axial wind speeds for which to pre-calculate pitch and torque control settings')
    turbine = VarTree(AeroelasticHAWTVT_CCBlade_floris(), iotype='in')
    PreCalculatedPowerSpeedController = VarTree(PowerSpeedControllerPreCalculated(), iotype='out')

    def configure(self):

        self.add("PowerSpeedController", PowerSpeedController())
        self.add('driver', CaseIteratorDriver())
        # connect all needed CCBlade and control turbine properties
        attributes = ("nblades", "hub_height", "tilt_angle", "cone_angle", "hub_radius", "tip_radius",
                      "air_density","dynamic_viscosity","shear_exponent","radial_locations","chord","theta","precurve",
                      "precurveTip","airfoil_files","nSector","tiploss","hubloss","wakerotation","usecd","rotor_speed",
                      "turbineName", "rated_generator_speed", "gearbox_ratio", "rated_power", "generator_efficiency",
                      "opt_tsr", "cut_in_wind_speed", "cut_out_wind_speed", "rotor_area","minimum_generator_speed",
                      "transitional_generator_speed")
        for attribute in attributes:
            self.connect("turbine.%s" % attribute, "PowerSpeedController.turbineIn.%s" % attribute)

        # set controller options
        self.PowerSpeedController.turbineIn.yaw = 0.0
        self.PowerSpeedController.preCalculated = False
        self.PowerSpeedController.verbose = True

        # set up CaseIterator to calculate pitch and rotor_speed for range of wind speeds
        self.driver.add_parameter('PowerSpeedController.turbineIn.wind_speed_hub')
        self.connect('wind_speeds', 'driver.case_inputs.PowerSpeedController.turbineIn.wind_speed_hub')
        self.driver.add_response('PowerSpeedController.turbineOut.pitch')
        self.driver.add_response('PowerSpeedController.turbineOut.rotor_speed')
        self.connect('wind_speeds', 'PreCalculatedPowerSpeedController.wind_speeds')
        self.connect('driver.case_outputs.PowerSpeedController.turbineOut.pitch', 'PreCalculatedPowerSpeedController.pitch')
        self.connect('driver.case_outputs.PowerSpeedController.turbineOut.rotor_speed', 'PreCalculatedPowerSpeedController.rotor_speed')


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

    FLORIS_CCBlade_control.create_passthrough('floris.wind_speed', 'wind_speed')
    FLORIS_CCBlade_control.create_passthrough('floris.wind_direction', 'wind_direction')
    FLORIS_CCBlade_control.create_passthrough('floris.air_density', 'air_density')
    FLORIS_CCBlade_control.create_passthrough('control.verbose', 'controller_verbosity')
    FLORIS_CCBlade_control.create_passthrough('floris.verbose', 'floris_verbosity')
    FLORIS_CCBlade_control.create_passthrough('floris.power', 'power')
    FLORIS_CCBlade_control.create_passthrough('floris.wt_power', 'wt_power')

    return FLORIS_CCBlade_control



class CoupledFLORIS_CCBlade_control(GenericWindFarm):

    air_density = Float(iotype='in', units='kg/m**3', desc='density of air')
    controller_verbosity = Bool(False, iotype='in', desc='verbosity of controller, False is no output')
    floris_verbosity = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')
    coupledModel_verbosity = Bool(False, iotype='in', desc='verbosity of overall model, False is no output')

    def execute(self):
        if self.coupledModel_verbosity:
            print " construct coupled model..."

        self.coupledModel = constructCoupledFLORIS_CCBlade_control(len(self.wt_layout.wt_list), startpoint_wind_speed = self.wind_speed)
        self.coupledModel.control.listOfTurbinesIn = self.wt_layout.wt_list
        self.coupledModel.wind_speed = self.wind_speed
        self.coupledModel.air_density = self.air_density
        self.coupledModel.wind_direction = self.wind_direction
        self.coupledModel.controller_verbosity = self.controller_verbosity
        self.coupledModel.floris_verbosity = self.floris_verbosity

        if self.coupledModel_verbosity:
            print " ... done"

        if self.coupledModel_verbosity:
            np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print " ======== COUPLED MODEL: FLORIS-SE+ROTOR-SE+CONTROL ============="
            print " free-stream wind direction %0.3f deg" % self.wind_direction
            print " free-stream wind speed %0.3f m/s" % self.wind_speed

        self.coupledModel.run()
        self.power = self.coupledModel.power
        self.wt_power = self.coupledModel.wt_power

        if self.coupledModel_verbosity:
            print " no. of iterations in model evaluation: %s" % self.coupledModel.floris.exec_count

        if self.coupledModel_verbosity:
            print_power = self.wt_power/1000.0
            if len(self.wt_power) > 6:
                print " turbine powers: [ %0.3f  %0.3f  %0.3f  %0.3f  %0.3f ... ] MW" % (print_power[0], print_power[1], print_power[2], print_power[3], print_power[4])
            else:
                print " turbine powers: %s MW" % print_power
            print " total wind plant power: %0.3f MW" % (self.power/1000.0)
            print " ================================================================"








