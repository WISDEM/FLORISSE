from openmdao.main.api import Component, VariableTree
#from ccblade import CCAirfoil, CCBlade
from openmdao.lib.datatypes.api import Array, Float, Bool, Int, List, Str, VarTree
import numpy as np
from scipy import interp

class PowerSpeedControllerPreCalculated(VariableTree):

    def __init__(self, nTurbines):

        super(PowerSpeedControllerPreCalculated, self).__init__()

        self.add('wind_speeds', Array(np.array(nTurbines), iotype='in', units='m/s', desc='range of wind speeds'))
        self.add('pitch', Array(np.array(nTurbines), iotype='in', units='deg', desc='pitch control setting'))
        self.add('rotor_speed', Array(np.array(nTurbines), iotype='in', units='rpm', desc='rotor speed control setting'))

class turbine_ccBlade(VariableTree):
    nblades = Int(3, desc='number of blades')
    rotor_diameter = Float(units='m', desc='Rotor Diameter')
    hub_height = Float(units='m', desc='Hub height')
    tilt_angle = Float(units='deg', desc='Rotor tilt angle')
    cone_angle = Float(units='deg', desc='Rotor cone angle')
    hub_radius = Float(units='m', desc='Hub radius')

    # CCblade inputs
    tip_radius = Float(iotype='in', units='m', desc='tip radius') # Rtip

    precurveTip = Float(0.0, iotype='in', units='m', desc='precurve at tip')
    airfoil_files = List(Str, iotype='in', desc='names of airfoil file')
    nSector = Int(4, iotype='in', desc='number of sectors to divide rotor face into in computing thrust and power')

    #tiploss = Bool(True, iotype='in', desc='include Prandtl tip loss model')
    #hubloss = Bool(True, iotype='in', desc='include Prandtl hub loss model')
    #wakerotation = Bool(True, iotype='in', desc='include effect of wake rotation (i.e., tangential induction factor is nonzero)')
    #usecd = Bool(True, iotype='in', desc='use drag coefficient in computing induction factors')

    def __init__(self, nTurbines):

        super(turbine_ccBlade, self).__init__()

        self.add('radial_locations', Array(np.zeros(nTurbines), iotype='in', units='m', desc='radial locations where blade is defined (should be increasing and not go all the way to hub or tip)')) #r
        self.add('chord', Array(np.zeros(nTurbines), iotype='in', units='m', desc='chord length at each section'))
        self.add('theta', Array(np.zeros(nTurbines), iotype='in', units='deg', desc='twist angle at each section (positive decreases angle of attack)'))
        self.add('precurve', Array(np.zeros(nTurbines), iotype='in', units='m', desc='precurve at each section'))


class ccblade_CPCT(Component):
    # flow
    air_density = Float(1.225, iotype='in', units='kg/m**3', desc='density of air') # rho
    dynamic_viscosity = Float(1.81206e-5, iotype='in', units='kg/(m*s)', desc='dynamic viscosity of air') # mu
    shear_exponent = Float(0.2, iotype='in', desc='shear exponent') # shearExp

    def __init__(self, nTurbines):

        super(ccblade_CPCT, self).__init__()

        self.nTurbines = nTurbines

        self.add('turbine', VarTree(turbine_ccBlade(nTurbines), iotype='in'))

        self.add('wind_speed_hub', Array(np.zeros(nTurbines), iotype='in', units='m/s', desc='hub height wind speed')) # UHub

        # control DOFs
        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', desc='yaw error', units='deg'))
        self.add('pitch', Array(np.zeros(nTurbines), iotype='in', desc='blade pitch angle', units='deg'))
        self.add('rotor_speed', Array(np.zeros(nTurbines), iotype='in', desc='rotor speed', units='rpm')) #Omega

        # coefficients
        self.add('CP', Array(np.zeros(nTurbines), iotype='out', desc='power coefficient'))
        self.add('CT', Array(np.zeros(nTurbines), iotype='out', desc='thrust coefficient'))

    def execute(self):
        turbine = self.turbine
        self.CP = np.zeros_like(self.yaw)
        self.CT = np.zeros_like(self.yaw)
        for turbI, yaw in enumerate(self.yaw):

            afinit = CCAirfoil.initFromAerodynFile
            airfoil = [afinit(f) for f in turbine.airfoil_files]
            rotor = CCBlade(turbine.radial_locations, turbine.chord, turbine.theta, airfoil, turbine.hub_radius, turbine.tip_radius, turbine.nblades, self.air_density, self.dynamic_viscosity, turbine.cone_angle, turbine.tilt_angle, yaw, self.shear_exponent, turbine.hub_height, turbine.nSector)
            CP, CT, CQ = rotor.evaluate(np.array([self.wind_speed_hub[turbI]]), np.array([self.rotor_speed[turbI]]), np.array([self.pitch[turbI]]), coefficient=True)
            self.CP[turbI] = CP
            self.CT[turbI] = CT
        print 'CP %s' % self.CP
        print 'CT %s' % self.CT

class power_and_speed_controller_Interpolate(Component):

    def __init__(self, nTurbines):

        super(power_and_speed_controller_Interpolate, self).__init__()

        self.nTurbines = nTurbines

        self.add('PreCalculatedPowerSpeedController', VarTree(PowerSpeedControllerPreCalculated(nTurbines), iotype='in',
                                                              desc='pre-calculated control schedule'))

        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', desc='yaw error', units='deg'))
        self.add('wind_speed_hub', Array(np.zeros(nTurbines), iotype='in', units='m/s', desc='hub height wind speed')) # Uhub
        self.add('pitch', Array(np.zeros(nTurbines), iotype='out', desc='blade pitch angle', units='deg'))
        self.add('rotor_speed', Array(np.zeros(nTurbines), iotype='out', desc='rotor speed', units='rpm')) #Omega

    def execute(self):
        wind_speed_ax = np.cos(self.yaw*np.pi/180.0)*self.wind_speed_hub
        print '>>>>>>>  wind speed ax %s' % wind_speed_ax
        self.pitch = interp(wind_speed_ax, self.PreCalculatedPowerSpeedController.wind_speeds, self.PreCalculatedPowerSpeedController.pitch)
        self.rotor_speed = interp(wind_speed_ax, self.PreCalculatedPowerSpeedController.wind_speeds, self.PreCalculatedPowerSpeedController.rotor_speed)
        print 'pitch %s' % self.pitch
        print 'rotor speed %s' % self.rotor_speed

## ---- if you know wind speed to power and thrust, you can use these tools ----------------


class windSpeedToCPCT(VariableTree):
    # variable tree that defines known wind speed to CP/CT curve
    def __init__(self, datasize=0):

        super(windSpeedToCPCT, self).__init__()

        self.add('wind_speed', Array(np.zeros(datasize), iotype='in', units='m/s', desc='range of wind speeds'))
        self.add('CP', Array(np.zeros(datasize), iotype='out', desc='power coefficients'))
        self.add('CT', Array(np.zeros(datasize), iotype='out', desc='thrust coefficients'))


# class CPCT_Interpolate(Component):

#     pP = Float(3.0, iotype='in')

#     def __init__(self, nTurbines, datasize=0):

#         super(CPCT_Interpolate, self).__init__()

#         self.nTurbines = nTurbines
#         self.datasize = datasize
#         self.add('windSpeedToCPCT', VarTree(windSpeedToCPCT(datasize), iotype='in', desc='pre-calculated CPCT'))

#         self.add('yaw', Array(np.zeros(nTurbines), iotype='in', desc='yaw error', units='deg'))
#         self.add('wind_speed_hub', Array(np.zeros(nTurbines), iotype='in', units='m/s', desc='hub height wind speed')) # Uhub
#         self.add('CP', Array(np.zeros(nTurbines), iotype='out'))
#         self.add('CT', Array(np.zeros(nTurbines), iotype='out'))

#     def execute(self):
#         wind_speed_ax = np.cos(self.yaw*np.pi/180.0)**(self.pP/3.0)*self.wind_speed_hub
#         # use interpolation on precalculated CP-CT curve
#         wind_speed_ax = np.maximum(wind_speed_ax, self.windSpeedToCPCT.wind_speed[0])
#         wind_speed_ax = np.minimum(wind_speed_ax, self.windSpeedToCPCT.wind_speed[-1])
#         self.CP = interp(wind_speed_ax, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CP)
#         self.CT = interp(wind_speed_ax, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CT)

#         # normalize on incoming wind speed to correct coefficients for yaw
#         self.CP = self.CP * np.cos(self.yaw*np.pi/180.0)**self.pP
#         self.CT = self.CT * np.cos(self.yaw*np.pi/180.0)**2
#         # print 'in CPCT interp, wind_speed_hub = ', self.wind_speed_hub

class CPCT_Interpolate(Component):

    pP = Float(3.0, iotype='in')

    def __init__(self, nTurbines, datasize=0):

        #print 'in CPCT_Interpolate init'

        super(CPCT_Interpolate, self).__init__()

        self.nTurbines = nTurbines
        self.datasize = datasize
        self.add('windSpeedToCPCT', VarTree(windSpeedToCPCT(datasize), iotype='in', desc='pre-calculated CPCT'))

        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', desc='yaw error', units='deg'))
        self.add('wind_speed_hub', Array(np.zeros(nTurbines), iotype='in', units='m/s', desc='hub height wind speed')) # Uhub
        self.add('CP', Array(np.zeros(nTurbines), iotype='out'))
        self.add('CT', Array(np.zeros(nTurbines), iotype='out'))

    def execute(self):

        #print 'in CPCT_Interpolate'

        wind_speed_ax = np.cos(self.yaw*np.pi/180.0)**(self.pP/3.0)*self.wind_speed_hub
        # use interpolation on precalculated CP-CT curve
        wind_speed_ax = np.maximum(wind_speed_ax, self.windSpeedToCPCT.wind_speed[0])
        wind_speed_ax = np.minimum(wind_speed_ax, self.windSpeedToCPCT.wind_speed[-1])
        self.CP = interp(wind_speed_ax, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CP)
        self.CT = interp(wind_speed_ax, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CT)

        # normalize on incoming wind speed to correct coefficients for yaw
        self.CP = self.CP * np.cos(self.yaw*np.pi/180.0)**self.pP
        self.CT = self.CT * np.cos(self.yaw*np.pi/180.0)**2
        # print 'in CPCT interp, wind_speed_hub = ', self.wind_speed_hub

    def list_deriv_vars(self):
        return ('yaw', 'wind_speed_hub'), ('CP', 'CT')

    def provideJ(self):

        #print 'in CPCT_Interpolate - provideJ'
        
         # standard central differencing
        # set step size for finite differencing
        h = 1e-6

        # calculate upper and lower function values
        wind_speed_ax_high_yaw = np.cos((self.yaw+h)*np.pi/180.0)**(self.pP/3.0)*self.wind_speed_hub
        wind_speed_ax_low_yaw = np.cos((self.yaw-h)*np.pi/180.0)**(self.pP/3.0)*self.wind_speed_hub
        wind_speed_ax_high_wind = np.cos(self.yaw*np.pi/180.0)**(self.pP/3.0)*(self.wind_speed_hub+h)
        wind_speed_ax_low_wind = np.cos(self.yaw*np.pi/180.0)**(self.pP/3.0)*(self.wind_speed_hub-h)

        # use interpolation on precalculated CP-CT curve
        wind_speed_ax_high_yaw = np.maximum(wind_speed_ax_high_yaw, self.windSpeedToCPCT.wind_speed[0])
        wind_speed_ax_low_yaw = np.maximum(wind_speed_ax_low_yaw, self.windSpeedToCPCT.wind_speed[0])
        wind_speed_ax_high_wind = np.maximum(wind_speed_ax_high_wind, self.windSpeedToCPCT.wind_speed[0])
        wind_speed_ax_low_wind = np.maximum(wind_speed_ax_low_wind, self.windSpeedToCPCT.wind_speed[0])

        wind_speed_ax_high_yaw = np.minimum(wind_speed_ax_high_yaw, self.windSpeedToCPCT.wind_speed[-1])
        wind_speed_ax_low_yaw = np.minimum(wind_speed_ax_low_yaw, self.windSpeedToCPCT.wind_speed[-1])
        wind_speed_ax_high_wind = np.minimum(wind_speed_ax_high_wind, self.windSpeedToCPCT.wind_speed[-1])
        wind_speed_ax_low_wind = np.minimum(wind_speed_ax_low_wind, self.windSpeedToCPCT.wind_speed[-1])

        CP_high_yaw = interp(wind_speed_ax_high_yaw, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CP)
        CP_low_yaw = interp(wind_speed_ax_low_yaw, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CP)
        CP_high_wind = interp(wind_speed_ax_high_wind, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CP)
        CP_low_wind = interp(wind_speed_ax_low_wind, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CP)

        CT_high_yaw = interp(wind_speed_ax_high_yaw, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CT)
        CT_low_yaw = interp(wind_speed_ax_low_yaw, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CT)
        CT_high_wind = interp(wind_speed_ax_high_wind, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CT)
        CT_low_wind = interp(wind_speed_ax_low_wind, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CT)

        # normalize on incoming wind speed to correct coefficients for yaw
        CP_high_yaw = CP_high_yaw * np.cos((self.yaw+h)*np.pi/180.0)**self.pP
        CP_low_yaw = CP_low_yaw * np.cos((self.yaw-h)*np.pi/180.0)**self.pP
        CP_high_wind = CP_high_wind * np.cos((self.yaw)*np.pi/180.0)**self.pP
        CP_low_wind = CP_low_wind * np.cos((self.yaw)*np.pi/180.0)**self.pP

        CT_high_yaw = CT_high_yaw * np.cos((self.yaw+h)*np.pi/180.0)**2
        CT_low_yaw = CT_low_yaw * np.cos((self.yaw-h)*np.pi/180.0)**2
        CT_high_wind = CT_high_wind * np.cos((self.yaw)*np.pi/180.0)**2
        CT_low_wind = CT_low_wind * np.cos((self.yaw)*np.pi/180.0)**2

        # compute derivative via central differencing and arrange in sub-matrices of the Jacobian
        dCP_dyaw = np.eye(self.nTurbines)*(CP_high_yaw-CP_low_yaw)/(2.0*h)
        dCP_dwind = np.eye(self.nTurbines)*(CP_high_wind-CP_low_wind)/(2.0*h)
        dCT_dyaw = np.eye(self.nTurbines)*(CT_high_yaw-CT_low_yaw)/(2.0*h)
        dCT_dwind = np.eye(self.nTurbines)*(CT_high_wind-CT_low_wind)/(2.0*h)

        # compile full Jacobian from sub-matrices
        dCP = np.hstack((dCP_dyaw, dCP_dwind))
        dCT = np.hstack((dCT_dyaw, dCT_dwind))
        J = np.vstack((dCP, dCT))

        return J