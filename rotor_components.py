from openmdao.main.api import Component, VariableTree
from openmdao.lib.datatypes.api import Array, Float, Bool, Int, List, Str, VarTree
import numpy as np
from scipy import interp

## Components that interpolate a predefined CP/CT curve and apply a yaw correction

class windSpeedToCPCT(VariableTree):
    # variable tree that defines known wind speed to CP/CT curve
    def __init__(self, datasize=0):

        super(windSpeedToCPCT, self).__init__()

        self.add('wind_speed', Array(np.zeros(datasize), iotype='in', units='m/s', desc='range of wind speeds'))
        self.add('CP', Array(np.zeros(datasize), iotype='out', desc='power coefficients'))
        self.add('CT', Array(np.zeros(datasize), iotype='out', desc='thrust coefficients'))

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
        
        print "first speed, first CT", self.windSpeedToCPCT.wind_speed[0], self.windSpeedToCPCT.CT[0], self.windSpeedToCPCT.CT[1]
        print "last speed, last CT", self.windSpeedToCPCT.wind_speed[-1], self.windSpeedToCPCT.CP[-1]
        print "pP: ", self.pP
        wind_speed_ax = np.cos(self.yaw*np.pi/180.0)**(self.pP/3.0)*self.wind_speed_hub
        # use interpolation on precalculated CP-CT curve
        wind_speed_ax = np.maximum(wind_speed_ax, self.windSpeedToCPCT.wind_speed[0])
        wind_speed_ax = np.minimum(wind_speed_ax, self.windSpeedToCPCT.wind_speed[-1])
        self.CP = interp(wind_speed_ax, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CP)
        self.CT = interp(wind_speed_ax, self.windSpeedToCPCT.wind_speed, self.windSpeedToCPCT.CT)

        # normalize on incoming wind speed to correct coefficients for yaw
        self.CP = self.CP * np.cos(self.yaw*np.pi/180.0)**self.pP
        self.CT = self.CT * np.cos(self.yaw*np.pi/180.0)**2
        
        print "in rotor, Cp = ", self.CP
        print "in rotor, Ct = ", self.CT
        print "in rotor, pP = ", self.pP
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
