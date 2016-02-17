from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

from rotor_components import CPCT_Interpolate, windSpeedToCPCT
import numpy as np

import cPickle as pickle

from scipy.io import loadmat

d = loadmat('matWakeAnalysis/turbineData.mat')

datasize = len(d['windSpeedI'].flatten())

windSpeedToCPCT = windSpeedToCPCT(datasize = datasize)
windSpeedToCPCT.wind_speed = d['windSpeedI'].flatten()
windSpeedToCPCT.CP = d['rotCpI'].flatten()
windSpeedToCPCT.CT = d['rotCtI'].flatten()


rangeYaw = np.arange(-50.,50.,1)
windspeedRange = np.linspace(2.,windSpeedToCPCT.wind_speed[-1],1000)
power = np.zeros((rangeYaw.size, windspeedRange.size))
rotor_diameter = 77.0
A = np.pi*(rotor_diameter/2)**2
air_density = 1.1716

CPCT_Interpolate = CPCT_Interpolate(nTurbines = 1, datasize = datasize)
CPCT_Interpolate.windSpeedToCPCT = windSpeedToCPCT
CPCT_Interpolate.pP = 1.88
gen_eff = 0.95

for yawI, yaw in enumerate(rangeYaw): 
        for wind_speedI, wind_speed in enumerate(windspeedRange):
            CPCT_Interpolate.wind_speed_hub = np.array([wind_speed])
            CPCT_Interpolate.yaw = np.array([yaw])
            CPCT_Interpolate.run()
            power[yawI, wind_speedI] = gen_eff * 0.5 * air_density * CPCT_Interpolate.CP * wind_speed**3 * A

rangeYaw, windspeedRange = np.meshgrid(rangeYaw, windspeedRange)

fig = plt.figure(figsize=(8, 6), facecolor='white')
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(rangeYaw, windspeedRange, power.transpose()/1.e6, rstride=10, cstride=4, linewidth=0.2, cmap = cm.coolwarm, rasterized=True)
ax.azim=-45.
ax.elev=23.
ax.set_xlabel('yaw (deg)')
ax.set_ylabel('wind speed (m/s)')
ax.set_zlabel('turbine power production (MW)')
ax.autoscale(enable=True, axis='both', tight=True)

fig.tight_layout()
pos1 = ax.get_position() # get the original position 
pos2 = [pos1.x0 - 0.05, pos1.y0,  pos1.width / 1.0, pos1.height / 1.0] 
ax.set_position(pos2) # set a new position

# fig.savefig('YawWindSpeedvsPower.png')
# fig.savefig('YawWindSpeedvsPower.pdf')

rangeYaw = np.arange(-50.,50.,1)
windspeedRange2 = np.linspace(4.,windSpeedToCPCT.wind_speed[-1],1000)
thrustC = np.zeros((rangeYaw.size, windspeedRange2.size))

for yawI, yaw in enumerate(rangeYaw): 
        for wind_speedI, wind_speed in enumerate(windspeedRange2):
            CPCT_Interpolate.wind_speed_hub = np.array([wind_speed])
            CPCT_Interpolate.yaw = np.array([yaw])
            CPCT_Interpolate.run()
            thrustC[yawI, wind_speedI] = CPCT_Interpolate.CT

rangeYaw, windspeedRange2 = np.meshgrid(rangeYaw, windspeedRange2)

fig = plt.figure(figsize=(8, 6), facecolor='white')
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(rangeYaw, windspeedRange2, thrustC.transpose(), rstride=10, cstride=4, linewidth=0.2, cmap = cm.coolwarm, rasterized=True)
ax.azim=45.
ax.elev=23.
ax.autoscale(enable=True, axis='both', tight=True)
ax.set_xlabel('yaw (deg)')
ax.set_ylabel('wind speed (m/s)')
ax.set_zlabel('rotor thrust coefficient (-)')
fig.tight_layout()
pos1 = ax.get_position() # get the original position 
pos2 = [pos1.x0 + 0.05, pos1.y0,  pos1.width / 1.0, pos1.height / 1.0] 
ax.set_position(pos2) # set a new position
#fig.savefig('YawWindSpeedvsCT.png')
#fig.savefig('YawWindSpeedvsCT.pdf')

plt.show()
