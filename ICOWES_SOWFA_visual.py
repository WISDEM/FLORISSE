import numpy as np
from matplotlib import pyplot as plt
from scipy.io import loadmat

import cPickle as pickle

from Original_assembly import floris_assembly_opt_AEP

from Parameters import FLORISParameters

ICOWESdata = loadmat('YawPosResults.mat')
yawrange = ICOWESdata['yaw'][0]

optimize_position = False
optimize_yaw = False
use_rotor_components = True

myFloris = floris_assembly_opt_AEP(nTurbines=2., nDirections=1, optimize_yaw=optimize_yaw,
                                   optimize_position=optimize_position, use_rotor_components=use_rotor_components,
                                   datasize=datasize)

myFloris.parameters = FLORISParameters()

rotorDiameter = 126.4
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0
CP = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
CT = 4.0*axialInduction*(1.0-axialInduction)
generator_efficiency = 0.944

myFloris.parameters.CPcorrected = True
myFloris.parameters.CTcorrected = True
myFloris.parameters.FLORISoriginal = False

# Define turbine characteristics
myFloris.axialInduction = np.array([axialInduction, axialInduction])
myFloris.rotorDiameter = np.array([rotorDiameter, rotorDiameter])
myFloris.rotorArea = np.array([rotorArea, rotorArea])

myFloris.generator_efficiency = np.array([generator_efficiency, generator_efficiency])

# Define site measurements
myFloris.windrose_directions = np.array([30.])
wind_speed = 8.1    # m/s
myFloris.windrose_speeds = np.ones_like(myFloris.windrose_directions)*wind_speed
myFloris.air_density = 1.1716

if use_rotor_components:
    myFloris.initVelocitiesTurbines = np.ones_like(myFloris.windrose_directions)*wind_speed
    myFloris.curve_CP = NREL5MWCPCT.CP
    myFloris.curve_CT = NREL5MWCPCT.CT
    myFloris.curve_wind_speed = NREL5MWCPCT.wind_speed
    myFloris.parameters.ke = 0.05
    myFloris.parameters.kd = 0.17
    myFloris.parameters.aU = 12.0
    myFloris.parameters.bU = 1.3
    myFloris.parameters.initialWakeAngle = 3.0
    myFloris.parameters.useaUbU = True
    myFloris.parameters.useWakeAngle = True
    myFloris.parameters.adjustInitialWakeDiamToYaw = False
else:
    myFloris.Ct = np.array([CT, CT])
    myFloris.Cp = np.array([CP, CP])

FLORISpower = list()
FLORISgradient = list()

for yaw1 in yawrange:

    # Defube turbine locations and orientation
    myFloris.turbineX = np.array([1118.1, 1881.9])
    myFloris.turbineY = np.array([1279.5, 1720.5])

    myFloris.yaw = np.array([yaw1, 0.0])

    # Call FLORIS
    myFloris.run()

    FLORISpower.append(myFloris.floris_power_0.wt_power)

FLORISpower = np.array(FLORISpower)
SOWFApower = np.array([ICOWESdata['yawPowerT1'][0],ICOWESdata['yawPowerT2'][0]]).transpose()/1000.

fig, axes = plt.subplots(ncols = 2, sharey = True)
axes[0].plot(yawrange.transpose(), FLORISpower[:,0], 'r-', yawrange.transpose(), SOWFApower[:,0], 'ro')
axes[0].plot(yawrange.transpose(), FLORISpower[:,1], 'b-', yawrange.transpose(), SOWFApower[:,1], 'bo')
axes[0].plot(yawrange.transpose(), FLORISpower[:,0]+FLORISpower[:,1], 'k-', yawrange.transpose(), SOWFApower[:,0]+SOWFApower[:,1], 'ko')

error_turbine2 = np.sum(np.abs(FLORISpower[:, 1] - SOWFApower[:, 1]))

posrange = ICOWESdata['pos'][0]

myFloris.yaw = np.array([0.0, 0.0])
FLORISpower = list()
for pos2 in posrange:
    # Defube turbine locations and orientation
    effUdXY = 0.523599

    Xinit = np.array([1118.1, 1881.9])
    Yinit = np.array([1279.5, 1720.5])
    XY = np.array([Xinit, Yinit]) + np.dot(np.array([[np.cos(effUdXY),-np.sin(effUdXY)], [np.sin(effUdXY),np.cos(effUdXY)]]), np.array([[0., 0], [0,pos2]]))
    myFloris.turbineX = XY[0,:]
    myFloris.turbineY = XY[1,:]

    yaw = np.array([0.0, 0.0])

    # Call FLORIS
    myFloris.run()

    FLORISpower.append(myFloris.floris_power_0.wt_power)

FLORISpower = np.array(FLORISpower)
SOWFApower = np.array([ICOWESdata['posPowerT1'][0],ICOWESdata['posPowerT2'][0]]).transpose()/1000.

error_turbine2 += np.sum(np.abs(FLORISpower[:,1] - SOWFApower[:,1]))

print error_turbine2

axes[1].plot(posrange, FLORISpower[:,0], 'r-', posrange, SOWFApower[:,0], 'ro')
axes[1].plot(posrange, FLORISpower[:,1], 'b-', posrange, SOWFApower[:,1], 'bo')
axes[1].plot(posrange, FLORISpower[:,0]+FLORISpower[:,1], 'k-', posrange, SOWFApower[:,0]+SOWFApower[:,1], 'ko')

# visualize

x = np.linspace(0,3000,200)
y = np.linspace(0,3000,200)
xx, yy = np.meshgrid(x, y)
ws_positionX = xx.flatten()
ws_positionY = yy.flatten()
myFlorisV = floris_assembly_opt_AEP(nTurbines=2., nDirections=1, optimize_yaw=optimize_yaw,
                                   optimize_position=optimize_position, use_rotor_components=use_rotor_components,
                                   datasize=datasize, nSamples =len(ws_positionX))
myFlorisV.ws_positionX = ws_positionX
myFlorisV.ws_positionY = ws_positionY

myFlorisV.parameters = myFloris.parameters
myFlorisV.yaw = myFloris.yaw
myFlorisV.axialInduction = myFloris.axialInduction
myFlorisV.rotorDiameter = myFloris.rotorDiameter
myFlorisV.rotorArea = myFloris.rotorArea
myFlorisV.generator_efficiency = myFloris.generator_efficiency

myFlorisV.turbineX = myFloris.turbineX
myFlorisV.turbineY = myFloris.turbineY
myFlorisV.curve_CP = myFloris.curve_CP
myFlorisV.curve_CT = myFloris.curve_CT
myFlorisV.curve_wind_speed = myFloris.curve_wind_speed
myFlorisV.initVelocitiesTurbines = myFloris.initVelocitiesTurbines
myFlorisV.windrose_directions = myFloris.windrose_directions
myFlorisV.windrose_speeds = myFloris.windrose_speeds
myFlorisV.air_density = myFloris.air_density

myFlorisV.run()
velocities = myFlorisV.ws_array_0
velocities = velocities.reshape(len(y), len(x))

fig, (ax1) = plt.subplots(nrows=1)
im = ax1.pcolormesh(x, y, velocities, cmap='coolwarm')
cbar = plt.colorbar(im, orientation='vertical')
cbar.set_label('Wind Speed (m/s)', rotation=270, labelpad=30)
ax1.set_aspect('equal')
ax1.autoscale(tight=True)

plt.show()