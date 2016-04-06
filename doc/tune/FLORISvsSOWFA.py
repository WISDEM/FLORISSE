# This example script compares FLORIS predictions with steady-state SOWFA data as obtained 
# throught the simulations described in:
#   

import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt

from scipy.io import loadmat
import pickle

from florisse.floris import AEPGroup
from florisse import config

from openmdao.api import Problem

# set temp option to use unified floris
config.floris_single_component = True

# Load steady-state power data from SOWFA 
ICOWESdata = loadmat('YawPosResults.mat')

# visualization: define resolution
resolution = 75

# Define turbine characteristics
rotorDiameter = 126.4
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0 # used only for initialization
generator_efficiency = 0.944
hub_height = 90.0
NREL5MWCPCT = pickle.load(open('NREL5MWCPCT_dict.p'))
datasize = NREL5MWCPCT['CP'].size
turbineXinit = np.array([1118.1, 1881.9])
turbineYinit = np.array([1279.5, 1720.5])

prob = Problem(root=AEPGroup(nTurbines=2, nDirections=1, differentiable=True, use_rotor_components=True,
                                   datasize=datasize, nSamples=resolution*resolution))

prob.setup()

# load turbine properties into FLORIS
prob['gen_params:windSpeedToCPCT_CP'] = NREL5MWCPCT['CP']
prob['gen_params:windSpeedToCPCT_CT'] = NREL5MWCPCT['CT']
prob['gen_params:windSpeedToCPCT_wind_speed'] = NREL5MWCPCT['wind_speed']
prob['axialInduction'] = np.array([axialInduction, axialInduction])
prob['rotorDiameter'] = np.array([rotorDiameter, rotorDiameter])
prob['hubHeight'] = np.array([hub_height, hub_height])
prob['generator_efficiency'] = np.array([generator_efficiency, generator_efficiency])
prob['turbineX'] = turbineXinit
prob['turbineY'] = turbineYinit


# Define site measurements
windDirection = 270.-0.523599*180./np.pi
prob['windDirections'] = np.array([windDirection])
wind_speed = 8.1    # m/s
prob['windSpeeds'] = np.array([wind_speed])
prob['air_density'] = 1.1716
wind_frequency = 1.0
prob['windFrequencies'] = np.array([wind_frequency])

# prob.initVelocitiesTurbines = np.ones_like(prob.windrose_directions)*wind_speed

windDirectionPlot = 270. - windDirection
# visualization:
#    generate points downstream slice
y_cut = np.linspace(-rotorDiameter, rotorDiameter, resolution)
z_cut = np.linspace(-hub_height, rotorDiameter, resolution)
yy, zz = np.meshgrid(y_cut, z_cut)
xx = np.ones(yy.shape) * 3.5*rotorDiameter
position = np.array([xx.flatten(), yy.flatten(), zz.flatten()])
rotationMatrix = np.array([(np.cos(windDirectionPlot*np.pi/180.), -np.sin(windDirectionPlot*np.pi/180.), 0.),
                                   (np.sin(windDirectionPlot*np.pi/180.), np.cos(windDirectionPlot*np.pi/180.), 0.),
                                   (0., 0., 1.)])
positionF = np.dot(rotationMatrix, position) + np.dot(np.array([(prob['turbineX'][0], prob['turbineY'][0], hub_height)]).transpose(), np.ones((1, np.size(position, 1))))

#      generate points hub-height
x = np.linspace(750, 2400, resolution)
y = np.linspace(750, 2400, resolution)
xx, yy = np.meshgrid(x, y)
wsPositionX = xx.flatten()
wsPositionY = yy.flatten()
wsPositionZ = np.ones(wsPositionX.shape)*hub_height

# SWEEP TURBINE YAW
FLORISpower = list()
yawrange = ICOWESdata['yaw'][0]

velocities = list()
velocities_cut = list()

for yaw1 in yawrange:

    prob['yaw0'] = np.array([yaw1, 0.0])

    # Call FLORIS horizontal slice
    prob['wsPositionX'] = np.copy(wsPositionX)
    prob['wsPositionY'] = np.copy(wsPositionY)
    prob['wsPositionZ'] = np.copy(wsPositionZ)
    prob.run()
    FLORISpower.append(np.array(prob['wtPower0']))
    velocities.append(np.array(prob['wsArray0']))

   # Call FLORIS cut-through slice
    prob['wsPositionX'] = np.copy(positionF[0])
    prob['wsPositionY'] = np.copy(positionF[1])
    prob['wsPositionZ'] = np.copy(positionF[2])
    prob.run()
    velocities_cut.append(np.array(prob['wsArray0']))

# plot slices
velocities = np.array(velocities)
vmin = np.min(velocities)
vmax = np.max(velocities)
velocities_cut = np.array(velocities_cut)

fig, axes = plt.subplots(ncols=int(np.ceil(len(yawrange)/2.)), nrows=4, figsize=(23, 12))
fig.suptitle("FLORIS flow-field prediction at hub-height and wake cut-through at 3.5D, for yaw sweep")

axes1 = list(axes[0])+list(axes[2])
axes2 = list(axes[1])+list(axes[3])

for i in range(len(yawrange)):
    vel = velocities[i].flatten()
    vel = vel.reshape(len(y), len(x))
    ax1 = axes1[i]
    im = ax1.pcolormesh(x, y, vel, cmap='coolwarm', vmin=vmin, vmax=vmax)
    ax1.set_aspect('equal')
    ax1.set_xticks(np.arange(800, 2800, 400))
    ax1.set_yticks(np.arange(800, 2800, 400))
    ax1.autoscale(tight=True)
    ax1.set_title('front turbine yawed %d deg' % yawrange[i])

    vel = velocities_cut[i].flatten()
    vel = vel.reshape(len(z_cut), len(y_cut))
    ax2 = axes2[i]
    im = ax2.pcolormesh(y_cut, z_cut, vel, cmap='coolwarm', vmin=vmin, vmax=vmax)
    ax2.set_aspect('equal')
    ax2.autoscale(tight=True)
    ax2.invert_xaxis()

cbar = plt.colorbar(im, orientation = 'horizontal', ticks=[vmin,(vmin+vmax)/2,vmax])
cbar.set_label('wind speed (m/s)')
axes1[-1].axis('off')
axes2[-1].axis('off')
mpl.rcParams.update({'font.size': 8})
plt.tight_layout()
fig.subplots_adjust(top=0.95)

FLORISpower = np.array(FLORISpower)
SOWFApower = np.array([ICOWESdata['yawPowerT1'][0],ICOWESdata['yawPowerT2'][0]]).transpose()/1000.

figPower, axesPower = plt.subplots(ncols=2, sharey=True)
axesPower[0].plot(yawrange.transpose(), FLORISpower[:, 0], 'r-', yawrange.transpose(), SOWFApower[:, 0], 'ro')
axesPower[0].plot(yawrange.transpose(), FLORISpower[:, 1], 'b-', yawrange.transpose(), SOWFApower[:, 1], 'bo')
axesPower[0].plot(yawrange.transpose(), FLORISpower[:, 0]+FLORISpower[:, 1], 'k-', yawrange.transpose(), SOWFApower[:, 0]+SOWFApower[:, 1], 'ko')
axesPower[0].set_xlabel('yaw front turbine 1 (deg)')
axesPower[0].set_ylabel('power (kW)')
axesPower[0].legend(['front turbine FLORIS', 'front turbine SOWFA', 'back turbine FLORIS',  'back turbine SOWFA', 'total FLORIS', 'total SOWFA'])

# SWEEP TURBINE POSITIONS
posrange = ICOWESdata['pos'][0]
prob['yaw0'] = np.array([0.0, 0.0])
FLORISpower = list()

velocities = list()
velocities_cut = list()

for pos2 in posrange:
    
    # Define turbine locations and orientation
    effUdXY = 0.523599
    XY = np.array([turbineXinit, turbineYinit]) + np.dot(np.array([[np.cos(effUdXY),-np.sin(effUdXY)], [np.sin(effUdXY),np.cos(effUdXY)]]), np.array([[0., 0], [0,pos2]]))
    prob['turbineX'] = XY[0, :]
    prob['turbineY'] = XY[1, :]

    # Call FLORIS horizontal slice
    prob['wsPositionX'] = np.copy(wsPositionX)
    prob['wsPositionY'] = np.copy(wsPositionY)
    prob['wsPositionZ'] = np.copy(wsPositionZ)
    prob.run()
    FLORISpower.append(np.array(prob['wtPower0']))
    velocities.append(np.array(prob['wsArray0']))


    # Call FLORIS cut-through slice
    prob['wsPositionX'] = np.copy(positionF[0])
    prob['wsPositionY'] = np.copy(positionF[1])
    prob['wsPositionZ'] = np.copy(positionF[2])
    print "CUT: "
    prob.run()
    print "CUT END"
    velocities_cut.append(np.array(prob['wsArray0']))
    print prob['wsArray0'][len(prob['wsArray0'])/2]

# plot powers
FLORISpower = np.array(FLORISpower)
SOWFApower = np.array([ICOWESdata['posPowerT1'][0], ICOWESdata['posPowerT2'][0]]).transpose()/1000.
# print FLORISpower
axesPower[1].plot(posrange, FLORISpower[:, 0], 'r-', posrange, SOWFApower[:, 0], 'ro')
axesPower[1].plot(posrange, FLORISpower[:, 1], 'b-', posrange, SOWFApower[:, 1], 'bo')
axesPower[1].plot(posrange, FLORISpower[:, 0]+FLORISpower[:, 1], 'k-', posrange, SOWFApower[:, 0]+SOWFApower[:, 1], 'ko')
axesPower[1].set_xlabel('back turbine displacement (m)')
axesPower[1].set_ylabel('power (kW)')

# plot slices
velocities = np.array(velocities)
vmin = np.min(velocities)
vmax = np.max(velocities)
velocities_cut = np.array(velocities_cut)

fig, axes = plt.subplots(ncols=int(np.ceil(len(posrange)/2.)), nrows=4, figsize=(23, 12))
fig.suptitle("FLORIS flow-field prediction at hub-height and wake cut-through at 3.5D, for yaw sweep")

axes1 = list(axes[0])+list(axes[2])
axes2 = list(axes[1])+list(axes[3])

for i in range(len(posrange)):
    vel = velocities[i].flatten()
    vel = vel.reshape(len(y), len(x))
    ax1 = axes1[i]
    im = ax1.pcolormesh(x, y, vel, cmap='coolwarm', vmin=vmin, vmax=vmax)
    ax1.set_aspect('equal')
    ax1.autoscale(tight=True)
    ax1.set_xticks(np.arange(800, 2800, 400))
    ax1.set_yticks(np.arange(800, 2800, 400))
    ax1.set_title('back turbine moved %d m' % posrange[i])

    vel = velocities_cut[i].flatten()
    vel = vel.reshape(len(z_cut), len(y_cut))
    ax2 = axes2[i]
    im = ax2.pcolormesh(y_cut, z_cut, vel, cmap='coolwarm', vmin=vmin, vmax=vmax)
    ax2.set_aspect('equal')
    ax2.autoscale(tight=True)
    ax2.invert_xaxis()

cbar = plt.colorbar(im, orientation='horizontal', ticks=[vmin, (vmin+vmax)/2, vmax])
cbar.set_label('wind speed (m/s)')
axes1[-1].axis('off')
axes2[-1].axis('off')
mpl.rcParams.update({'font.size': 8})
plt.tight_layout()
fig.subplots_adjust(top=0.95)


if __name__ == "__main__":
    plt.show()
else:
    plt.show(block=False)