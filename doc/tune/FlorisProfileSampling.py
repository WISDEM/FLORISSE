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

# visualization: define resolution
resolution = 75

# Define turbine characteristics
rotor_diameter = 126.4
rotorArea = np.pi*rotor_diameter*rotor_diameter/4.0
axial_induction = 1.0/3.0 # used only for initialization
generatorEfficiency = 0.944
hub_height = 90.0
NREL5MWCPCT = pickle.load(open('NREL5MWCPCT_dict.p'))
datasize = NREL5MWCPCT['CP'].size
turbineXinit = np.array([0.0])
turbineYinit = np.array([0.0])

nTurbines = turbineXinit.size

prob = Problem(root=AEPGroup(nTurbines=2, nDirections=1, differentiable=True, use_rotor_components=True,
                             datasize=datasize, nSamples=resolution))

prob.setup()

# load turbine properties into FLORIS
prob['gen_params:windSpeedToCPCT_CP'] = NREL5MWCPCT['CP']
prob['gen_params:windSpeedToCPCT_CT'] = NREL5MWCPCT['CT']
prob['gen_params:windSpeedToCPCT_wind_speed'] = NREL5MWCPCT['wind_speed']
prob['floris_params:cos_spread'] = 1E12
prob['axialInduction'] = np.ones(nTurbines)*axial_induction
prob['rotorDiameter'] = np.ones(nTurbines)*rotor_diameter
prob['hubHeight'] = np.array([hub_height, hub_height])
prob['generatorEfficiency'] = np.array([generatorEfficiency, generatorEfficiency])
prob['turbineX'] = turbineXinit
prob['turbineY'] = turbineYinit


# define site measurements
windDirection = 270.
prob['windDirections'] = np.array([windDirection])
wind_speed = 8.1    # m/s
prob['windSpeeds'] = np.array([wind_speed])
prob['air_density'] = 1.1716
wind_frequency = 1.0
prob['windFrequencies'] = np.array([wind_frequency])

# generate points at hub-height where velocities are to be sampled
x = np.ones(resolution)*7.0*rotor_diameter
y = np.linspace(-3.0*rotor_diameter, 3.0*rotor_diameter, resolution)

wsPositionX = x
wsPositionY = y
wsPositionZ = np.ones(wsPositionX.shape)*hub_height

# set yaw to zero
prob['yaw0'] = np.array([0.0, 0.0])

# Define turbine locations
prob['turbineX'] = turbineXinit
prob['turbineY'] = turbineYinit

# set sampling points in problem
prob['wsPositionX'] = np.copy(wsPositionX)
prob['wsPositionY'] = np.copy(wsPositionY)
prob['wsPositionZ'] = np.copy(wsPositionZ)

# Call FLORIS
prob.run()

# extract velocities at sampling points from problem
velocities = prob['wsArray0']

# visualize
plt.figure()
plt.plot(y/rotor_diameter, velocities/wind_speed)
plt.xlabel('Crosswind offset ($\Delta X / D_r$)')
plt.ylabel('Normalized Velocity ($V/V_0$)')
plt.show()




if __name__ == "__main__":
    plt.show()
else:
    plt.show(block=False)