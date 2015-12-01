import numpy as np
from matplotlib import pyplot as plt
from scipy.io import loadmat
from openmdao.api import Problem
import cPickle as pickle

from floris_openmdao1 import AEPGroupFLORIS

from Parameters import FLORISParameters

ICOWESdata = loadmat('YawPosResults.mat')
yawrange = ICOWESdata['yaw'][0]

optimize_position = False
optimize_yaw = False
use_rotor_components = True

# NREL5MWCPCT = pickle.load(open('NREL5MWCPCT.p'))
# datasize = NREL5MWCPCT.CP.size
myFloris = Problem(root=AEPGroupFLORIS(nTurbines=2., nDirections=1, resolution=0.0))
myFloris.setup()

# myFloris.parameters = FLORISParameters()

rotorDiameter = 126.4
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0
CP = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
CT = 4.0*axialInduction*(1.0-axialInduction)
generator_efficiency = 0.944

myFloris['floris_params:CPcorrected'] = False
myFloris['floris_params:CTcorrected'] = False
myFloris['floris_params:FLORISoriginal'] = True

# Define turbine characteristics
myFloris['axialInduction'] = np.array([axialInduction, axialInduction])
myFloris['rotorDiameter'] = np.array([rotorDiameter, rotorDiameter])
# myFloris['rotorArea'] = np.array([rotorArea, rotorArea])

myFloris['generator_efficiency'] = np.array([generator_efficiency, generator_efficiency])

# Define site measurements
myFloris['windDirections'] = np.array([240.])
myFloris['wind_speed'] = 8.1    # m/s
# myFloris['windrose_speeds'] = np.ones_like(myFloris['windDirections'])*myFloris['wind_speed']
myFloris['air_density'] = 1.1716


# # assign initial values to design variables
#     prob['turbineX'] = turbineX
#     prob['turbineY'] = turbineY
#     prob['yaw'] = yaw
#
#     # assign values to constant inputs (not design variables)
#     prob['rotorDiameter'] = rotorDiameter
#     prob['axialInduction'] = axialInduction
#     prob['generator_efficiency'] = generator_efficiency
#     prob['wind_speed'] = wind_speed
#     prob['air_density'] = air_density
#     prob['wind_direction'] = wind_direction
#     prob['Ct_in'] = Ct
#     prob['Cp_in'] = Cp
#     prob['floris_params:FLORISoriginal'] = True
#     prob['floris_params:CPcorrected'] = False
#     prob['floris_params:CTcorrected'] = False

# if use_rotor_components:
#     myFloris.wind_speed = 8.1
#     myFloris.initVelocitiesTurbines = np.ones_like(myFloris.windrose_directions)*myFloris.wind_speed
#     myFloris.curve_CP = NREL5MWCPCT.CP
#     myFloris.curve_CT = NREL5MWCPCT.CT
#     myFloris.curve_wind_speed = NREL5MWCPCT.wind_speed
# else:
myFloris['Ct_in'] = np.array([CT, CT])
myFloris['Cp_in'] = np.array([CP, CP])

FLORISpower = list()
FLORISgradient = list()

for yaw1 in yawrange:

    # Defube turbine locations and orientation
    myFloris['turbineX'] = np.array([1118.1, 1881.9])
    myFloris['turbineY'] = np.array([1279.5, 1720.5])

    myFloris['yaw0'] = np.array([yaw1, 0.0])

    # Call FLORIS
    myFloris.run()
    # print 'power = ', myFloris.root.dir0.unknowns['wt_power']
    FLORISpower.append(list(myFloris.root.dir0.unknowns['wt_power']))
    # print FLORISpower
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
    myFloris['turbineX'] = XY[0,:]
    myFloris['turbineY'] = XY[1,:]

    yaw = np.array([0.0, 0.0])



    # Call FLORIS
    myFloris.run()
    print 'power = ', myFloris.root.dir0.unknowns['wt_power']
    FLORISpower.append(list(myFloris.root.dir0.unknowns['wt_power']))

FLORISpower = np.array(FLORISpower)
SOWFApower = np.array([ICOWESdata['posPowerT1'][0],ICOWESdata['posPowerT2'][0]]).transpose()/1000.

error_turbine2 += np.sum(np.abs(FLORISpower[:,1] - SOWFApower[:,1]))

print error_turbine2

axes[1].plot(posrange, FLORISpower[:,0], 'r-', posrange, SOWFApower[:,0], 'ro')
axes[1].plot(posrange, FLORISpower[:,1], 'b-', posrange, SOWFApower[:,1], 'bo')
axes[1].plot(posrange, FLORISpower[:,0]+FLORISpower[:,1], 'k-', posrange, SOWFApower[:,0]+SOWFApower[:,1], 'ko')

plt.show()

import numpy as np
from matplotlib import pyplot as plt
from scipy.io import loadmat
from openmdao.api import Problem
import cPickle as pickle

from floris_openmdao1 import AEPGroupFLORIS

from Parameters import FLORISParameters

ICOWESdata = loadmat('YawPosResults.mat')
yawrange = ICOWESdata['yaw'][0]

optimize_position = False
optimize_yaw = False
use_rotor_components = True

# NREL5MWCPCT = pickle.load(open('NREL5MWCPCT.p'))
# datasize = NREL5MWCPCT.CP.size
myFloris = Problem(root=AEPGroupFLORIS(nTurbines=2., nDirections=1, resolution=0.0))
myFloris.setup()

# myFloris.parameters = FLORISParameters()

rotorDiameter = 126.4
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0
CP = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
CT = 4.0*axialInduction*(1.0-axialInduction)
generator_efficiency = 0.944

myFloris['floris_params:CPcorrected'] = False
myFloris['floris_params:CTcorrected'] = False
myFloris['floris_params:FLORISoriginal'] = True

# Define turbine characteristics
myFloris['axialInduction'] = np.array([axialInduction, axialInduction])
myFloris['rotorDiameter'] = np.array([rotorDiameter, rotorDiameter])
# myFloris['rotorArea'] = np.array([rotorArea, rotorArea])

myFloris['generator_efficiency'] = np.array([generator_efficiency, generator_efficiency])

# Define site measurements
myFloris['windDirections'] = np.array([240.])
myFloris['wind_speed'] = 8.    # m/s
# myFloris['windrose_speeds'] = np.ones_like(myFloris['windDirections'])*myFloris['wind_speed']
myFloris['air_density'] = 1.1716


# # assign initial values to design variables
#     prob['turbineX'] = turbineX
#     prob['turbineY'] = turbineY
#     prob['yaw'] = yaw
#
#     # assign values to constant inputs (not design variables)
#     prob['rotorDiameter'] = rotorDiameter
#     prob['axialInduction'] = axialInduction
#     prob['generator_efficiency'] = generator_efficiency
#     prob['wind_speed'] = wind_speed
#     prob['air_density'] = air_density
#     prob['wind_direction'] = wind_direction
#     prob['Ct_in'] = Ct
#     prob['Cp_in'] = Cp
#     prob['floris_params:FLORISoriginal'] = True
#     prob['floris_params:CPcorrected'] = False
#     prob['floris_params:CTcorrected'] = False

# if use_rotor_components:
#     myFloris.wind_speed = 8.1
#     myFloris.initVelocitiesTurbines = np.ones_like(myFloris.windrose_directions)*myFloris.wind_speed
#     myFloris.curve_CP = NREL5MWCPCT.CP
#     myFloris.curve_CT = NREL5MWCPCT.CT
#     myFloris.curve_wind_speed = NREL5MWCPCT.wind_speed
# else:
myFloris['Ct_in'] = np.array([CT, CT])
myFloris['Cp_in'] = np.array([CP, CP])

FLORISpower = list()
FLORISgradient = list()

for yaw1 in yawrange:

    # Defube turbine locations and orientation
    myFloris['turbineX'] = np.array([1118.1, 1881.9])
    myFloris['turbineY'] = np.array([1279.5, 1720.5])

    myFloris['yaw0'] = np.array([yaw1, 0.0])

    # Call FLORIS
    myFloris.run()

    FLORISpower.append(list(myFloris.root.dir0.unknowns['wt_power']))

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
    print XY
    myFloris['turbineX'] = XY[0,:]
    myFloris['turbineY'] = XY[1,:]

    myFloris['yaw0'] = np.array([0.0, 0.0])



    # Call FLORIS
    myFloris.run()
    # print 'power = ', myFloris.root.dir0.unknowns['wt_power']
    FLORISpower.append(list(myFloris.root.dir0.unknowns['wt_power']))
FLORISpower = np.array(FLORISpower)
SOWFApower = np.array([ICOWESdata['posPowerT1'][0],ICOWESdata['posPowerT2'][0]]).transpose()/1000.

error_turbine2 += np.sum(np.abs(FLORISpower[:,1] - SOWFApower[:,1]))

print error_turbine2

axes[1].plot(posrange, FLORISpower[:,0], 'r-', posrange, SOWFApower[:,0], 'ro')
axes[1].plot(posrange, FLORISpower[:,1], 'b-', posrange, SOWFApower[:,1], 'bo')
axes[1].plot(posrange, FLORISpower[:,0]+FLORISpower[:,1], 'k-', posrange, SOWFApower[:,0]+SOWFApower[:,1], 'ko')

plt.show()

FLORISeffu = list()
y = np.linspace(-1.5*rotorDiameter, 1.5*rotorDiameter, 100)

for i in range(0, 100):
    myFloris['windDirections'] = np.array([270])
    X = np.array([0, 20.*rotorDiameter])
    Y = np.array([0, y[i]])
    myFloris['turbineX'] = X
    myFloris['turbineY'] = Y
    myFloris['yaw0'] = np.array([0.0, 0.0])

    # Call FLORIS
    myFloris.run()

    FLORISeffu.append(list(myFloris.root.dir0.unknowns['velocitiesTurbines']))

FLORISeffu = np.array(FLORISeffu)
plt.figure()
plt.plot(y/rotorDiameter, FLORISeffu[:, 1])
plt.show()

FLORISindiam = list()
res = 1000
x = np.linspace(-0.25*rotorDiameter, 20.0*rotorDiameter, res)
for i in range(0, res):
    myFloris['windDirections'] = np.array([270])
    X = np.array([0, x[i]])
    Y = np.array([0, myFloris['wakeCentersYT'][2]])
    print myFloris['wakeCentersYT'][2]
    # Y = np.array([0, 0])
    myFloris['turbineX'] = X
    myFloris['turbineY'] = Y
    myFloris['yaw0'] = np.array([0.0, 0.0])

    # Call FLORIS
    myFloris.run()

    FLORISindiam.append(list(myFloris.root.dir0.unknowns['velocitiesTurbines']))

FLORISindiam = np.array(FLORISindiam)
plt.figure()
plt.plot(x/rotorDiameter, FLORISindiam[:, 1])
plt.show()
