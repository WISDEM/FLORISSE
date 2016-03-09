from matplotlib import pyplot as plt
import cPickle as pickle
import numpy as np
from Circle_assembly import floris_assembly_opt_AEP
from scipy.io import loadmat

import time

useSubset = True

# optimize yaw for all or a subset of NREL 5MW turbines in Princess Amalia Wind Farm configuration, for one wind speed and direction
baselinePowers = list()
optPowers = list()
increasePercentages = list()
optYaws = list()

AmaliaLocationsAndHull = loadmat('Amalia_locAndHull.mat')
if useSubset:
    turbineX = AmaliaLocationsAndHull['turbineX'].flatten()[:10]
    turbineY = AmaliaLocationsAndHull['turbineY'].flatten()[:10]
    maxiter = 100
else:
    turbineX = AmaliaLocationsAndHull['turbineX'].flatten()
    turbineY = AmaliaLocationsAndHull['turbineY'].flatten()
    maxiter = 5

nTurbines= len(turbineX)

rotorDiameter = 126.4
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0 # used only for initialization
generator_efficiency = 0.944
hub_height = 90.0
NREL5MWCPCT = pickle.load(open('NREL5MWCPCT.p'))
datasize = NREL5MWCPCT.CP.size

windSpeed = 8.0
windDirectionsWindRose = np.array([0.])
windDirectionsFLORIS = 270.-windDirectionsWindRose
airDensity = 1.1716

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

for windDirectionI, windDirection in enumerate(windDirectionsFLORIS):
    
    print 'wind direction %f deg' % windDirectionsWindRose[windDirectionI]

    baselineFloris = floris_assembly_opt_AEP(nTurbines=nTurbines, nDirections=1, optimize_yaw=False, optimize_position=False, datasize=datasize, nSamples = 0, nSpeeds = 1)
    baselineFloris.windrose_directions    = np.array([windDirection]);
    baselineFloris.initVelocitiesTurbines = np.copy(np.ones_like(baselineFloris.windrose_directions)*windSpeed)
    baselineFloris.windrose_speeds        = windSpeed
    baselineFloris.air_density            = airDensity
    baselineFloris.curve_wind_speed       = NREL5MWCPCT.wind_speed
    baselineFloris.curve_CP               = NREL5MWCPCT.CP
    baselineFloris.curve_CT               = NREL5MWCPCT.CT
    baselineFloris.rotorDiameter          = np.ones(nTurbines)*rotorDiameter
    baselineFloris.rotorArea              = np.ones(nTurbines)*rotorArea
    baselineFloris.turbineX               = turbineX
    baselineFloris.turbineY               = turbineY
    baselineFloris.axialInduction         = np.ones(nTurbines)*axialInduction # values used for initialization only
    baselineFloris.hubHeight              = np.ones(nTurbines)*hub_height
    baselineFloris.generator_efficiency   = np.ones(nTurbines)*generator_efficiency

    baselineFloris.run()
    baselinePower = np.sum(baselineFloris.floris_power_0.wt_power)
    baselinePowers.append(baselinePower)
    
    print 'baseline power %f kW' % baselinePower
    
    optFloris = floris_assembly_opt_AEP(nTurbines=nTurbines, nDirections=1, optimize_yaw=True, optimize_position=False, datasize=datasize, nSamples = 0, nSpeeds = 1, maxiter = maxiter)

    optFloris.windrose_directions    = np.array([windDirection]);
    optFloris.initVelocitiesTurbines = np.copy(np.ones_like(baselineFloris.windrose_directions)*windSpeed)
    optFloris.windrose_speeds        = windSpeed
    optFloris.air_density            = airDensity
    optFloris.curve_wind_speed       = NREL5MWCPCT.wind_speed
    optFloris.curve_CP               = NREL5MWCPCT.CP
    optFloris.curve_CT               = NREL5MWCPCT.CT
    optFloris.rotorDiameter          = np.ones(nTurbines)*rotorDiameter
    optFloris.rotorArea              = np.ones(nTurbines)*rotorArea
    optFloris.turbineX               = turbineX
    optFloris.turbineY               = turbineY
    optFloris.axialInduction         = np.ones(nTurbines)*axialInduction # values used for initialization only
    optFloris.hubHeight              = np.ones(nTurbines)*hub_height
    optFloris.generator_efficiency   = np.ones(nTurbines)*generator_efficiency

    tic = time.time()
    optFloris.run()
    toc = time.time()

    optPower = np.sum(optFloris.floris_power_0.wt_power)
    optPowers.append(optPower)

    increasePercentage = 100 * (optPower - baselinePower) / baselinePower
    increasePercentages.append(increasePercentage)
    
    optYaw = optFloris.yaw_0
    optYaws.append(optYaw)

    print 'optimal yaw %s deg' % optYaw 
    print 'optimal power %f kW' % optPower

    print 'increase %f %%' % increasePercentage
    print '----------------------'

baselinePowers = np.array(baselinePowers)
optPowers = np.array(optPowers)
increasePercentages = np.array(increasePercentages)

print('FLORIS yaw optimization took %.03f sec.' % (toc-tic))

# # VISUALIZE THE FLOW (full farm) for highest increase case - baseline and optimized

# get the case with highest increase
caseMaxIncrease = np.argmax(increasePercentages)
optYawMaxIncrease = optYaws[caseMaxIncrease]
windDirectionFLORISmaxIncrease = windDirectionsFLORIS[caseMaxIncrease]
windDirectionsWindRosemaxIncrease = windDirectionsWindRose[caseMaxIncrease]

# ... define sampling field
resolution = 300
xSamples = np.linspace(np.min(turbineX-200.), np.max(turbineX+200.), resolution)
ySamples = np.linspace(np.min(turbineY-200.), np.max(turbineY+200.), resolution)
nSamples = resolution**2
xxSamples, yySamples = np.meshgrid(xSamples, ySamples)
zzSamples = np.ones(xxSamples.shape)*optFloris.hubHeight[0]
nSamples = len(zzSamples.flatten())

# setup the plot figure
fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize=(20,10))

fig.suptitle('wind direction %.1f deg' % windDirectionsWindRosemaxIncrease)

# .. setup FLORIS for visualization case
yaws = list()

visualFloris = floris_assembly_opt_AEP(nTurbines=nTurbines, nDirections=1, optimize_yaw=False, optimize_position=False, datasize=datasize, nSamples = nSamples, nSpeeds = 1)

visualFloris.windrose_directions    = np.array([windDirection]);
visualFloris.initVelocitiesTurbines = np.copy(np.ones_like(baselineFloris.windrose_directions)*windSpeed)
visualFloris.windrose_speeds        = windSpeed
visualFloris.air_density            = airDensity
visualFloris.curve_wind_speed       = NREL5MWCPCT.wind_speed
visualFloris.curve_CP               = NREL5MWCPCT.CP
visualFloris.curve_CT               = NREL5MWCPCT.CT
visualFloris.rotorDiameter          = np.ones(nTurbines)*rotorDiameter
visualFloris.rotorArea              = np.ones(nTurbines)*rotorArea
visualFloris.turbineX               = turbineX
visualFloris.turbineY               = turbineY
visualFloris.axialInduction         = np.ones(nTurbines)*axialInduction # values used for initialization only
visualFloris.hubHeight              = np.ones(nTurbines)*hub_height
visualFloris.generator_efficiency   = np.ones(nTurbines)*generator_efficiency
visualFloris.ws_positionX           = xxSamples.flatten()
visualFloris.ws_positionY           = yySamples.flatten()
visualFloris.ws_positionZ           = zzSamples.flatten()

# run visual baseline case (default zero yaw)
print 'running visualization'
tic = time.time()
visualFloris.run()
toc = time.time()
velocitiesBaseline = np.copy(visualFloris.ws_array_0.reshape(len(ySamples), len(xSamples)))
yaws.append(np.copy(visualFloris.yaw))

# run visual optimal yawcase
visualFloris.yaw = optYawMaxIncrease
yaws.append(np.copy(visualFloris.yaw))
visualFloris.run()
print 'done'
print('FLORIS visualization took %.03f sec.' % (toc-tic))
velocitiesOpt = np.copy(visualFloris.ws_array_0.reshape(len(ySamples), len(xSamples)))

vmax = np.max([velocitiesBaseline.max(),velocitiesOpt.max()])
vmin = np.min([velocitiesBaseline.min(),velocitiesOpt.min()])
axes[0].pcolormesh(xSamples, ySamples, velocitiesBaseline, cmap='coolwarm', vmin=vmin, vmax=vmax)
axes[0].set_title('baseline')
im = axes[1].pcolormesh(xSamples, ySamples, velocitiesOpt, cmap='coolwarm', vmin=vmin, vmax=vmax)
axes[1].set_title('optimized yaw')

for axI, ax in enumerate(axes):
    ax.set_aspect('equal')
    ax.autoscale(tight=True)
    ax.set_xlabel('Easting (m)')
    ax.set_ylabel('Northing (m)')

    for turbI in range(nTurbines):
        # plot the rotors
        yaw = yaws[axI][turbI]
        turbineX = visualFloris.turbineX[turbI]
        turbineY = visualFloris.turbineY[turbI]
        rotorDiameter = visualFloris.rotorDiameter[turbI]
        
        rotorAbsAngle = (windDirectionFLORISmaxIncrease + yaw)*np.pi/180.
        rotationMatrix = np.array([(np.cos(rotorAbsAngle), -np.sin(rotorAbsAngle)),(np.sin(rotorAbsAngle), np.cos(rotorAbsAngle))])
        rotorX = np.array([0,0])
        rotorY = np.array([-rotorDiameter/2,rotorDiameter/2])
        rotor = np.dot(rotationMatrix,np.array([rotorX,rotorY])) + np.array(([turbineX,turbineX],[turbineY,turbineY]))
        ax.plot(rotor[0,],rotor[1,],'k-')

plt.tight_layout()

if useSubset:
    cax = fig.add_axes([0.85, 0.1, 0.02, 0.8])
    cb = fig.colorbar(im, cax=cax, orientation = 'vertical')
else:
    cax = fig.add_axes([0.3, 0.065, 0.4, 0.02])
    cb = fig.colorbar(im, cax=cax, orientation = 'horizontal')
cb.set_label('Wind Speed (m/s)', fontsize = 10)


if __name__ == "__main__":
    plt.show()
else:
    plt.show(block=False)
