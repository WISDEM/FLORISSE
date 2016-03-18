# This example script compares FLORIS predictions with steady-state SOWFA data as obtained 
# throught the simulations described in:
#   

import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt

from scipy.io import loadmat
import pickle

from Parameters import FLORISParameters
from Circle_assembly import floris_assembly_opt_AEP

# select what plots to generate
just_SOWFA = True
plot_prefix = ""

# Load steady-state power data from SOWFA 
ICOWESdata = loadmat('YawPosResults.mat')

# visualization: define resolution
resolution = 0

# Define turbine characteristics
rotorDiameter = 126.4
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0 # used only for initialization
generator_efficiency = 0.944
hub_height = 90.0
NREL5MWCPCT = pickle.load(open('NREL5MWCPCT.p'))
datasize = NREL5MWCPCT.CP.size
turbineXinit = np.array([1118.1, 1881.9])
turbineYinit = np.array([1279.5, 1720.5])
nTurbines = len(turbineXinit)

myFloris = floris_assembly_opt_AEP(nTurbines=2, nDirections=1, optimize_yaw=False,
                                   optimize_position=False,
                                   datasize=datasize, nSamples = resolution*resolution)

# use default FLORIS parameters
myFloris.parameters = FLORISParameters()

# load turbine properties into FLORIS
myFloris.curve_wind_speed = NREL5MWCPCT.wind_speed
myFloris.curve_CP = NREL5MWCPCT.CP
myFloris.curve_CT = NREL5MWCPCT.CT
myFloris.axialInduction = np.array([axialInduction, axialInduction])
myFloris.rotorDiameter = np.array([rotorDiameter, rotorDiameter])
myFloris.rotorArea = np.array([rotorArea, rotorArea])
myFloris.hubHeight = np.array([hub_height, hub_height])
myFloris.generator_efficiency = np.array([generator_efficiency, generator_efficiency])
myFloris.turbineX = turbineXinit
myFloris.turbineY = turbineYinit

# Define site measurements
windDirection = 30.
myFloris.windrose_directions = np.array([windDirection])
wind_speed = 8.1    # m/s
myFloris.windrose_speeds = wind_speed
myFloris.air_density = 1.1716

myFloris.initVelocitiesTurbines = np.ones_like(myFloris.windrose_directions)*wind_speed

# visualization:

# SWEEP TURBINE YAW
FLORISpower = list()
yawrange = ICOWESdata['yaw'][0]

for yaw1 in yawrange:
    print "yaw in: ", yaw1*np.pi/180.
    myFloris.yaw = np.array([yaw1, 0.0])

    myFloris.run()
    FLORISpower.append(myFloris.floris_power_0.wt_power)
# quit()
FLORISpower = np.array(FLORISpower)
SOWFApower = np.array([ICOWESdata['yawPowerT1'][0],ICOWESdata['yawPowerT2'][0]]).transpose()/1000.

figPower, axesPower = plt.subplots(ncols = 2, sharey = True)
axesPower[0].plot(yawrange.transpose(), FLORISpower[:,0], 'r-', yawrange.transpose(), SOWFApower[:,0], 'ro')
axesPower[0].plot(yawrange.transpose(), FLORISpower[:,1], 'b-', yawrange.transpose(), SOWFApower[:,1], 'bo')
axesPower[0].plot(yawrange.transpose(), FLORISpower[:,0]+FLORISpower[:,1], 'k-', yawrange.transpose(), SOWFApower[:,0]+SOWFApower[:,1], 'ko')
axesPower[0].set_xlabel('yaw front turbine 1 (deg)')
axesPower[0].set_ylabel('power (kW)')
axesPower[0].legend(['front turbine FLORIS', 'front turbine SOWFA', 'back turbine FLORIS',  'back turbine SOWFA', 'total FLORIS', 'total SOWFA'])

# SWEEP TURBINE POSITIONS
posrange = ICOWESdata['pos'][0]
myFloris.yaw = np.array([0.0, 0.0])
FLORISpower = list()


for pos2 in posrange:
    
    # Define turbine locations and orientation
    effUdXY = 0.523599
    XY = np.array([turbineXinit, turbineYinit]) + np.dot(np.array([[np.cos(effUdXY),-np.sin(effUdXY)], [np.sin(effUdXY),np.cos(effUdXY)]]), np.array([[0., 0], [0,pos2]]))
    myFloris.turbineX = XY[0,:]
    myFloris.turbineY = XY[1,:]

    
    myFloris.run()
    FLORISpower.append(myFloris.floris_power_0.wt_power)

   

# plot powers
FLORISpower = np.array(FLORISpower)
SOWFApower = np.array([ICOWESdata['posPowerT1'][0],ICOWESdata['posPowerT2'][0]]).transpose()/1000.

axesPower[1].plot(posrange, FLORISpower[:,0], 'r-', posrange, SOWFApower[:,0], 'ro')
axesPower[1].plot(posrange, FLORISpower[:,1], 'b-', posrange, SOWFApower[:,1], 'bo')
axesPower[1].plot(posrange, FLORISpower[:,0]+FLORISpower[:,1], 'k-', posrange, SOWFApower[:,0]+SOWFApower[:,1], 'ko')
axesPower[1].set_xlabel('back turbine displacement (m)')
axesPower[1].set_ylabel('power (kW)')
lgd = axesPower[0].legend(loc='lower left', bbox_to_anchor=(0.0, 1.05),
          fancybox=False, shadow=False, ncol=2)
plt.tight_layout()

if not plot_prefix == "":
    plt.savefig(plot_prefix+"SOWFA.pdf", bbox_extra_artists=(lgd,))#, bbox_inches='tight')
# plt.savefig("masterSowfaFloris.pdf")

# ############
if not just_SOWFA:
    fig, axes = plt.subplots(ncols=2, nrows=1, sharey=False)

    posrange = np.linspace(-3.*rotorDiameter, 30.*rotorDiameter, num=1000)
    yaw = np.array([0.0, 0.0])
    wind_direction = 0.

    myFloris.yaw = yaw
    myFloris.windrose_directions = np.array([wind_direction])
    
    FLORISpower = list()
    FLORISvelocity = list()

    for pos2 in posrange:

        # assign values to yaw
        myFloris.turbineX = np.array([0., pos2])
        myFloris.turbineY = np.array([0.0, 0.0])

        # run the problem at given conditions
        myFloris.run()
#         quit()
        # print np.sqrt((turbineY[0]-turbineY[1])**2+(turbineX[0]-turbineX[1])**2)/rotorDiameter # print downwind distance

        FLORISpower.append(list(myFloris.floris_power_0.wt_power))
        FLORISvelocity.append(list(myFloris.floris_power_0.velocitiesTurbines))

    FLORISpower = np.array(FLORISpower)
    FLORISvelocity = np.array(FLORISvelocity)

    axes[1].plot(posrange/rotorDiameter, FLORISpower[:, 1], '#7CFC00', label='FLORIS model')
    axes[1].plot(np.array([7, 7]), np.array([0, 1800]), '--k', label='Tuning point')
    axes[1].set_xlabel('x/D')
    axes[1].set_ylabel('Power (kW)')
    axes[1].legend(loc=4)

    axes[0].plot(posrange/rotorDiameter, FLORISvelocity[:, 1], '#7CFC00', label='FLORIS model')
    axes[0].plot(np.array([7, 7]), np.array([2, 9]), '--k', label='Tuning point')
    axes[0].set_xlabel('x/D')
    axes[0].set_ylabel('Valocity (m/s)')
    axes[0].legend(loc=4)
    # plt.show()
    if not plot_prefix == "":
        plt.savefig(plot_prefix+"DownwindVelocity.pdf", bbox_extra_artists=(lgd,))#, bbox_inches='tight')
#     plt.savefig("masterPowerVelocityDownwindCorrected.pdf")

    plt.tight_layout()

    FLORIScenters = list()
    FLORISdiameters = list()
    FLORISoverlap = list()
    # prob['wakeCentersYT'][2] = 0.
    for pos2 in posrange:

        # assign values to yaw
        myFloris.turbineX = np.array([0., pos2])
        myFloris.turbineY = np.array([0.0, 0.0])
        
        # run the problem at given conditions
        myFloris.run()

        # print np.sqrt((turbineY[0]-turbineY[1])**2+(turbineX[0]-turbineX[1])**2)/rotorDiameter # print downwind distance
        # print prob['wakeCentersYT']
        print "wakeCentersYT", myFloris.floris_wcent_wdiam_0.wakeCentersYT
        print "wakeDiametersT: ", myFloris.floris_wcent_wdiam_0.wakeDiametersT 
        print "wakeOverlapTRel: ", myFloris.floris_overlap_0.wakeOverlapTRel 
        wakeCentersYT = myFloris.floris_wcent_wdiam_0.wakeCentersYT
        wakeDiametersT = myFloris.floris_wcent_wdiam_0.wakeDiametersT
        wakeOverlapTRel = myFloris.floris_overlap_0.wakeOverlapTRel
        FLORIScenters.append(list(wakeCentersYT))

        FLORISdiameters.append(list(wakeDiametersT))
        FLORISoverlap.append(list(wakeOverlapTRel))
        # print prob['velocitiesTurbines0']

        # print prob['wakeOverlapTRel'][6:]

    FLORIScenters = np.array(FLORIScenters)
    FLORISdiameters = np.array(FLORISdiameters)
    FLORISoverlap = np.array(FLORISoverlap)


    fig, axes = plt.subplots(ncols=2, nrows=2, sharey=False, sharex=False)
    # plot wake center
    axes[0, 0].plot(posrange/rotorDiameter, FLORIScenters[:,2], 'k', label='Wake Center')
    # axes[0, 0].set_xlabel('x/D')
    axes[0, 0].set_ylabel('Position')
    axes[0, 0].legend(loc=1)

    # plot wake diameters
    axes[0, 1].plot(posrange/rotorDiameter, FLORISdiameters[:, 3*nTurbines+0]/rotorDiameter,
                    'b', label='Near Wake')
    axes[0, 1].plot(posrange/rotorDiameter, FLORISdiameters[:, 3*nTurbines+nTurbines]/rotorDiameter,
                    'r', label='Far Wake')
    axes[0, 1].plot(posrange/rotorDiameter, FLORISdiameters[:, 3*nTurbines+2*nTurbines]/rotorDiameter,
                    'y', label='Mixing Zone')
    # axes[0, 1].set_xlabel('x/D')
    axes[0, 1].set_ylabel('Wake Diameter / Rotor Diameter')
    axes[0, 1].legend(loc=2)
    axes[0, 1].set_ylim([-1., 5.])

    # plot wake relative overlap
    axes[1, 0].plot(posrange/rotorDiameter, FLORISoverlap[:, 3*nTurbines+0],
                    'b', label='Near Wake')
    axes[1, 0].plot(posrange/rotorDiameter, FLORISoverlap[:, 3*nTurbines+nTurbines],
                    'r', label='Far Wake')
    axes[1, 0].plot(posrange/rotorDiameter, FLORISoverlap[:, 3*nTurbines+2*nTurbines],
                    'y', label='Mixing Zone')
    axes[1, 0].set_xlabel('x/D')
    axes[1, 0].set_ylabel('Relative Overlap')
    axes[1, 0].legend(loc=0)
    axes[1, 0].set_ylim([-0.1, 1.1])


    posrange = np.linspace(-3.*rotorDiameter, 7.*rotorDiameter, num=300)
    yaw = np.array([0.0, 0.0])
    wind_direction = 0.

    myFloris.yaw = yaw
    myFloris.windrose_directions = np.array([wind_direction])

    FLORISpower = list()

    for pos2 in posrange:

        # assign values to yaw
        myFloris.turbineX = np.array([0., 5*rotorDiameter])
        myFloris.turbineY = np.array([0.0, pos2])

        # run the problem at given conditions
        myFloris.run()

        # print np.sqrt((turbineY[0]-turbineY[1])**2+(turbineX[0]-turbineX[1])**2)/rotorDiameter # print downwind distance

        FLORISpower.append(list(myFloris.floris_power_0.wt_power))

    FLORISpower = np.array(FLORISpower)

    axes[1, 1].plot(posrange/rotorDiameter, FLORISpower[:, 1], '#7CFC00')
    axes[1, 1].set_xlabel('x/D')
    axes[1, 1].set_ylabel('Power (kW)')
    axes[1, 1].legend(loc=4)
    if not plot_prefix == "":
        plt.savefig(plot_prefix+"WakeProfile.pdf", bbox_extra_artists=(lgd,))#, bbox_inches='tight')
#     plt.savefig("masterWakeProfile.pdf")
#####################
plt.show()



# if __name__ == "__main__":
#     plt.show()
# else:
#     plt.show(block=False)
