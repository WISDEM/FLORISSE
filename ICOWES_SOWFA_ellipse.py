import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import numpy as np
from scipy.io import loadmat

import cPickle as pickle

from Parameters import FLORISParameters

from Ellipse_assembly import floris_assembly_opt_AEP

ICOWESdata = loadmat('YawPosResults.mat')
yawrange = ICOWESdata['yaw'][0]

optimize_position = False
optimize_yaw = False
use_rotor_components = True

NREL5MWCPCT = pickle.load(open('NREL5MWCPCT.p'))
datasize = NREL5MWCPCT.CP.size

rotorDiameter = 126.4
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0
CP = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
CT = 4.0*axialInduction*(1.0-axialInduction)
generator_efficiency = 0.944
hub_height = 90.0

resolution = 75

x = np.linspace(800,2200,resolution)
y = np.linspace(1000,2000,resolution)
xx, yy = np.meshgrid(x, y)
ws_positionX = xx.flatten()
ws_positionY = yy.flatten()
ws_positionZ = np.ones(ws_positionX.shape)*hub_height

def tryFLORISyaw(aU,bU,visual):

    if visual:
        myFloris = floris_assembly_opt_AEP(nTurbines=2., nDirections=1, optimize_yaw=optimize_yaw,
                                       optimize_position=optimize_position, use_rotor_components=use_rotor_components,
                                       datasize=datasize, nSamples =len(ws_positionX))
    else:
        myFloris = floris_assembly_opt_AEP(nTurbines=2., nDirections=1, optimize_yaw=optimize_yaw,
                                       optimize_position=optimize_position, use_rotor_components=use_rotor_components,
                                       datasize=datasize)

    myFloris.parameters = FLORISParameters()

    myFloris.parameters.CPcorrected = True
    myFloris.parameters.CTcorrected = True
    myFloris.parameters.FLORISoriginal = False

    # Define turbine characteristics
    myFloris.axialInduction = np.array([axialInduction, axialInduction])
    myFloris.rotorDiameter = np.array([rotorDiameter, rotorDiameter])
    myFloris.rotorArea = np.array([rotorArea, rotorArea])
    myFloris.hubHeight = np.array([hub_height, hub_height])
    myFloris.wakeSkew = np.zeros(2)
    myFloris.generator_efficiency = np.array([generator_efficiency, generator_efficiency])

    # Define site measurements
    myFloris.windrose_directions = np.array([30.])
    wind_speed = 8.1    # m/s
    myFloris.windrose_speeds = np.ones_like(myFloris.windrose_directions)*wind_speed
    myFloris.air_density = 1.1716

    myFloris.initVelocitiesTurbines = np.ones_like(myFloris.windrose_directions)*wind_speed
    myFloris.curve_CP = NREL5MWCPCT.CP
    myFloris.curve_CT = NREL5MWCPCT.CT
    myFloris.curve_wind_speed = NREL5MWCPCT.wind_speed
    myFloris.parameters.ke = 0.051
    myFloris.parameters.kd = 0.17

    myFloris.parameters.aU = aU
    myFloris.parameters.bU = bU
       
    myFloris.parameters.initialWakeAngle = 3.0
    myFloris.parameters.useaUbU = True
    myFloris.parameters.useWakeAngle = True
    myFloris.parameters.adjustInitialWakeDiamToYaw = True

    FLORISpower = list()
    FLORISgradient = list()

    # Defube turbine locations and orientation
    myFloris.turbineX = np.array([1118.1, 1881.9])
    myFloris.turbineY = np.array([1279.5, 1720.5])

    # generate points downstream slices
    y_cut = np.linspace(-rotorDiameter,rotorDiameter,resolution)
    z_cut = np.linspace(-hub_height,rotorDiameter,resolution)
    yy, zz = np.meshgrid(y_cut, z_cut)
    xx = np.ones(yy.shape) * 3.5*rotorDiameter
    position = np.array([xx.flatten(),yy.flatten(),zz.flatten()])
    windDirection = myFloris.windrose_directions[0]*np.pi/180.
    rotationMatrix = np.array([(np.cos(windDirection), -np.sin(windDirection), 0.),
                                       (np.sin(windDirection), np.cos(windDirection),0.),
                                       (0., 0., 1.)])
    positionF = np.dot(rotationMatrix, position) + np.dot(np.array([(myFloris.turbineX[0],myFloris.turbineY[0], hub_height)]).transpose(),np.ones((1,np.size(position,1))))


    velocities = list()
    velocities_cut = list()
    for yaw1 in yawrange:

        myFloris.yaw = np.array([yaw1, 0.0])

        if visual:
            # Call FLORIS horizontal slice
            myFloris.ws_positionX = np.copy(ws_positionX)
            myFloris.ws_positionY = np.copy(ws_positionY)
            myFloris.ws_positionZ = np.copy(ws_positionZ)

        myFloris.run()

        FLORISpower.append(myFloris.floris_power_0.wt_power)

        if visual:
            velocities.append(np.copy(myFloris.ws_array_0))

           # Call FLORIS cut-through slice

            myFloris.ws_positionX = np.copy(positionF[0])
            myFloris.ws_positionY = np.copy(positionF[1])
            myFloris.ws_positionZ = np.copy(positionF[2])

            myFloris.run()

            velocities_cut.append(np.copy(myFloris.ws_array_0))

    if visual:
        velocities = np.array(velocities)
        vmin = np.min(velocities)
        vmax = np.max(velocities)
        velocities_cut = np.array(velocities_cut)

        fig, axes = plt.subplots(ncols=int(np.ceil(len(yawrange)/2.)), nrows=4, figsize=(23,10))

        axes1 = list(axes[0])+list(axes[2])
        axes2 = list(axes[1])+list(axes[3])

        for i in range(len(yawrange)):
            vel = velocities[i].flatten()
            vel = vel.reshape(len(y), len(x))
            ax1 = axes1[i]
            im = ax1.pcolormesh(x, y, vel, cmap='coolwarm', vmin=vmin, vmax=vmax)
            ax1.set_aspect('equal')
            ax1.autoscale(tight=True)

            vel = velocities_cut[i].flatten()
            vel = vel.reshape(len(z_cut), len(y_cut))
            ax2 = axes2[i]
            im = ax2.pcolormesh(y_cut, z_cut, vel, cmap='coolwarm', vmin=vmin, vmax=vmax)
            ax2.set_aspect('equal')
            ax2.autoscale(tight=True)
            ax2.invert_xaxis()

        divider = make_axes_locatable(axes1[-1])
        cax = divider.append_axes("bottom", "5%", pad="20%")
        plt.colorbar(im, cax=cax, orientation = 'horizontal', ticks=[vmin,vmax])
        axes1[-1].axis('off')
        axes2[-1].axis('off')

        mpl.rcParams.update({'font.size': 8})

        plt.tight_layout()

    FLORISpower = np.array(FLORISpower)
    SOWFApower = np.array([ICOWESdata['yawPowerT1'][0],ICOWESdata['yawPowerT2'][0]]).transpose()/1000.

    if visual:
        fig, axes = plt.subplots(ncols = 2, sharey = True)
        axes[0].plot(yawrange.transpose(), FLORISpower[:,0], 'r-', yawrange.transpose(), SOWFApower[:,0], 'ro')
        axes[0].plot(yawrange.transpose(), FLORISpower[:,1], 'b-', yawrange.transpose(), SOWFApower[:,1], 'bo')
        axes[0].plot(yawrange.transpose(), FLORISpower[:,0]+FLORISpower[:,1], 'k-', yawrange.transpose(), SOWFApower[:,0]+SOWFApower[:,1], 'ko')

    error_turbine2 = np.sum(np.abs(FLORISpower[:, 1] - SOWFApower[:, 1]))

    print error_turbine2

    if visual:
        return error_turbine2, myFloris, axes
    else:
        return error_turbine2, myFloris
       
findBestaU = False

if findBestaU:

    aUrange = np.linspace(8.,20.0,20)
    bUrange = np.linspace(0.1,2.0,30)

    errors = np.zeros((len(aUrange),len(bUrange)))

    for aUi, aU in enumerate(aUrange):
        for bUi, bU in enumerate(bUrange):
            errors[aUi,bUi],_ = tryFLORISyaw(aU=aU, bU=bU, visual=False)

    fig = plt.figure()
    plt.pcolormesh(aUrange-(aUrange[1]-aUrange[0])/2.,bUrange-(bUrange[1]-bUrange[0])/2.,errors.transpose())
    plt.colorbar()
    plt.xlabel('aU')
    plt.xlabel('bU')
    [best_aUi,best_bUi] = np.unravel_index(np.argmin(errors), errors.shape)
    aU = aUrange[best_aUi]
    print('best aU %f' % aU)
    bU = bUrange[best_bUi]
    print('best bU %f' % bU)
    plt.plot(aU,bU,'w*')

    error, myFloris, axes = tryFLORISyaw(aU=aU, bU=bU, visual=True)

else:
    aU = 13.052632
    bU = 0.755172
    error, myFloris, axes = tryFLORISyaw(aU=aU, bU=bU, visual=True)

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

axes[1].plot(posrange, FLORISpower[:,0], 'r-', posrange, SOWFApower[:,0], 'ro')
axes[1].plot(posrange, FLORISpower[:,1], 'b-', posrange, SOWFApower[:,1], 'bo')
axes[1].plot(posrange, FLORISpower[:,0]+FLORISpower[:,1], 'k-', posrange, SOWFApower[:,0]+SOWFApower[:,1], 'ko')

plt.show()