import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import griddata, interp1d, interp2d
import cPickle as pickle
import pyOpt

from Parameters import FLORISParameters
from Ellipse_assembly import floris_assembly_opt_AEP
from EEOver.skewedToRotatedEllipse import skewedToRotatedEllipse, rotatedToSkewedEllipse

# add shear and plot full cut-through
myFlorisFull = floris_assembly_opt_AEP(nTurbines=1, nDirections=1, optimize_yaw=optimize_yaw,
                           optimize_position=optimize_position, use_rotor_components=use_rotor_components,
                           datasize=datasize, nSamples = len(xx.flatten()), shearProfileSize = len(shearProfileZnorm))
myFlorisFull.shearProfileUnorm = shearProfileUnorm
myFlorisFull.ws_positionX = xx.flatten()
myFlorisFull.ws_positionY = yy.flatten()
myFlorisFull.ws_positionZ = zz.flatten()
myFlorisFull.parameters = myFloris.parameters
myFlorisFull.wakeSkew = myFloris.wakeSkew
myFlorisFull.shearProfileZnorm = myFloris.shearProfileZnorm
myFlorisFull.turbineX = myFloris.turbineX
myFlorisFull.turbineY = myFloris.turbineY
myFlorisFull.windrose_speeds = myFloris.windrose_speeds
myFlorisFull.windrose_directions = myFloris.windrose_directions
myFlorisFull.air_density = myFloris.air_density
myFlorisFull.axialInduction = myFloris.axialInduction
myFlorisFull.rotorDiameter = myFloris.rotorDiameter
myFlorisFull.rotorArea = myFloris.rotorArea
myFlorisFull.hubHeight = myFloris.hubHeight
myFlorisFull.generator_efficiency = myFloris.generator_efficiency
myFlorisFull.initVelocitiesTurbines = myFloris.initVelocitiesTurbines
myFlorisFull.curve_CP = myFloris.curve_CP
myFlorisFull.curve_CT = myFloris.curve_CT
myFlorisFull.curve_wind_speed = myFloris.curve_wind_speed
myFlorisFull.yaw = myFloris.yaw

myFlorisFull.run()

velocities = myFlorisFull.ws_array_0       
velocities = velocities.reshape(xx.shape)

axes = [figFinalFitAxes[downstreamDistI, stabilityCaseI*2], figFinalFitAxes[downstreamDistI, stabilityCaseI*2+1]]
axes[1].contourf(deltay/rotorDiameter,deltaz/rotorDiameter, velocities, 51, cmap='coolwarm', vmin=vmin, vmax=vmax)

axes[1].set_title('FLORIS, %s, at %dD downstream' % (stabilityCasesLabels[stabilityCaseI], downstreamDist))
axes[0].set_title('SOWFA, %s, at %dD downstream' % (stabilityCasesLabels[stabilityCaseI], downstreamDist))

im = axes[0].contourf(deltay/rotorDiameter,deltaz/rotorDiameter,profileSOWFA['uMean'],51, cmap='coolwarm', vmin=vmin, vmax=vmax)

if(downstreamDistI==len(downstreamDistRange)-1 and stabilityCaseI==0):
    axes[0].set_xlabel('distance to hub in rotor diameters')
    axes[0].set_ylabel('distance to hub in rotor diameters')
for ax in axes:
    ax.set_aspect('equal')
    #ax.set_ylim([-100.,100.])
    #ax.set_xlim([-220.,220.])
axes[1].plot((diag_y-turbineY)/rotorDiameter,(diag_z-hub_height)/rotorDiameter,'k--')
axes[0].plot((diag_y-turbineY)/rotorDiameter,(diag_z-hub_height)/rotorDiameter,'k--')


figFinalFit.tight_layout(rect=(0, 0.05, 1, 1))
cax = figFinalFit.add_axes([0.35,0.025,0.3,0.015])
cb = plt.colorbar(im, cax=cax, orientation='horizontal')
cb.set_label('Wind Speed (m/s)')
