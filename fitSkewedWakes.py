import matplotlib as mpl
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp1d, interp2d
import cPickle as pickle
import pyOpt

from Parameters import FLORISParameters
from Ellipse_assembly import floris_assembly_opt_AEP
from EEOver.skewedToRotatedEllipse import skewedToRotatedEllipse, rotatedToSkewedEllipse

optimize_position = False
optimize_yaw = False
use_rotor_components = True

d = loadmat('matWakeAnalysis/turbineData.mat')
datasize = len(d['windSpeedI'].flatten())

rotorDiameter = 77.
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0 # just used for initialization
generator_efficiency = 0.95
hub_height = 80.0
wind_speed = 6.5    # m/s
turbineX = np.array([1000.])
turbineY = np.array([1250.])

downstreamProfiles = True
savefigs = True
showfigs = False 
mpl.rcParams.update({'font.size': 8})
downstreamDistRange = [1,3,5,7,9]

stabilityCases = ['slightlyStable', 'stronglyStable']
stabilityCasesLabels = ['slightly stable', 'strongly stable']

fittedSkewFactors = np.zeros((len(downstreamDistRange), len(stabilityCases)))
calcSkewFactors = np.zeros((len(downstreamDistRange), len(stabilityCases)))

figShear, figShearAxes = plt.subplots(1,len(stabilityCases), sharex=True, sharey=True, figsize=(7,3))
figContourCalculatedSkew, figContourCalculatedSkewAxes = plt.subplots(len(downstreamDistRange), len(stabilityCases), sharex=True, sharey=True, figsize=(15,15))
figDiagWakeCutThrough, figDiagWakeCutThroughAxes = plt.subplots(len(downstreamDistRange), len(stabilityCases), sharex=True, sharey=True, figsize=(10,15))
figFinalFit, figFinalFitAxes = plt.subplots(len(downstreamDistRange), len(stabilityCases)*2, sharex=True, sharey=True, figsize=(20,15))

# Possibly use log-law and correction for 
# https://books.google.com/books?id=jyP2CAAAQBAJ&pg=PA56&lpg=PA56&dq=Parametric+relations+for+the+atmospheric+boundary+layer&source=bl&ots=Y3QXzQPYSm&sig=u7BYsuobOHmatj9zVwl4-CwOR0M&hl=en&sa=X&ved=0CDIQ6AEwAmoVChMIxN60yY-dyQIVSBgeCh3LvAHC#v=onepage&q=Parametric%20relations%20for%20the%20atmospheric%20boundary%20layer&f=false
# http://www.sjsu.edu/people/frank.freedman/courses/metr130/s1/Met130_Spring11_Lect3.pdf

for stabilityCaseI, stabilityCase in enumerate(stabilityCases):
    alphaEstimated = False
    for downstreamDistI, downstreamDist in enumerate(downstreamDistRange):
        profileSOWFA = loadmat('matWakeAnalysis/%s_%dD.mat' % (stabilityCase, downstreamDist))
        vmin = np.min(profileSOWFA['uMean'])
        vmax = np.max(profileSOWFA['uMean'])

        xx = profileSOWFA['x']
        yy = profileSOWFA['y']+turbineY
        zz = profileSOWFA['z']   

        # define shear properties (get predefined full shear profile from data, or use power law)
        shearZh = hub_height
        shearProfileUnorm = 0.5*(profileSOWFA['uMean'][:,0]+profileSOWFA['uMean'][:,-1])/wind_speed
        shearProfileZnorm = zz[:,0]/shearZh
        usePredefinedShearProfile = True

        if downstreamDistI == 0:
            figShearAxes[stabilityCaseI].plot(shearProfileZnorm, shearProfileUnorm, 'r')
            figShearAxes[stabilityCaseI].set_xlabel('z / hub height')
            figShearAxes[stabilityCaseI].set_ylabel('average axial wind speed / hub-height wind speed')
            figShearAxes[stabilityCaseI].set_title(stabilityCasesLabels[stabilityCaseI])

        if not usePredefinedShearProfile:
            # estimate shear coefficient of power law
            # only use first profile, then check results in next 
            if not alphaEstimated:
                alphaEstimated = True
                # 1: use top velocity for first estimate of alpha
                top = np.argmax(shearProfileZnorm)
                z_top = shearProfileZnorm[top]
                UtopNorm = shearProfileUnorm[top]
                
                shearCoefficientAlpha = np.log(UtopNorm)/np.log(z_top)
                # 2: use shear profile from SOWFA data (furthest away from hub) and search alpha with least error
                alphaRange = np.linspace(0.1*shearCoefficientAlpha,10*shearCoefficientAlpha,200)
                RMSerrors = np.zeros(len(alphaRange))
                for alphaI, alpha in enumerate(alphaRange):
                    shearProfileLogTry = shearProfileZnorm**alpha
                    RMSerrors[alphaI] = np.sqrt(((shearProfileUnorm - shearProfileLogTry) ** 2).mean())
                alphaI = np.argmin(RMSerrors)
                shearCoefficientAlpha = alphaRange[alphaI]
                shearProfileLog = shearProfileZnorm**shearCoefficientAlpha

            # 3: plot profile
            plt.plot(shearProfileZnorm, shearProfileLog, 'b--')
            plt.title('shear coefficient alpha %f' % shearCoefficientAlpha)

        # for wake detection, remove shear profile from data
        profileSOWFAnoShear = profileSOWFA['uMean']/(np.tile(shearProfileUnorm,(np.size(profileSOWFA['uMean'],1),1)).transpose())
        fig = plt.figure()
        im = plt.contourf(yy-turbineY,zz-hub_height,profileSOWFAnoShear,51, cmap='coolwarm', vmin=vmin, vmax=vmax)
        cax = fig.add_axes([0.2,0.95,0.6,0.02])
        cb = plt.colorbar(im, cax=cax, orientation='horizontal')
        
        # convert to polar coordinates
        def cart2pol(x, y):
            rho = np.sqrt(x**2 + y**2)
            phi = np.arctan2(y, x)
            return rho, phi

        nContours = 20
        # find contours
        figContour, axContour = plt.subplots(1,1,figsize=(20,10))
        cs = plt.contour(yy-turbineY,zz-hub_height,profileSOWFAnoShear,nContours, cmap='coolwarm', vmin=vmin, vmax=vmax)
        cb = plt.colorbar(cs, shrink=0.8, extend='both')
        select_contour = np.argmin(np.abs(cs.levels - (np.min(cs.levels) + (wind_speed-np.min(cs.levels))*0.4)))
        plt.axis('equal')

        figPhiRho, axPhiRho = plt.subplots(1,1,figsize=(20,10))

        R_ellipse_max = 0.
        found = False

        print downstreamDist
        print cs.levels[select_contour]

        for contour in cs.collections[select_contour].get_paths():
         
            path = contour.vertices
            x = path[:,0]
            y = path[:,1]
            
            x_center = np.median(x)
            y_center = np.median(y)
            
            dx = x# - x_center
            dy = y# - y_center
            R_ellipse = np.max(dy)
            rho_old, phi_old = cart2pol(dx, dy)

            phi_old[phi_old<0] = phi_old[phi_old<0]+2*np.pi
            indsort = np.argsort(phi_old)
            phi_old = phi_old[indsort]
            rho_old = rho_old[indsort]

            if np.abs(phi_old[-1] - phi_old[0]) > 0.9*2*np.pi :
                axPhiRho.plot(phi_old,rho_old,'.-')
                if R_ellipse > R_ellipse_max:
                    R_ellipse_max = R_ellipse
                    found = True
                    x_select = x
                    y_select = y

        if found:
            print 'profile found'
            x = x_select
            y = y_select
            x_center = np.median(x)
            y_center = np.median(y)
            dx = x - x_center
            dy = y - y_center
            R_ellipse = np.max(dy)
            rho_old, phi_old = cart2pol(dx, dy)
            phi_old[phi_old<0] = phi_old[phi_old<0]+2*np.pi
            indsort = np.argsort(phi_old)
            phi_old = phi_old[indsort]
            rho_old = rho_old[indsort]

            axContour.plot(x,y,'k-')
            axContour.plot(x_center,y_center,'kx')   

            # resample, use uniform distibution on phi
            fint = interp1d(phi_old, rho_old, kind='linear')
            phi = np.linspace(phi_old[0],phi_old[-1], 1000)
            rho = fint(phi)
            x = np.cos(phi)*rho+x_center
            y = np.sin(phi)*rho+y_center

        else:
            print 'profile not found!'

        plt.plot(x,y,'b-')
        plt.plot(x_center,y_center,'bx')
        if savefigs:
            figContour.savefig('matWakeAnalysis/fitResults/%s_D%d_figContour.png' % (stabilityCase,downstreamDist))

        max_rho_i = np.argmax(rho)
        max_rho = rho[max_rho_i]
        max_rho_phi = phi[max_rho_i]

        def skewedCirclePolar(c, R, phi):
            ws = R*2
            hs = R*2
            r = np.sqrt(hs**2*ws**2/(c**2*hs**2*np.sin(phi)**2 + c*hs**2*np.sin(2*phi) + hs**2*np.cos(phi)**2 + ws**2*np.sin(phi)**2))/2
            return r


        # find left upper point in contour
        contour_left_upper_i = np.argmax(rho[phi<np.pi])

        # find skew from left upper point
        phi_ellipseLU = phi[contour_left_upper_i]
        skew_ellipseLU = rotatedToSkewedEllipse(R_ellipse*2, R_ellipse*2, phi_ellipseLU)

        # show other ellipses with similar skew
        resolution = 9
        c_range = np.linspace(skew_ellipseLU-0.5,skew_ellipseLU+0.5,resolution-1)
        errors = np.zeros(len(c_range))

        figContoursSkew,axes = plt.subplots(3, resolution/3, figsize=(20,10))
        axes = axes.flatten()

        figSkewRhoErrors,axes2 = plt.subplots(3, resolution/3, sharex=True, sharey=True, figsize=(20,10))
        axes2 = axes2.flatten()

        for ci,c in enumerate(c_range):
            r_ellipse = skewedCirclePolar(c, R_ellipse, phi)
            x_ellipse = r_ellipse*np.cos(phi)+x_center
            y_ellipse = r_ellipse*np.sin(phi)+y_center
            axes[ci].contour(yy-turbineY,zz-hub_height,profileSOWFAnoShear,nContours, cmap='coolwarm', vmin=vmin, vmax=vmax)
            axes[ci].set_aspect('equal')

            # # uses error at extreme points of fitted ellipse
            rotatedEllipse_phi, rotatedEllipse_w, rotatedEllipse_h = skewedToRotatedEllipse(R_ellipse*2, R_ellipse*2, c)
            i_leftUpper = np.argmin(np.abs(phi-(rotatedEllipse_phi+np.pi/2)))
            i_rightLower = np.argmin(np.abs(phi-(rotatedEllipse_phi+np.pi/2+np.pi)))
            errors[ci] = np.sqrt((x_ellipse[i_leftUpper]-x[contour_left_upper_i])**2 + (y_ellipse[i_leftUpper]-y[contour_left_upper_i])**2)# + np.sqrt((x_ellipse[i_rightLower]-x[i_rightLower])**2 + (y_ellipse[i_rightLower]-y[i_rightLower])**2)
            # errors[ci] = np.sqrt(np.mean((((r_ellipse-rho)/rho)**2))) # RMS r
            axes2[ci].plot((((r_ellipse-rho)/rho)**2))
            axes[ci].plot(x_center,y_center,'kx')
            axes[ci].plot(x,y,'k')
            axes[ci].plot(x_ellipse,y_ellipse,'b-')

            axes[ci].plot(x_ellipse[i_leftUpper],y_ellipse[i_leftUpper],'bx')
            axes[ci].plot(x[contour_left_upper_i],y[contour_left_upper_i],'kx')

            axes[ci].set_title('skew = %.2f, rotation = %.2f' % (c, (rotatedEllipse_phi*180./np.pi)))

        best_ci = np.argmin(errors)
        best_c = c_range[best_ci]

        best_R = R_ellipse

        r_ellipse = skewedCirclePolar(best_c, best_R, phi)
        x_ellipse = r_ellipse*np.cos(phi)+x_center
        y_ellipse = r_ellipse*np.sin(phi)+y_center
        axes[best_ci].plot(x_ellipse,y_ellipse,'g-')

        r_ellipseLU = skewedCirclePolar(skew_ellipseLU, R_ellipse, phi)
        x_ellipseLU = r_ellipseLU*np.cos(phi)+x_center
        y_ellipseLU = r_ellipseLU*np.sin(phi)+y_center

        axes[-1].contour(yy-turbineY,zz-hub_height,profileSOWFAnoShear,nContours, cmap='coolwarm', vmin=vmin, vmax=vmax)
        axes[-1].set_aspect('equal')
        axes[-1].plot(x_ellipseLU,y_ellipseLU,'g--')
        axes[-1].plot(x_center,y_center,'kx')
        axes[-1].plot(x,y,'k-')
        axes[-1].set_title('skew = %.2f, rotation = %.2f' % (skew_ellipseLU, (phi_ellipseLU*180./np.pi)))
        axes2[-1].plot((((r_ellipseLU-rho)/rho)**2))

        figContourCalculatedSkewAxes[downstreamDistI, stabilityCaseI].contour(yy-turbineY,zz-hub_height,profileSOWFAnoShear,nContours, cmap='coolwarm', vmin=vmin, vmax=vmax)
        figContourCalculatedSkewAxes[downstreamDistI, stabilityCaseI].set_aspect('equal')
        figContourCalculatedSkewAxes[downstreamDistI, stabilityCaseI].plot(x_ellipseLU,y_ellipseLU,'g--')
        figContourCalculatedSkewAxes[downstreamDistI, stabilityCaseI].plot(x_center,y_center,'kx')
        figContourCalculatedSkewAxes[downstreamDistI, stabilityCaseI].plot(x,y,'k-')
        figContourCalculatedSkewAxes[downstreamDistI, stabilityCaseI].set_title('%s, at %dD downstr., ellipse skew=%.2f, rotation=%.2f deg' % (stabilityCasesLabels[stabilityCaseI], downstreamDistRange[downstreamDistI], skew_ellipseLU, (phi_ellipseLU*180./np.pi)))

        fittedSkewFactors[downstreamDistI, stabilityCaseI] = best_c
        calcSkewFactors[downstreamDistI, stabilityCaseI] = skew_ellipseLU

        if savefigs:
            figContoursSkew.savefig('matWakeAnalysis/fitResults/%s_D%d_figContoursSkew.png' % (stabilityCase,downstreamDist))
            figSkewRhoErrors.savefig('matWakeAnalysis/fitResults/%s_D%d_figSkewRhoErrors.png' % (stabilityCase,downstreamDist))    

        figSkewError = plt.figure(figsize=(20,10))
        plt.plot(c_range,errors)
        plt.xlabel('skew')
        plt.ylabel('error')
        if savefigs:
            figSkewError.savefig('matWakeAnalysis/fitResults/%s_D%d_figSkewError.png' % (stabilityCase,downstreamDist))

        deltay = yy-turbineY
        deltaz = zz-hub_height


        # for getting SOWFA diagonal profile along skew angle without shear
        diag_z = np.linspace(np.minimum(1.5*R_ellipse, np.max(deltaz)), np.maximum(-1.5*R_ellipse, np.min(deltaz)), 25)
        diag_y = diag_z/np.tan(phi_ellipseLU)
        diag_x = np.ones(len(diag_z)) * downstreamDist*rotorDiameter

        fUsowfa = interp2d(deltay[0,:].flatten(), deltaz[:,0].flatten(), profileSOWFAnoShear)
        diagU_SOWFA = np.array([fUsowfa(diag_y[i], diag_z[i]) for i in range(len(diag_y))]).flatten()

        diag_x = np.ones(len(diag_x))*turbineX + diag_x
        diag_y = np.ones(len(diag_y))*turbineY + diag_y
        diag_z = np.ones(len(diag_z))*hub_height + diag_z

        # set up floris to calculate diagonal
        myFloris = floris_assembly_opt_AEP(nTurbines=1, nDirections=1, optimize_yaw=optimize_yaw,
                                           optimize_position=optimize_position, use_rotor_components=use_rotor_components,
                                           datasize=datasize, nSamples = diag_z.size, shearProfileSize = len(shearProfileZnorm))

        myFloris.parameters = FLORISParameters()

        myFloris.parameters.CPcorrected = True
        myFloris.parameters.CTcorrected = True
        myFloris.parameters.FLORISoriginal = False

        # Define turbine characteristics
        myFloris.axialInduction = np.array([axialInduction])
        myFloris.rotorDiameter = np.array([rotorDiameter])
        myFloris.rotorArea = np.array([rotorArea])
        myFloris.hubHeight = np.array([hub_height])
        myFloris.generator_efficiency = np.array([generator_efficiency])

        # Define site measurements
        myFloris.windrose_directions = np.array([0.])
        myFloris.windrose_speeds = np.ones_like(myFloris.windrose_directions)*wind_speed
        myFloris.air_density = 1.225
       
        myFloris.initVelocitiesTurbines = np.ones_like(myFloris.windrose_directions)*wind_speed
        myFloris.curve_CP = d['rotCpI'].flatten()
        myFloris.curve_CT = d['rotCtI'].flatten()
        myFloris.curve_wind_speed = d['windSpeedI'].flatten()
        myFloris.parameters.ke = 0.051
        myFloris.parameters.kd = 0.17
        myFloris.parameters.aU = 13.052632
        myFloris.parameters.bU = 0.755172

        myFloris.parameters.me = np.array([-0.5, 0.25, 1.0])
        myFloris.parameters.MU = np.array([  0.7,   1. ,  5.5 ])
           
        myFloris.parameters.initialWakeAngle = 1.0 # 3.0
        myFloris.parameters.useaUbU = True
        myFloris.parameters.useWakeAngle = True
        myFloris.parameters.adjustInitialWakeDiamToYaw = True

        # Defube turbine locations and orientation
        myFloris.turbineX = turbineX
        myFloris.turbineY = turbineY

        myFloris.yaw = np.array([0.0])

        myFloris.ws_positionX = diag_x
        myFloris.ws_positionY = diag_y
        myFloris.ws_positionZ = diag_z

        myFloris.parameters.usePredefinedShearProfile = usePredefinedShearProfile

        # run model with fitted skew, but without shear
        myFloris.wakeSkew = np.array([skew_ellipseLU])
        myFloris.shearProfileZnorm = shearProfileZnorm
        myFloris.shearProfileUnorm = np.ones(len(shearProfileZnorm))
        print usePredefinedShearProfile
        if not usePredefinedShearProfile:
            myFloris.parameters.shearCoefficientAlpha = shearCoefficientAlpha
        myFloris.parameters.shearZh = shearZh

        myFloris.parameters.initialWakeAngle = 5.0
        myFloris.parameters.me=np.array([-0.5, 0.25, 1.0])
        
        if stabilityCase == 'slightlyStable':
            myFloris.parameters.MU = np.array([1.1, 1.4, 5.5])
            #myFloris.parameters.MU = np.array([0.7, 1.0, 5.5])
        elif stabilityCase == 'stronglyStable':
            myFloris.parameters.MU = np.array([1.6, 1.8, 5.5])
            #myFloris.parameters.MU = np.array([0.7, 1.0, 5.5])

        myFloris.run()
        diagU_FLORIS = myFloris.ws_array_0

        figDiagWakeCutThroughAxes[downstreamDistI, stabilityCaseI].plot(diag_y, diagU_SOWFA)
        figDiagWakeCutThroughAxes[downstreamDistI, stabilityCaseI].plot(diag_y, diagU_FLORIS)

        

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


figSkewVsDist, axfigSkewVsDist = plt.subplots(1,1)
for stabilityCaseI, stabilityCase in enumerate(stabilityCases):
    #axfigSkewVsDist.plot(downstreamDistRange, fittedSkewFactors[:,stabilityCaseI], 'g.-', label='fitted skew %s' % stabilityCasesLabels[stabilityCaseI], linewidth=(stabilityCaseI+1))
    axfigSkewVsDist.plot(downstreamDistRange, calcSkewFactors[:,stabilityCaseI], 'g.--', label=stabilityCasesLabels[stabilityCaseI], linewidth=(stabilityCaseI+1))
    axfigSkewVsDist.set_xlabel('distance downstream from rotor in rotor diameters')
    axfigSkewVsDist.set_ylabel('estimated skew factor')
axfigSkewVsDist.legend(loc='best')

if savefigs:
    figShear.savefig('matWakeAnalysis/fitResults/figShear')
    figContourCalculatedSkew.savefig('matWakeAnalysis/fitResults/figContourCalculatedSkew')
    figDiagWakeCutThrough.savefig('matWakeAnalysis/fitResults/figDiagWakeCutThrough')
    figFinalFit.savefig('matWakeAnalysis/fitResults/figFinalFit')
    figSkewVsDist.savefig('matWakeAnalysis/fitResults/figSkewVsDist.png')
else:
    print 'figures saved in matWakeAnalysis/fitResults'

if showfigs:
    if __name__ == "__main__":
        plt.show()
    else:
        plt.show(block=False)
