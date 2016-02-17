import matplotlib.pyplot as plt
from utilities.readVTK import *
from utilities.readTurbineArrayPropertiesFAST import readTurbineArrayPropertiesFAST
from utilities.readSuperCONOUT import readSuperCONOUT
import cPickle as pickle
import numpy as np
from scipy.interpolate import griddata, interp1d, interp2d
from scipy import interp

from Parameters import FLORISParameters
from Ellipse_assembly import floris_assembly_opt_AEP
from EEOver.skewedToRotatedEllipse import skewedToRotatedEllipse, rotatedToSkewedEllipse

from Ellipse_assembly import floris_assembly_opt_AEP

# convert to polar coordinates
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return rho, phi

turbineLocation = [1118.1, 1279.5]
rotorDiameter = 126.0
hubHeight = 90.1
windDirection = 30.
windSpeed = 8.0 # m/s

upstreamProfiles = dict()

stabilityCases = ['neutral', 'stable']
casesICOWES3 = ['r01_c1_y30_t-5_R0_I0', 'r03_c3_y30_t-5_R0_I0']
deltaTimeRange = np.arange(705,1000,5) # time range over which to average power ( same as averaging range slides - see ICOWES3_averageProfiles.py )

downStreamDiamRange = np.array([1,3,5])

# linear coefficient: skew c = cCoefficient * distance to rotor (in rotor diameters)
cCoefficientsD = [0.15/5., 1.0/5.]

figInflowProfs, axesInflowProfs = plt.subplots(nrows=1, ncols=3, figsize = (15.,5.))
figFullFlow, axesFullFlow = plt.subplots(nrows=len(stabilityCases), ncols=1, figsize=(15.,10.))
figSlicesSOWFAvsFLORIS, axesSlicesSOWFAvsFLORIS = plt.subplots(nrows=len(downStreamDiamRange), ncols=len(stabilityCases)*2, figsize=(20.,3.*len(downStreamDiamRange)))
figSlicesContour, axesSlicesContour = plt.subplots(nrows=len(downStreamDiamRange), ncols=len(stabilityCases), figsize=(15.,5.*len(downStreamDiamRange)))

colors = ['g','b']

plotPhiRho = False
plotErrors = False

fittedSkewFactors = np.zeros((len(downStreamDiamRange), len(stabilityCases)))
calcSkewFactors = np.zeros((len(downStreamDiamRange), len(stabilityCases)))

nContours = 10

vmin = windSpeed/3.
vmax = windSpeed*1.25

levels = np.linspace(vmin, vmax, nContours)

veerSlope = np.zeros(len(stabilityCases))

powersSOWFA = dict()
powersFLORIS = dict()
velocityHubheight = dict()

nTurbines = 2;
figPowerTime, axesPowerTime = plt.subplots(nrows = nTurbines, ncols = len(stabilityCases), figsize=(15.,10.), sharex=True, sharey=True)

for stabilityCaseI, case in enumerate(stabilityCases):

    caseICOWES3 = casesICOWES3[stabilityCaseI]
    caseFolder = '/scratch/pfleming/runs/ICOWES_3/runCases/' + caseICOWES3

    # load and average turbine power data
    SCOfile = caseFolder + '/CONOUT/superCONOUT.csv'
    SCO = readSuperCONOUT(SCOfile)

    # get average power for each turbine over deltaTimeRange in kW
    powersSOWFA[case] = np.array([np.mean([SCO.data[turb]['Power (W)'][i] for i in range(len(SCO.time)) if (SCO.time[i]>=deltaTimeRange[0] and SCO.time[i] <= deltaTimeRange[-1])]) for turb in range(SCO.nTurbines)])/1000.
    for turb in range(SCO.nTurbines):
        axesPowerTime[turb,stabilityCaseI].plot(SCO.time, SCO.data[turb]['Power (W)']/1000., 'b', label='SOWFA')
        axesPowerTime[turb,stabilityCaseI].plot([deltaTimeRange[0], deltaTimeRange[-1]],[powersSOWFA[case][turb],powersSOWFA[case][turb]],'b--', label='SOWFA Time-Averaged')
        axesPowerTime[turb,stabilityCaseI].set_title('Turbine %d, %s' % ((turb+1), case))
        axesPowerTime[turb,stabilityCaseI].set_ylabel('Generator Power (kW)')
    axesPowerTime[turb,stabilityCaseI].set_xlabel('Time (s)')

    # load vertical slice data 
    sdiFolder = caseFolder + '/sliceDataInstant'
    fileSlice = 'U_slice_2'
    print 'load vertical slice data'
    d = pickle.load(file(sdiFolder + '/' + fileSlice + '.avg_pickle'))
    print 'done'
    cellData = d['cellData']; cellCenters = d['cellCenters'];

    # change reference frame and convert to down-/cross-stream coordinates
    refFrame = np.array([(cellCenters[:,0]-1118.1)/(np.cos(windDirection*np.pi/180.)*rotorDiameter),(cellCenters[:,2]-hubHeight)/rotorDiameter]).transpose()
    rotMat = np.array([[np.cos(-windDirection*np.pi/180.), -np.sin(-windDirection*np.pi/180.), 0.], [np.sin(-windDirection*np.pi/180.), np.cos(-windDirection*np.pi/180.), 0.], [0., 0., 1.]])
    cellDataCrossDown = np.dot(rotMat,cellData.transpose()).transpose()

    def interpolant(samplePoints):
        sampledData = griddata(refFrame, cellDataCrossDown, samplePoints, method='nearest')
        return sampledData

    # get hub-height velocity 5D upstream (used as input to FLORIS to represent free-stream inflow velocity)
    velocityHubheight[case] = interpolant((-5., 0.))[0]

    # get veer profile at 5D in front of rotor
    resolution = 200
    zUpstream = np.linspace(-hubHeight/rotorDiameter, 2., resolution)
    xUpstream = -5.*np.ones(zUpstream.shape)
    velocitiesUpstream = interpolant((xUpstream, zUpstream))
    
    axesInflowProfs[0].plot(zUpstream, velocitiesUpstream[:,0], label=case, color = colors[stabilityCaseI])
    axesInflowProfs[1].plot(zUpstream, velocitiesUpstream[:,1], label=case, color = colors[stabilityCaseI])
    axesInflowProfs[2].plot(zUpstream, velocitiesUpstream[:,2], label=case, color = colors[stabilityCaseI])

    # get linear fit of crosswind component over rotor plane
    zUpstreamRotorPlane = np.linspace(-0.5, 0.5, resolution)
    xUpstreamRotorPlane = -5.*np.ones(zUpstreamRotorPlane.shape)
    velocitiesUpstreamRotorPlane = interpolant((xUpstreamRotorPlane, zUpstreamRotorPlane))

    crosswindLinear = np.polyfit(zUpstreamRotorPlane, velocitiesUpstreamRotorPlane[:,1].flatten(),1)
    velocitiesUpstreamRotorPlaneLin = crosswindLinear[0]*zUpstreamRotorPlane + crosswindLinear[1]
    axesInflowProfs[1].plot(zUpstreamRotorPlane, velocitiesUpstreamRotorPlaneLin, '--', color = colors[stabilityCaseI])

    upstreamProfiles[case] = {'zUpstream': zUpstream, 'velocitiesUpstream': velocitiesUpstream, 'crosswindLinear': crosswindLinear} # save the profiles

    # visualize the flow field SOWFA vertical cut-through
    resolution = 500
    x = np.linspace(-turbineLocation[0]/(np.cos(windDirection*np.pi/180.)*rotorDiameter), 10., resolution)
    z = np.linspace(-hubHeight/rotorDiameter, 2., resolution)
    xMesh, zMesh = np.meshgrid(x, z)

    velocitiesMesh = interpolant((xMesh.flatten(), zMesh.flatten()))
    absVelocitiesMesh = np.sqrt((velocitiesMesh**2).sum(axis=1)).reshape(resolution,resolution)

    if len(stabilityCases) == 1:
        axFullFlow = axesFullFlow
    else:
        axFullFlow = axesFullFlow[stabilityCaseI]

    im = axFullFlow.pcolormesh(x, z, absVelocitiesMesh, cmap='coolwarm')
    axFullFlow.set_aspect('equal')
    axFullFlow.autoscale(tight=True)
    axFullFlow.set_title(case)
    axFullFlow.set_xlabel('distance to hub, in rotor diameters')
    axFullFlow.set_ylabel('distance to hub,\nin rotor diameters')

    axesInflowProfs[0].set_ylabel('$U_x$ (m/s)')
    axesInflowProfs[0].set_xlabel('vertical distance to hub, in rotor diameters')
    axesInflowProfs[1].set_ylabel('$U_y$ (m/s)')
    axesInflowProfs[1].set_xlabel('vertical distance to hub, in rotor diameters')
    axesInflowProfs[1].legend(loc='best')
    axesInflowProfs[2].set_ylabel('$U_z$ (m/s)')
    axesInflowProfs[2].set_xlabel('vertical distance to hub, in rotor diameters')

    figInflowProfs.tight_layout()

    veerSlope[stabilityCaseI] = crosswindLinear[0]

print velocityHubheight

# fit slope of cross-wind speed (veerSlope) to skew factors (cCoefficientsD)
veerToSkewCoefDLinear = np.polyfit(veerSlope, cCoefficientsD,1)

# read in turbine positions
tapfFile = caseFolder + '/constant/turbineArrayPropertiesFAST'
TAPFproperties = readTurbineArrayPropertiesFAST(tapfFile)

# load NREL 5MW characteristics
NREL5MWCPCT = pickle.load(open('NREL5MWCPCT.p'))
datasize = NREL5MWCPCT.CP.size
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0 # used for initialization only
generator_efficiency = 0.944

for stabilityCaseI, case in enumerate(stabilityCases):

    caseICOWES3 = casesICOWES3[stabilityCaseI]

    slices = dict()

    caseICOWES3s = caseICOWES3+'_slices'
        
    slices[case] = dict()

    for downStreamDiamI, downStreamDiam in enumerate(downStreamDiamRange):

        if len(stabilityCases) == 1:
            axContour = axesSlicesContour[downStreamDiamI]
        elif len(downStreamDiamRange) == 1:
            axContour = axesSlicesContour[caseI]
        else:
            axContour = axesSlicesContour[downStreamDiamI, stabilityCaseI]

        if len(downStreamDiamRange) == 1:
            axSlicesSOWFA = axesSlicesSOWFAvsFLORIS[caseI]
            axSlicesFLORIS = axesSlicesSOWFAvsFLORIS[caseI]
        else:
            axSlicesSOWFA = axesSlicesSOWFAvsFLORIS[downStreamDiamI, stabilityCaseI*2]
            axSlicesFLORIS = axesSlicesSOWFAvsFLORIS[downStreamDiamI, stabilityCaseI*2+1]

        print('D ' + str(downStreamDiam))
        slices[case][downStreamDiam] = dict()
        sdiFolder = '/scratch/pfleming/runs/ICOWES_3/runCases/%s/sliceDataInstant' % caseICOWES3s
        fileSlice = 'U_slide_%dDdwnstrT1' % downStreamDiam
        print 'load wake cut-through data'
        d = pickle.load(file(sdiFolder + '/' + fileSlice + '.avg_pickle'))
        print 'done'
        cellData = d['cellData']; cellCenters = d['cellCenters'];
        sdiFolder = '/scratch/pfleming/runs/ICOWES_3/runCases/%s/sliceDataInstant' % caseICOWES3s

        # change reference frame and convert to down-/cross-stream coordinates
        rotMat = np.array([[np.cos(-windDirection*np.pi/180.), -np.sin(-windDirection*np.pi/180.), 0.], [np.sin(-windDirection*np.pi/180.), np.cos(-windDirection*np.pi/180.), 0.], [0., 0., 1.]])
        turbineLocationH = np.append(turbineLocation, hubHeight)
        refFrame = (np.dot(rotMat, (cellCenters-np.kron(np.ones((cellCenters.shape[0],1)), np.array(turbineLocationH))).transpose()).transpose())/rotorDiameter
        refFrame = refFrame[:,1:]
        cellDataCrossDown = np.dot(rotMat,cellData.transpose()).transpose()

        def interpolant(samplePoints):
            sampledData = griddata(refFrame, cellDataCrossDown, samplePoints, method='nearest')
            return sampledData

        # plot cut-through
        resolution = 100
        yCutThrough = np.linspace(-3., 3., resolution)
        zCutThrough = np.linspace(-hubHeight/rotorDiameter, 2., resolution)
        yMesh, zMesh = np.meshgrid(yCutThrough, zCutThrough)
        print 'interpolating'
        velocitiesMesh = interpolant((yMesh.flatten(), zMesh.flatten()))
        print 'done'
        axialVelocitiesMesh = velocitiesMesh[:,0].reshape(resolution,resolution)

        # remove shear profile using upstream shear profile, based on velocity profile -3 diameters left of rotor
        shearProfile = (interpolant((-3.*np.ones(zCutThrough.shape), zCutThrough))[:,0]/interpolant((-5., 0))[0] + interpolant((3.*np.ones(zCutThrough.shape), zCutThrough))[:,0]/interpolant((5., 0))[0])/2

        def interpolantShearProfile(samplePoints):
            sampledData = griddata(zCutThrough, shearProfile, samplePoints, method='linear')
            return sampledData

        shearMesh = interpolantShearProfile(zMesh.flatten()).reshape(resolution,resolution)
        axialVelocitiesMeshShearRemoved = axialVelocitiesMesh/shearMesh

        cs = axContour.contour(yCutThrough, zCutThrough, axialVelocitiesMeshShearRemoved, levels = levels, cmap='coolwarm', vmin=vmin, vmax=vmax)
        axContour.set_aspect('equal')
        axContour.autoscale(tight=True)
        axContour.set_title('%s, %dD downstream' % (case, downStreamDiam))

        pc = axSlicesSOWFA.pcolormesh(yCutThrough, zCutThrough, axialVelocitiesMeshShearRemoved, cmap='coolwarm', vmin=vmin, vmax=vmax)
        axSlicesSOWFA.set_aspect('equal')
        axSlicesSOWFA.autoscale(tight=True)
        axSlicesSOWFA.set_title('SOWFA, %s, %dD downstream' % (case, downStreamDiam))

        slices[case][downStreamDiam] = {'yCutThrough': yCutThrough, 'zCutThrough': zCutThrough, 'axialVelocitiesMeshShearRemoved': axialVelocitiesMeshShearRemoved, 'axialVelocitiesMesh': axialVelocitiesMesh}

        #select_contour = np.argmin(np.abs(cs.levels - (np.min(cs.levels) + (windSpeed-np.min(cs.levels))*0.7)))
        select_contour = np.argmin(np.abs(cs.levels - 7.))
        found = False

        if plotPhiRho:
            figPhiRho, axPhiRho = plt.subplots(1,1,figsize=(20,10))

        R_ellipse_max = 0.

        for contour in cs.collections[select_contour].get_paths():
            path = contour.vertices
            x = path[:,0]
            y = path[:,1]
            
            x_center = np.median(x)
            y_center = np.median(y)
            
            dx = x# - x_center
            dy = y# - y_center
            R_ellipse = np.max(dy)
            print R_ellipse
            rho_old, phi_old = cart2pol(dx, dy)

            phi_old[phi_old<0] = phi_old[phi_old<0]+2*np.pi
            indsort = np.argsort(phi_old)
            phi_old = phi_old[indsort]
            rho_old = rho_old[indsort]

            if np.abs(phi_old[-1] - phi_old[0]) > 0.9*2*np.pi :
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
            print R_ellipse
            rho_old, phi_old = cart2pol(dx, dy)
            phi_old[phi_old<0] = phi_old[phi_old<0]+2*np.pi
            indsort = np.argsort(phi_old)
            phi_old = phi_old[indsort]
            rho_old = rho_old[indsort]

            axContour.plot(x,y,'k--')
            axContour.plot(x_center,y_center,'kx')   

            # resample, use uniform distibution on phi
            fint = interp1d(phi_old, rho_old, kind='linear')
            phi = np.linspace(phi_old[0],phi_old[-1], 1000)
            rho = fint(phi)
            x = np.cos(phi)*rho+x_center
            y = np.sin(phi)*rho+y_center
        else:
            print 'profile not found!'

        def skewedCirclePolar(c, R, phi):
            ws = R*2
            hs = R*2
            r = np.sqrt(hs**2*ws**2/(c**2*hs**2*np.sin(phi)**2 + c*hs**2*np.sin(2*phi) + hs**2*np.cos(phi)**2 + ws**2*np.sin(phi)**2))/2
            return r

        # first use manually fitted skew cCoefficientsD
        c = downStreamDiam * cCoefficientsD[stabilityCaseI];

        r_ellipse = skewedCirclePolar(c, R_ellipse, phi)
        x_ellipse = r_ellipse*np.cos(phi)+x_center
        y_ellipse = r_ellipse*np.sin(phi)+y_center
        axContour.plot(x_ellipse,y_ellipse,'g-')

        # # then use linear fit coefficients
        cCoefficientsDlin = veerToSkewCoefDLinear[0]*veerSlope[stabilityCaseI] + veerToSkewCoefDLinear[1]
        c = downStreamDiam * cCoefficientsDlin;
        r_ellipse = skewedCirclePolar(c, R_ellipse, phi)
        x_ellipse = r_ellipse*np.cos(phi)+x_center
        y_ellipse = r_ellipse*np.sin(phi)+y_center
        axContour.plot(x_ellipse,y_ellipse,'r--')

        # SET UP FLORIS MODEL
        
        # define sampling points
               
        ySampling, zSampling = np.meshgrid(yCutThrough, zCutThrough+hubHeight/rotorDiameter)
        xSampling = np.ones(ySampling.shape)*downStreamDiam;

        xSampling = xSampling * rotorDiameter;
        ySampling = ySampling * rotorDiameter;
        zSampling = zSampling * rotorDiameter;

        rotMat = np.array([[np.cos(windDirection*np.pi/180.), -np.sin(windDirection*np.pi/180.), 0.], [np.sin(windDirection*np.pi/180.), np.cos(windDirection*np.pi/180.), 0.], [0., 0., 1.]])

        samplingPoints = np.dot(rotMat, np.array((xSampling.flatten(),ySampling.flatten(),zSampling.flatten())))

        xSampling = samplingPoints[0] + TAPFproperties['x'][0]
        ySampling = samplingPoints[1] + TAPFproperties['y'][0]
        zSampling = samplingPoints[2]

        nSamples = len(zSampling)

        # define floris model with sampling
        myFloris = floris_assembly_opt_AEP(nTurbines=2, nDirections=1, optimize_yaw=False,
                                   optimize_position=False, use_rotor_components=True,
                                   datasize=datasize, nSamples = nSamples, shearProfileSize = 0)

        myFloris.parameters.CPcorrected = True
        myFloris.parameters.CTcorrected = True
        myFloris.parameters.FLORISoriginal = False

        myFloris.air_density            = 1.1716
        myFloris.axialInduction         = np.array([axialInduction, axialInduction])
        myFloris.rotorDiameter          = np.array([rotorDiameter, rotorDiameter])
        myFloris.rotorArea              = np.array([rotorArea, rotorArea])
        myFloris.hubHeight              = np.array([hubHeight, hubHeight])
        myFloris.curve_CP               = NREL5MWCPCT.CP
        myFloris.curve_CT               = NREL5MWCPCT.CT
        myFloris.curve_wind_speed       = NREL5MWCPCT.wind_speed
        myFloris.generator_efficiency   = np.array([generator_efficiency, generator_efficiency])
        myFloris.turbineX               = TAPFproperties['x']
        myFloris.turbineY               = TAPFproperties['y']
        myFloris.windrose_directions    = np.array([30.])
        myFloris.windrose_speeds        = np.ones_like(myFloris.windrose_directions) * velocityHubheight[case]
        myFloris.parameters.wakeSkewCoefficient = cCoefficientsDlin/rotorDiameter

        # calculate stability correction factors for CP and CT:
        CPinit = interp(velocityHubheight[case],myFloris.curve_wind_speed,myFloris.curve_CP)
        CPtrue = (powersSOWFA[case][0]*1000.)/(myFloris.generator_efficiency[0] * 0.5 * myFloris.air_density * myFloris.rotorArea[0] * velocityHubheight[case]**3.)
        myFloris.curve_CP = myFloris.curve_CP * CPtrue / CPinit
        myFloris.curve_CT = myFloris.curve_CT * CPtrue / CPinit

        myFloris.ws_positionX = xSampling
        myFloris.ws_positionY = ySampling
        myFloris.ws_positionZ = zSampling

        print 'running FLORIS..'
        myFloris.run()
        print 'done'

        velocitiesFLORIS = myFloris.ws_array_0.reshape(len(zCutThrough), len(yCutThrough))
        axSlicesFLORIS.pcolormesh(yCutThrough,zCutThrough,velocitiesFLORIS, cmap='coolwarm',vmin=vmin,vmax=vmax)
        axSlicesFLORIS.set_aspect('equal')
        axSlicesFLORIS.autoscale(tight=True)
        axSlicesFLORIS.set_title('FLORIS, %s, %dD downstream' % (case, downStreamDiam))

    powersFLORIS[case] = myFloris.floris_power_0.wt_power

    for turb in range(nTurbines):
        axesPowerTime[turb,stabilityCaseI].plot([deltaTimeRange[0], deltaTimeRange[-1]],[powersFLORIS[case][turb],powersFLORIS[case][turb]],'r--', label='FLORIS')
        axesPowerTime[turb,stabilityCaseI].legend(loc='best')


cb = figSlicesContour.colorbar(cs, orientation='horizontal', label = 'axial velocity, shear removed (m/s)')
cb.set_clim(vmin, vmax)

cb = figSlicesSOWFAvsFLORIS.colorbar(pc, orientation='horizontal', label = 'axial velocity, shear removed (m/s)')
cb.set_clim(vmin, vmax)

plt.figure()
plt.plot(veerSlope, cCoefficientsD, '*')
plt.plot(veerSlope, veerToSkewCoefDLinear[0]*veerSlope + veerToSkewCoefDLinear[1], '-')


if __name__ == '__main__':
    plt.show()
else:
    plt.show(block=False)
