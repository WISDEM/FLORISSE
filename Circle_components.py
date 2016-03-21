from openmdao.main.api import Component, VariableTree
from openmdao.lib.datatypes.api import Array, Bool, Float, VarTree
from Parameters import FLORISParameters
import numpy as np


class floris_windframe(Component):
    """ Calculates the locations of each turbine in the wind direction reference frame """

    # original variables
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')
    # position = Array(iotype='in', units='m', desc='position of turbines in original ref. frame')

    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    wind_direction = Float(iotype='in', units='deg', desc='overall wind direction for wind farm')

    def __init__(self, nTurbines, nSamples=0):

        super(floris_windframe, self).__init__()

        # Explicitly size input arrays
        self.add('turbineX', Array(np.zeros(nTurbines), iotype='in', \
                                   desc='x positions of turbines in original ref. frame'))
        self.add('turbineY', Array(np.zeros(nTurbines), iotype='in', \
                                   desc='y positions of turbines in original ref. frame'))

        # variables for verbosity
        self.add('Ct', Array(np.zeros(nTurbines), iotype='in'))
        self.add('Cp', Array(np.zeros(nTurbines), iotype='in', \
                             desc='power coefficient for all turbines'))
        self.add('axialInduction', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                                         desc='axial induction of all turbines'))
        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', \
                              desc='yaw of each turbine'))

        # for testing purposes only
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='out', units='m', \
                                    desc='x coordinates of turbines in wind dir. ref. frame'))
        self.add('turbineYw', Array(np.zeros(nTurbines), iotype='out', units='m', \
                                    desc='y coordinates of turbines in wind dir. ref. frame'))

        # variables for testing wind speed at various locations
        self.add('ws_positionX', Array(np.zeros([nSamples]), iotype='in', units='m', desc='X position of desired measurements in original ref. frame'))
        self.add('ws_positionY', Array(np.zeros([nSamples]), iotype='in', units='m', desc='Y position of desired measurements in original ref. frame'))
        self.add('ws_positionZ', Array(np.zeros([nSamples]), iotype='in', units='m', desc='Z position of desired measurements in original ref. frame'))

        # Explicitly size output arrays
        self.add('wsw_position', Array(np.zeros([3, nSamples]), iotype='out', units='m', desc='position of desired measurements in wind ref. frame'))

    def execute(self):

        Vinf = self.wind_speed
        windDirection = self.wind_direction*np.pi/180.0
        
        #variables to satisfy verbosity
        axialInd = self.axialInduction
        Cp = self.Cp
        Ct = self.Ct
        CTcorrected = self.parameters.CTcorrected
        # yaw = self.yaw*np.pi/180

        # get rotor coefficients, and apply corrections if necesary
        # Cp = np.hstack(self.wt_layout.wt_array(attr='CP'))
        # if CTcorrected == False:
        #     Ct = Ct * (np.cos(yaw)**2)


        if self.verbose:
            np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print "wind direction %s deg" % [windDirection*180.0/np.pi]
            print "free-stream wind speed %s" % Vinf
            print "axial induction turbines %s" % axialInd
            print "C_P turbines %s" % Cp
            print "C_T turbines %s" % Ct
            # print "yaw turbines %s" % yaw

        # get turbine positions and velocity sampling positions
        turbineX = self.turbineX
        turbineY = self.turbineY

        if len(self.ws_positionX)>0:
            velX = self.ws_positionX
            velY = self.ws_positionY
            velZ = self.ws_positionZ
        else:
            velX = np.zeros([0, 0])
            velY = np.zeros([0, 0])
            velZ = np.zeros([0, 0])

        # convert to downwind-crosswind coordinates
        rotationMatrix = np.array([(np.cos(-windDirection), -np.sin(-windDirection)),
                                   (np.sin(-windDirection), np.cos(-windDirection))])
        turbineLocations = np.dot(rotationMatrix, np.array([turbineX, turbineY]))
        # print turbineLocations
        self.turbineXw = np.zeros(turbineX.size)
        self.turbineYw = np.zeros(turbineX.size)
        self.turbineXw = turbineLocations[0]
        self.turbineYw = turbineLocations[1]

        if velX.size>0:
            locations = np.dot(rotationMatrix,np.array([velX,velY]))
            velX = locations[0]
            velY = locations[1]

        self.wsw_position = np.array([velX, velY, velZ])


class floris_wcent_wdiam(Component):
    """ Calculates the center and diameter of each turbine wake at each other turbine """

    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    def __init__(self, nTurbines, nSamples=0):
        super(floris_wcent_wdiam, self).__init__()

        # Explicitly size input arrays
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='in', \
                                    desc='x coordinates of turbines in wind dir. ref. frame'))
        self.add('turbineYw', Array(np.zeros(nTurbines), iotype='in', \
                                    desc='y coordinates of turbines in wind dir. ref. frame'))
        self.add('yaw', Array(np.zeros(nTurbines), iotype='in', desc='yaw of each turbine'))
        self.add('rotorDiameter', Array(np.zeros(nTurbines), dtype='float', iotype='in', \
                                        desc='rotor diameter of each turbine'))
        self.add('hubHeight', Array(np.zeros(nTurbines), dtype='float', iotype='in', units='m', \
                desc='hub heights of all turbines'))
        self.add('Ct', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                             desc='thrust coefficient of each turbine'))

        self.add('wsw_position', Array(np.zeros([3, nSamples]), iotype='in', units='m', desc='positions where measurements are desired in the windframe'))


        # Explicitly size output arrays
        self.add('wakeCentersYT', Array(np.zeros(nTurbines*nTurbines), iotype='out', dtype='float', \
                                        desc='wake center y position at each turbine'))
        self.add('wakeDiametersT', Array(np.zeros(nTurbines*nTurbines*3), iotype='out', dtype='float', \
                                         desc='wake diameter of each zone of each wake at each turbine'))
        self.add('wakeDiameters', Array(np.zeros([nSamples, nTurbines, 3]), iotype='out', dtype='float', desc='wake diameter of each zone of each wake at each turbine'))
        self.add('wakeCentersY', Array(np.zeros([nSamples, nTurbines]), iotype='out', units='m', desc='Y positions of wakes at measurement points'))
        self.add('wakeCentersZ', Array(np.zeros([nSamples, nTurbines]), iotype='out', units='m', desc='Z positions of wakes at measurement points'))

    def execute(self):

         # rename inputs and outputs
        # pP = self.parameters.pP
        kd = self.parameters.kd
        ke = self.parameters.ke
        initialWakeDisplacement = self.parameters.initialWakeDisplacement
        initialWakeAngle = self.parameters.initialWakeAngle
        rotorDiameter = self.rotorDiameter
        hubHeight = self.hubHeight
        Ct = self.Ct
        CTcorrected = self.parameters.CTcorrected
        keCorrCT = self.parameters.keCorrCT
        baselineCT = self.parameters.baselineCT
        me = self.parameters.me

        adjustInitialWakeDiamToYaw = self.parameters.adjustInitialWakeDiamToYaw
        useWakeAngle = self.parameters.useWakeAngle

        bd = self.parameters.bd

        turbineXw = self.turbineXw
        turbineYw = self.turbineYw
        yaw = self.yaw*np.pi/180.0
        nTurbines = turbineXw.size

        velX = self.wsw_position[0][:]
        nSamples = np.size(velX)

        if CTcorrected == False:
            Ct = Ct * (np.cos(yaw*np.pi/180.)**2)

        # calculate y-location of wake centers
        wakeCentersY = np.zeros((nSamples, nTurbines))
        wakeCentersZ = np.zeros((nSamples, nTurbines))
        wakeCentersYT_mat = np.zeros((nTurbines, nTurbines))
        for turb in range(0, nTurbines):
            wakeAngleInit = 0.5 * np.sin(yaw[turb]) * Ct[turb]
            if useWakeAngle:
                wakeAngleInit += initialWakeAngle*np.pi/180.0
            for loc in range(0, nSamples):  # at velX-locations
                deltax = np.maximum(velX[loc]-turbineXw[turb], 0)
                factor = (2.0*kd*deltax/rotorDiameter[turb])+1.0
                wakeCentersY[loc, turb] = turbineYw[turb]+initialWakeDisplacement # initial displacement for no yaw (positive to the left looking downstream)
                displacement = (wakeAngleInit*(15.0*(factor**4.0)+(wakeAngleInit**2.0))/((30.0*kd*(factor**5.0))/rotorDiameter[turb]))-(wakeAngleInit*rotorDiameter[turb]*(15.0+(wakeAngleInit**2.0))/(30.0*kd)) # yaw-induced deflection
                if not useWakeAngle:
                    displacement += bd*deltax
                wakeCentersY[loc, turb] = wakeCentersY[loc, turb] + displacement
                wakeCentersZ[loc, turb] = hubHeight[turb]

            for turbI in range(0, nTurbines):  # at turbineX-locations
                deltax = np.maximum(turbineXw[turbI]-turbineXw[turb], 0.0)
                factor = (2.0*kd*deltax/rotorDiameter[turb])+1.0
                wakeCentersYT_mat[turbI, turb] = turbineYw[turb]
                wakeCentersYT_mat[turbI, turb] = wakeCentersYT_mat[turbI, turb]+initialWakeDisplacement # initial displacement for no yaw (positive to the left looking downstream)
                                
                displacement = (wakeAngleInit*(15.0*(factor**4.0)+(wakeAngleInit**2.0))/
                ((30.0*kd*(factor**5.0))/rotorDiameter[turb]))- \
                (wakeAngleInit*rotorDiameter[turb]*(15.0+(wakeAngleInit**2.0))/(30.0*kd)) # yaw-induced wake center displacement
                
                wakeCentersYT_mat[turbI, turb] = wakeCentersYT_mat[turbI, turb] + displacement
                
        # adjust k_e to C_T, adjusted to yaw
        ke = ke + keCorrCT*(Ct-baselineCT) # FT = Ct*0.5*rho*A*(U*cos(yaw))^2, hence, thrust decreases with cos^2
                                                           #   Should ke increase directly with thrust? ==>No - Turbulence characteristics in wind-turbine wakes, A. Crespo"'*, J. Hern'andez b

        # calculate wake zone diameters at velX-locations
        wakeDiameters = np.zeros((nSamples,nTurbines, 3))
        wakeDiametersT_mat = np.zeros((nTurbines, nTurbines, 3))
        for turb in range(0, nTurbines):
            if adjustInitialWakeDiamToYaw:
                wakeDiameter0 = rotorDiameter[turb] * np.cos(yaw[turb]) # CHANGE: initial wake diameter at rotor adjusted to yaw
            else:
                wakeDiameter0 = rotorDiameter[turb]
            for loc in range(0, nSamples):  # at velX-locations
                deltax = velX[loc]-turbineXw[turb]
                for zone in range(0,3):
                    wakeDiameters[loc, turb, zone] = wakeDiameter0 + 2*ke[turb]*me[zone]*np.maximum(deltax, 0)
            for turbI in range(0, nTurbines):  # at turbineX-locations
                deltax = turbineXw[turbI]-turbineXw[turb]
                for zone in range(0, 3):
                    wakeDiametersT_mat[turbI, turb, zone] = np.maximum(wakeDiameter0 + 2*ke[turb]*me[zone]*deltax, 0)


        wakeDiametersT_vec = np.zeros(3*nTurbines*nTurbines)
        for i in range(0, nTurbines):
            wakeDiametersT_vec[(3*nTurbines*i):(3*nTurbines*i+nTurbines)] = wakeDiametersT_mat[i, :, 0]
            wakeDiametersT_vec[3*nTurbines*i+nTurbines:3*nTurbines*i+2*nTurbines] = wakeDiametersT_mat[i, :, 1]
            wakeDiametersT_vec[3*nTurbines*i+2*nTurbines:3*nTurbines*i+3*nTurbines] = wakeDiametersT_mat[i, :, 2]

        wakeCentersYT_vec = np.zeros(nTurbines*nTurbines)
        for i in range(0, nTurbines):
            wakeCentersYT_vec[(nTurbines*i):(nTurbines*i+nTurbines)] = wakeCentersYT_mat[i, :]


        self.wakeCentersYT = wakeCentersYT_vec
        self.wakeDiametersT = wakeDiametersT_vec
        self.wakeDiameters = wakeDiameters
        self.wakeCentersY = wakeCentersY
        self.wakeCentersZ = wakeCentersZ


class floris_overlap(Component):
    """ Calculates the overlap between each turbine rotor and the existing turbine wakes """

    def __init__(self, nTurbines):
        super(floris_overlap, self).__init__()

        # Explicitly size input arrays
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='in', units='m', ignore_deriv=True, \
                                    desc='X positions of turbines wrt the wind direction'))
        self.add('turbineYw', Array(np.zeros(nTurbines), iotype='in', units='m', \
                                          desc='Y positions of turbines wrt the wind direction'))
        self.add('rotorDiameter', Array(np.zeros(nTurbines), iotype='in', units='m', \
                                              desc='diameters of all turbine rotors'))
        self.add('wakeCentersYT', Array(np.zeros(nTurbines*nTurbines), iotype='in', units='m', \
                                              desc='Y positions of all wakes at each turbine'))
        self.add('wakeDiametersT', Array(np.zeros(nTurbines*nTurbines*3), iotype='in', units='m',\
                                               desc='diameters of all turbines wake zones'))

        # Explicitly size output arrays
        self.add('wakeOverlapTRel', Array(np.zeros(nTurbines*nTurbines*3), iotype='out', \
                                                desc='relative wake zone overlap to rotor area'))
        self.add('rotorArea', Array(np.zeros(nTurbines), iotype='in', units='m*m', desc='Area of each turbine rotor'))


    def execute(self):

        nTurbines = self.turbineYw.size

        wakeDiametersT_mat = np.zeros((nTurbines, nTurbines, 3))
        wakeCentersYT_mat = np.zeros((nTurbines, nTurbines))

        # convert the input vector to the array used for calculations
        for i in range(0, nTurbines):
                wakeDiametersT_mat[i, :, 0] = self.wakeDiametersT[3*nTurbines*i:3*nTurbines*i+nTurbines]
                wakeDiametersT_mat[i, :, 1] = self.wakeDiametersT[3*nTurbines*i+nTurbines:3*nTurbines*i+2*nTurbines]
                wakeDiametersT_mat[i, :, 2] = self.wakeDiametersT[3*nTurbines*i+2*nTurbines:3*nTurbines*i+3*nTurbines]

        # convert the input vector to the array used for calculations
        for i in range(0, nTurbines):
                wakeCentersYT_mat[i, :] = self.wakeCentersYT[nTurbines*i:nTurbines*i+nTurbines]


        # calculate overlap areas at rotors
        # wakeOverlapT(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake
        # of turbine TURB with rotor of turbine TURBI

        rotorArea = np.pi*self.rotorDiameter**2/4.

        wakeOverlapT = calcOverlapAreas(self.turbineXw, self.turbineYw, self.rotorDiameter, wakeDiametersT_mat, wakeCentersYT_mat)

        # make overlap relative to rotor area (maximum value should be 1)
        # wakeOverlapTRel = wakeOverlapT
        wakeOverlapTRel_mat = wakeOverlapT
        for turb in range(0, nTurbines):
            # wakeOverlapTRel[turb] = wakeOverlapTRel[turb]/self.rotorArea[turb]
            wakeOverlapTRel_mat[turb] = wakeOverlapTRel_mat[turb]/rotorArea[turb]


        # convert matrix format to vector format (all are of type ndarray)
        wakeOverlapTRel_vec = np.zeros(3*nTurbines**2)
        # print 'shape of woTRel is %s' %wakeOverlapTRel_vec.shape
        for i in range(0, nTurbines):
            wakeOverlapTRel_vec[(3*nTurbines*i):(3*nTurbines*i+nTurbines)] = wakeOverlapTRel_mat[i, :, 0]
            wakeOverlapTRel_vec[3*nTurbines*i+nTurbines:3*nTurbines*i+2*nTurbines] = wakeOverlapTRel_mat[i, :, 1]
            wakeOverlapTRel_vec[3*nTurbines*i+2*nTurbines:3*nTurbines*i+3*nTurbines] = wakeOverlapTRel_mat[i, :, 2]

        # self.wakeOverlapTRel = wakeOverlapTRel
        self.wakeOverlapTRel = wakeOverlapTRel_vec


class floris_power(Component):
    """ Calculates the turbine power and effective wind speed for each turbine """

    # original variables in Pieter's OpenMDAO stand-alone version of FLORIS
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    # Flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    air_density = Float(iotype='in', units='kg/(m*m*m)', desc='air density in free stream')

    # output variables added so I don't have to use WISDEM while developing gradients
    power = Float(iotype='out', units='kW', desc='total power output of the wind farm')


    def __init__(self, nTurbines, nSamples=0):
        super(floris_power, self).__init__()

        # Explicitly size input arrays
        # input variables added so I don't have to use WISDEM while developing gradients
        self.add('rotorDiameter', Array(np.zeros(nTurbines), dtype='float', iotype='in', units='m', \
                                        desc='rotor diameters of all turbines'))
        self.add('axialInduction', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                                         desc='axial induction of all turbines'))
        self.add('Ct', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                             desc='Thrust coefficient for all turbines'))
        self.add('Cp', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                             desc='power coefficient for all turbines'))
        self.add('generator_efficiency', Array(np.zeros(nTurbines), iotype='in', dtype='float', \
                                               desc='generator efficiency of all turbines'))
        self.add('turbineXw', Array(np.zeros(nTurbines), iotype='in', dtype='float', units='m', \
                                    desc='X positions of turbines in the wind direction reference frame'))
        self.add('wakeCentersYT',  Array(np.zeros([nTurbines, nTurbines]), iotype='in', units='m', \
                                         desc='centers of the wakes at each turbine'))
        self.add('wakeDiametersT', Array(np.zeros([nTurbines, nTurbines, 3]), iotype='in', units='m', \
                                         desc='diameters of each of the wake zones for each of the wakes \
                                         at each turbine'))
        self.add('wakeOverlapTRel', Array(np.zeros([nTurbines, nTurbines, 3]), iotype='in', units='m', \
                                          desc='ratios of wake overlap area per zone to rotor area'))
        self.add('yaw', Array(np.zeros([nTurbines]), iotype='in'))
        # input variables added so I don't have to use WISDEM while developing gradients
        self.add('rotorArea', Array(np.zeros(nTurbines), iotype='in', dtype='float', units='m*m', desc='rotor area of all turbines'))
        self.add('wsw_position', Array(np.zeros([3, nSamples]), iotype='in', units='m', desc='positions where measurements are desired in the windframe'))
        self.add('wakeDiameters', Array(np.zeros([nSamples, nTurbines, 3]), iotype='in', units='m', desc='diameter of wake zones at measurement points'))
        self.add('wakeCentersY', Array(np.zeros([nSamples, nTurbines]), iotype='in', units='m', desc='Y positions of wakes at measurement points'))
        self.add('wakeCentersZ', Array(np.zeros([nSamples, nTurbines]), iotype='in', units='m', desc='Z positions of wakes at measurement points'))

        # Explicitly size output arrays
        self.add('velocitiesTurbines', Array(np.zeros(nTurbines), iotype='out', units='m/s'))
        self.add('wt_power', Array(np.zeros(nTurbines), iotype='out', units='kW'))

        self.add('ws_array', Array(np.zeros(nSamples), iotype='out', units='m/s', desc='wind speed at measurement locations'))


    def execute(self):

        turbineXw = self.turbineXw
        nTurbines = turbineXw.size

        # wakeOverlapTRel = self.wakeOverlapTRel
        wakeOverlapTRel = np.zeros((nTurbines, nTurbines, 3))

        # convert the input vector to the array used for calculations
        for i in range(0, nTurbines):
            wakeOverlapTRel[i, :, 0] = self.wakeOverlapTRel[3*nTurbines*i:3*nTurbines*i+nTurbines]
            wakeOverlapTRel[i, :, 1] = self.wakeOverlapTRel[3*nTurbines*i+nTurbines:3*nTurbines*i+2*nTurbines]
            wakeOverlapTRel[i, :, 2] = self.wakeOverlapTRel[3*nTurbines*i+2*nTurbines:3*nTurbines*i+3*nTurbines]

        ke = self.parameters.ke
        keCorrArray = self.parameters.keCorrArray
        keCorrCT = self.parameters.keCorrCT
        baselineCT = self.parameters.baselineCT
        CTcorrected = self.parameters.CTcorrected
        CPcorrected = self.parameters.CPcorrected
        pP = self.parameters.pP
        Ct = self.Ct
        Vinf = self.wind_speed
        turbineXw = self.turbineXw
        axialInduction = self.axialInduction
        rotorDiameter = self.rotorDiameter
        rotorArea = np.pi*self.rotorDiameter**2/4.
        rho = self.air_density
        generator_efficiency = self.generator_efficiency
        yaw = self.yaw*np.pi/180.
        Cp = self.Cp
        MU = self.parameters.MU
        aU = self.parameters.aU
        bU = self.parameters.bU
        useaUbU = self.parameters.useaUbU
        shearCoefficientAlpha = self.parameters.shearCoefficientAlpha
        shearZh = self.parameters.shearZh

        velX = self.wsw_position[0][:]
        velY = self.wsw_position[1][:]
        velZ = self.wsw_position[2][:]
        nSamples = np.size(velX)
        wakeCentersY = self.wakeCentersY
        wakeCentersZ = self.wakeCentersZ
        wakeDiameters = self.wakeDiameters

        axialIndProvided = self.parameters.axialIndProvided

        if CTcorrected == False:
            Ct = Ct * (np.cos(yaw)**2)

        if CPcorrected == False:
            Cp = Cp * np.cos(yaw)**pP

        if axialIndProvided:
            axialInd = axialInduction
        else:
            axialInd = np.array([CTtoAxialInd(ct) for ct in Ct])

        # adjust k_e to C_T, adjusted to yaw
        ke = ke + keCorrCT*(Ct-baselineCT) # FT = Ct*0.5*rho*A*(U*cos(yaw))^2, hence, thrust decreases with cos^2
                                                           #   Should ke increase directly with thrust? ==>No - Turbulence characteristics in wind-turbine wakes, A. Crespo"'*, J. Hern'andez b

        # array effects with full or partial wake overlap:
        # use overlap area of zone 1 + 2 of upstream turbines to correct ke
        # Note: array effects only taken into account in calculating
        # velocity deficits, in order not to over-complicate code
        # (avoid loops in calculating overlaps)

        keArray = np.zeros(nTurbines)
        for turb in range(0, nTurbines):
            s = np.sum(wakeOverlapTRel[turb, :, 0]+wakeOverlapTRel[turb, :, 1])
            keArray[turb] = ke[turb]*(1+s*keCorrArray)

        # calculate velocities in full flow field (optional)
        self.ws_array = np.tile(Vinf, nSamples)

        # apply shear profile
        self.ws_array = self.ws_array*(velZ/shearZh)**shearCoefficientAlpha
        
        for turb in range(0, nTurbines):
            if useaUbU:
                mU = MU/np.cos(aU*np.pi/180+bU*yaw[turb]) # CHANGE: ke now only corrected with CT, which is already corrected with yaw
            else:
                mU = MU

            for loc in range(0, nSamples):
                deltax = velX[loc] - turbineXw[turb]
                deltay = velY[loc] - wakeCentersY[loc, turb]
                deltaz = velZ[loc] - wakeCentersZ[loc, turb]
                radiusLoc = np.sqrt(deltay**2+deltaz**2)
                axialIndAndNearRotor = 2*axialInd[turb]

                if deltax > 0 and radiusLoc < wakeDiameters[loc, turb, 0]/2.0:    # check if in zone 1
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*keArray[turb]*(mU[0])*np.maximum(0, deltax))), 2)
                elif deltax > 0 and radiusLoc < wakeDiameters[loc, turb, 1]/2.0:    # check if in zone 2
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*keArray[turb]*(mU[1])*np.maximum(0, deltax))), 2)
                elif deltax > 0 and radiusLoc < wakeDiameters[loc, turb, 2]/2.0:    # check if in zone 3
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*keArray[turb]*(mU[2])*np.maximum(0, deltax))), 2)
                elif deltax <= 0 and radiusLoc < rotorDiameter[turb]/2.0:     # check if axial induction zone in front of rotor
                    reductionFactor = axialIndAndNearRotor*(0.5+np.arctan(2.0*np.minimum(0, deltax)/(rotorDiameter[turb]))/np.pi)
                else:
                    reductionFactor = 0
                self.ws_array[loc] *= (1-reductionFactor)
        #print 'ws_array in floris_power is: ', self.ws_array
        # find effective wind speeds at downstream turbines, then predict power downstream turbine
        self.velocitiesTurbines = np.tile(Vinf, nTurbines)

        for turbI in range(0, nTurbines):

            # find overlap-area weighted effect of each wake zone
            wakeEffCoeff = 0
            for turb in range(0, nTurbines):

                wakeEffCoeffPerZone = 0
                deltax = turbineXw[turbI] - turbineXw[turb]

                if deltax > 0:
                    if useaUbU:
                        mU = MU/np.cos(aU*np.pi/180+bU*yaw[turb])
                    else:
                        mU = MU
                    for zone in range(0, 3):
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + np.power((rotorDiameter[turb])/(rotorDiameter[turb]+2*keArray[turb]*mU[zone]*deltax), 2.0) * wakeOverlapTRel[turbI, turb, zone]

                    wakeEffCoeff = wakeEffCoeff + np.power(axialInd[turb]*wakeEffCoeffPerZone, 2.0)

            wakeEffCoeff = (1 - 2 * np.sqrt(wakeEffCoeff))

            # multiply the inflow speed with the wake coefficients to find effective wind speed at turbine
            self.velocitiesTurbines[turbI] *= wakeEffCoeff

        if self.verbose:
            print "wind speed at turbines %s [m/s]" % self.velocitiesTurbines
            print "rotor area %s" % rotorArea
            print "rho %s" % rho
            print "generator_efficiency %s" % generator_efficiency

        # find turbine powers
        self.wt_power = np.power(self.velocitiesTurbines, 3.0) * (0.5*rho*rotorArea*Cp) * generator_efficiency
        self.wt_power /= 1000  # in kW

        if self.verbose:
            print "powers turbines %s [kW]" % self.wt_power

        self.power = np.sum(self.wt_power)
        

def CTtoAxialInd(CT):
    if CT > 0.96: # Glauert condition
        axial_induction = 0.143+np.sqrt(0.0203-0.6427*(0.889-CT))
    else:
        axial_induction = 0.5*(1-np.sqrt(1-CT))
    return axial_induction


def calcOverlapAreas(turbineX,turbineY,rotorDiameter,wakeDiameters,wakeCenters):
    """calculate overlap of rotors and wake zones (wake zone location defined by wake center and wake diameter)
    turbineX,turbineY is x,y-location of center of rotor

    wakeOverlap(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake of turbine TURB with rotor of downstream turbine
    TURBI"""

    nTurbines = turbineY.size

    wakeOverlap = np.zeros((nTurbines,nTurbines,3))

    for turb in range(0,nTurbines):
        for turbI in range(0,nTurbines):
            if turbineX[turbI] > turbineX[turb]:
                OVdYd = wakeCenters[turbI,turb]-turbineY[turbI]
                OVr = rotorDiameter[turbI]/2
                for zone in range(0,3):
                    OVR = wakeDiameters[turbI,turb,zone]/2
                    OVdYd = abs(OVdYd)
                    if OVdYd != 0:
                        OVL = (-np.power(OVr,2.0)+np.power(OVR,2.0)+np.power(OVdYd,2.0))/(2.0*OVdYd)
                    else:
                        OVL = 0

                    OVz = np.power(OVR,2.0)-np.power(OVL,2.0)

                    if OVz > 0:
                        OVz = np.sqrt(OVz)
                    else:
                        OVz = 0

                    if OVdYd < (OVr+OVR):
                        if OVL < OVR and (OVdYd-OVL) < OVr:
                            wakeOverlap[turbI,turb,zone] = np.power(OVR,2.0)*np.arccos(OVL/OVR) + np.power(OVr,2.0)*np.arccos((OVdYd-OVL)/OVr) - OVdYd*OVz
                        elif OVR > OVr:
                            wakeOverlap[turbI,turb,zone] = np.pi*np.power(OVr,2.0)
                        else:
                            wakeOverlap[turbI,turb,zone] = np.pi*np.power(OVR,2.0)
                    else:
                        wakeOverlap[turbI,turb,zone] = 0

    for turb in range(0,nTurbines):
        for turbI in range(0,nTurbines):
            wakeOverlap[turbI,turb,2] = wakeOverlap[turbI,turb,2]-wakeOverlap[turbI,turb,1]
            wakeOverlap[turbI,turb,1] = wakeOverlap[turbI,turb,1]-wakeOverlap[turbI,turb,0]

    return wakeOverlap
