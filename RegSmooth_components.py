from openmdao.main.api import Component, VariableTree
from openmdao.lib.datatypes.api import Array, Bool, Float, VarTree
import numpy as np
from Parameters import FLORISParameters


class floris_windframe(Component):
    """ Calculates the locations of each turbine in the wind direction reference frame """

    # original variables
    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')
    # position = Array(iotype='in', units='m', desc='position of turbines in original ref. frame')
    turbineX = Array(iotype='in', desc='x positions of turbines in original ref. frame')
    turbineY = Array(iotype='in', desc='y positions of turbines in original ref. frame')

    # variables for verbosity
    Ct = Array(iotype='in', dtype='float')
    Cp = Array(iotype='in', dtype='float', desc='power coefficient for all turbines')
    axialInduction = Array(iotype='in', dtype='float', desc='axial induction of all turbines')
    # yaw = Array(iotype='in', desc='yaw of each turbine')

    # variables for testing wind speed at various locations
    # ws_position = Array(iotype='in', units='m', desc='position of desired measurements in original ref. frame')
    # wsw_position = Array(iotype='out', units='m', deriv_ignore=True, desc='position of desired measurements in wind ref. frame')

    # flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    wind_direction = Float(iotype='in', units='deg', desc='overall wind direction for wind farm')

    # for testing purposes only
    turbineXw = Array(iotype='out', units='m', desc='x coordinates of turbines in wind dir. ref. frame')
    turbineYw = Array(iotype='out', units='m', desc='y coordinates of turbines in wind dir. ref. frame')

    def execute(self):

        print 'entering windframe'
        # print 'windframe turbineX = ', self.turbineX
        # print 'windframe turbineY = ', self.turbineY


        Vinf = self.wind_speed
        windDirection = self.wind_direction*np.pi/180.0

        #variables to satisfy verbosity
        axialInd = self.axialInduction
        Cp = self.Cp
        Ct = self.Ct
        # CTcorrected = self.parameters.CTcorrected
        # yaw = self.yaw*np.pi/180

        # get rotor coefficients, and apply corrections if necesary
        # # Cp = np.hstack(self.wt_layout.wt_array(attr='CP'))
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
        # position = self.position
        # turbineX = position[:, 0]
        # turbineY = position[:, 1]
        turbineX = self.turbineX
        turbineY = self.turbineY
        # print turbineX, turbineY

        # if self.ws_position.any():
        #     velX = self.ws_position[:, 0]
        #     velY = self.ws_position[:, 1]
        # else:
        #     velX = np.zeros([0, 0])
        #     velY = np.zeros([0, 0])

        # convert to downwind-crosswind coordinates
        rotationMatrix = np.array([(np.cos(-windDirection), -np.sin(-windDirection)),
                                   (np.sin(-windDirection), np.cos(-windDirection))])
        # print 'rotation matrix = ', rotationMatrix
        turbineLocations = np.dot(rotationMatrix, np.array([turbineX, turbineY]))
        # print turbineLocations
        self.turbineXw = np.zeros(turbineX.size)
        self.turbineYw = np.zeros(turbineX.size)
        # print self.turbineXw
        self.turbineXw = turbineLocations[0]
        self.turbineYw = turbineLocations[1]
        #print 'windframe.turbineX = %s' %self.turbineX
        # if velX.size>0:
        #     locations = np.dot(rotationMatrix, np.array([velX, velY]))
        #     velX = locations[0]
        #     velY = locations[1]

        # self.wsw_position = np.array([velX, velY])
        #print 'wsw_position in windframe is:', self.wsw_position
        #print 'ws_position in windframe is:', self.ws_position


        # print 'windframe turbineXw = ', self.turbineXw
        # print 'windframe turbineYw = ', self.turbineYw


class floris_wcent_wdiam(Component):
    """ Calculates the center and diameter of each turbine wake at each other turbine """

    parameters = VarTree(FLORISParameters(), iotype='in')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')
    turbineXw = Array(iotype='in', desc='x coordinates of turbines in wind dir. ref. frame')
    turbineYw = Array(iotype='in', desc='y coordinates of turbines in wind dir. ref. frame')
    yaw = Array(iotype='in', desc='yaw of each turbine')
    rotorDiameter = Array(dtype='float', iotype='in', desc='rotor diameter of each turbine')
    Ct = Array(iotype='in', dtype='float', desc='thrust coefficient of each turbine')

    # wsw_position = Array(iotype='in', units='m', desc='positions where measurements are desired in the windframe')

    wakeCentersYT = Array(iotype='out', dtype='float', desc='wake center y position at each turbine')
    wakeDiametersT = Array(iotype='out', dtype='float', desc='wake diameter of each zone of each wake at each turbine')
    # wakeDiameters = Array(iotype='out', dtype='float', desc='wake diameter of each zone of each wake at each turbine')
    # wakeCentersY = Array(iotype='out', units='m', desc='Y positions of wakes at measurement points')

    # p_near0 = Float(iotyp='out', desc='upwind location of diameter spline in rotor diameters')

    def execute(self):

        print 'entering wcent_wdiam'
        # print 'wake c and d turbineXw = ', self.turbineXw
        # print 'wake c and d turbineYw = ', self.turbineYw



         # rename inputs and outputs
        # pP = self.parameters.pP
        kd = self.parameters.kd
        ke = self.parameters.ke
        initialWakeDisplacement = self.parameters.initialWakeDisplacement
        initialWakeAngle = self.parameters.initialWakeAngle
        rotorDiameter = self.rotorDiameter
        Ct = self.Ct
        CTcorrected = self.parameters.CTcorrected
        keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        me = self.parameters.me


        turbineXw = self.turbineXw
        turbineYw = self.turbineYw
        yaw = self.yaw*np.pi/180.0
        # print yaw
        nTurbines = turbineXw.size


        # velX = self.wsw_position[0][:]
        # nLocations = np.size(velX)


        # set spline start and end locations (0 is start location, or upwind, 1 is end location, or downwind)
        p_center0 = 0.25
        p_center1 = 0.25
        p_unity = 0.25
        p_near0 = 1
        p_near1 = np.copy(p_unity)
        p_far0 = np.copy(p_unity)
        p_mix0 = np.copy(p_unity)
        p_mix1 = 0.25

        # if CTcorrected == False:
        #     Ct = Ct * (np.cos(yaw)**2)

        # calculate y-location of wake centers
        # wakeCentersY = np.zeros((nLocations, nTurbines))
        wakeCentersYT_mat = np.zeros((nTurbines, nTurbines))
        # print wakeCentersYT_mat
        for turb in range(0, nTurbines):
            wakeAngleInit = 0.5 * np.sin(yaw[turb]) * Ct[turb] + initialWakeAngle*np.pi/180.0 # change: 4*a*(1-a) --> C_T, initial angle
            # for loc in range(0, nLocations):  # at velX-locations
            #     deltax = np.maximum(velX[loc]-turbineXw[turb], 0)
            #     factor = (2.0*kd*deltax/rotorDiameter[turb])+1.0
            #     wakeCentersY[loc, turb] = turbineYw[turb]-initialWakeDisplacement # initial displacement for no yaw (positive to the right looking downstream)
            #     displacement = (wakeAngleInit*(15.0*(factor**4.0)+(wakeAngleInit**2.0))/((30.0*kd*(factor**5.0))/rotorDiameter[turb]))-(wakeAngleInit*rotorDiameter[turb]*(15.0+(wakeAngleInit**2.0))/(30.0*kd)) # yaw-induced deflection
            #     wakeCentersY[loc, turb] = wakeCentersY[loc, turb] + displacement
            #     # print "displacement %s" % displacement

            for turbI in range(0, nTurbines):  # at turbineX-locations
                # deltax = np.maximum(turbineXw[turbI]-turbineXw[turb], 0.0) #original
                deltax = turbineXw[turbI]-turbineXw[turb]
                factor = (2.0*kd*deltax/rotorDiameter[turb])+1.0
                if turbineXw[turb]+p_center1*rotorDiameter[turb] < turbineXw[turbI]:
                    wakeCentersYT_mat[turbI, turb] = turbineYw[turb]
                    wakeCentersYT_mat[turbI, turb] = wakeCentersYT_mat[turbI, turb]-initialWakeDisplacement # initial displacement for no yaw (positive to the right looking downstream)
                    displacement = (wakeAngleInit*(15.0*(factor**4.0)+(wakeAngleInit**2.0))/((30.0*kd*(factor**5.0))/rotorDiameter[turb]))-(wakeAngleInit*rotorDiameter[turb]*(15.0+(wakeAngleInit**2.0))/(30.0*kd)) # yaw-induced wake center displacement
                    wakeCentersYT_mat[turbI, turb] = wakeCentersYT_mat[turbI, turb] + displacement
                elif turbineXw[turb]+p_center1*rotorDiameter[turb] >= turbineXw[turbI] \
                        >= turbineXw[turb]-p_center0*rotorDiameter[turb]:
                    # set up spline values
                    x = turbineXw[turbI]                                    # point of interest
                    # print 'x = ', x
                    x0 = turbineXw[turb]-p_center0*rotorDiameter[turb]      # start (upwind) point of spline
                    x1 = turbineXw[turb]+p_center1*rotorDiameter[turb]      # end (downwind) point of spline
                    # print 'x0, x1 = ', x0, x1
                    y0 = turbineYw[turb]-initialWakeDisplacement           # y value at start point of spline
                    # print y0, x0
                    dy0 = 0.0                                              # slope of spline at start point
                    # print 'y0, dy0 = ', y0, dy0
                    dx_1 = x1-turbineXw[turb]
                    factor_1 = (2.0*kd*dx_1/rotorDiameter[turb])+1.0
                    # print 'dx_1, factor_1 = ', dx_1, factor_1
                    y1 = turbineYw[turb]
                    y1 += -initialWakeDisplacement # initial displacement for no yaw (positive to the right looking downstream)
                    # print 'y1 = ', y1
                    displacement = (wakeAngleInit*(15.0*(factor_1**4.0)+(wakeAngleInit**2.0))/((30.0*kd*(factor_1**5.0))/rotorDiameter[turb]))-(wakeAngleInit*rotorDiameter[turb]*(15.0+(wakeAngleInit**2.0))/(30.0*kd)) # yaw-induced wake center displacement
                    # print 'displacement = ', displacement
                    # print 'wakeAngleInit = ', wakeAngleInit
                    # print 'kd = ', kd
                    y1 += displacement                                      # y value at end point of spline

                    b = 2.0*kd/rotorDiameter[turb]
                    d = b*dx_1+1.0
                    dy1_yaw = -(5/(d**2)+(wakeAngleInit**2)/(3*(d**8)))+4.0/(d**2.0)
                    dy1 = wakeAngleInit*dy1_yaw    # slope of spline at end point
                    # if turbI == 1:
                    #     print dy1
                    # if turbI == 1 and turb == 0:
                    #     print x1, y1, dy1
                    # print 'right = %s' %dy1

                    # call spline function to determine center point of wake at point of interest.
                    wakeCentersYT_mat[turbI, turb], _ = Hermite_Spline(x, x0, x1, y0, dy0, y1, dy1)
                    # if x1+2 > turbineXw[turbI] > x1-2 or x0+2 > turbineXw[turbI] > x0-2:
                    #     wakeCentersYT_mat[turbI, turb] = 1495

                else:
                    wakeCentersYT_mat[turbI, turb] = turbineYw[turb] - initialWakeDisplacement

        # adjust k_e to C_T, adjusted to yaw
        ke += keCorrCT*(Ct-Region2CT) # FT = Ct*0.5*rho*A*(U*cos(yaw))^2, hence, thrust decreases with cos^2
                                                           #   Should ke increase directly with thrust? ==>No - Turbulence characteristics in wind-turbine wakes, A. Crespo"'*, J. Hern'andez b
        # print wakeCentersYT
        # calculate wake zone diameters at velX-locations
        # wakeDiameters = np.zeros((nLocations, nTurbines, 3))
        wakeDiametersT_mat = np.zeros((nTurbines, nTurbines, 3))
        dwakeDiametersT_dx = np.zeros((nTurbines, nTurbines, 3))

        for turb in range(0, nTurbines):
            wakeDiameter0 = rotorDiameter[turb] * np.cos(yaw[turb]) # CHANGE: initial wake diameter at rotor adjusted to yaw
            # for loc in range(0, nLocations):  # at velX-locations
            #     deltax = velX[loc]-turbineXw[turb]
            #     for zone in range(0, 3):
            #         wakeDiameters[loc, turb, zone] = wakeDiameter0 + 2*ke[turb]*me[zone]*np.maximum(deltax, 0)
            for turbI in range(0, nTurbines):  # at turbineX-locations
                deltax = turbineXw[turbI]-turbineXw[turb]

                # for zone in range(0, 3):
                #     wakeDiametersT_mat[turbI, turb, zone] = np.maximum(wakeDiameter0 + 2*ke[turb]*me[zone]*deltax, 0)
                for zone in range(0, 3):
                    if zone == 0:
                        if turbineXw[turb]+p_near1*rotorDiameter[turb] < turbineXw[turbI]:
                            wakeDiametersT_mat[turbI, turb, zone] = wakeDiameter0+2*ke[turb]*me[zone]*deltax
                            dwakeDiametersT_dx[turbI, turb, zone] = 2*ke[turb]*me[zone]
                        elif turbineXw[turb]+p_near1*rotorDiameter[turb] >= turbineXw[turbI] > turbineXw[turb]-p_unity*rotorDiameter[turb]:

                            x = turbineXw[turbI]                              # x position of interest
                            x1 = turbineXw[turb]-p_unity*rotorDiameter[turb]  # point where all zones have equal diameter
                            x2 = turbineXw[turb]+p_near1*rotorDiameter[turb]  # downwind end point of spline

                            # diameter at upwind point of spline
                            y1 = wakeDiameter0-2*ke[turb]*me[1]*p_unity*rotorDiameter[turb]
                            # derivative of diameter at upwind point of spline w.r.t downwind position
                            dy1_dx = 2*ke[turb]*me[1]

                            # diameter at downwind point of spline
                            y2 = wakeDiameter0+2*ke[turb]*me[zone]*(x2-turbineXw[turb])
                            # derivative of diameter at downwind point of spline w.r.t. downwind position
                            dy2_dx = 2*ke[turb]*me[zone]


                            # solve for the wake zone diameter and its derivative w.r.t. the downwind
                            # location at the point of interest
                            wakeDiametersT_mat[turbI, turb, zone], dwakeDiametersT_dx[turbI, turb, zone] = Hermite_Spline\
                                (x, x1, x2, y1, dy1_dx, y2, dy2_dx)
                            # if turbI == 1 and turb == 0:
                            #     print dy2_dx
                            #     print dwakeDiametersT_dx[turbI, turb, zone]

                        elif turbineXw[turb]-p_near0*rotorDiameter[turb] <= turbineXw[turbI] <= turbineXw[turb]:

                            x = turbineXw[turbI]                            # x position of interest
                            x0 = turbineXw[turb]-p_near0*rotorDiameter[turb]  # downwind end point of spline
                            x1 = turbineXw[turb]-p_unity*rotorDiameter[turb]  # point where all zones have equal diameter

                            # diameter at upwind point of spline
                            y0 = 0
                            # derivative of diameter at upwind point of spline
                            dy0_dx = 0

                            # diameter at upwind point of spline
                            y1 = wakeDiameter0-2*ke[turb]*me[1]*p_unity*rotorDiameter[turb]
                            # derivative of diameter at upwind point of spline
                            dy1_dx = 2*ke[turb]*me[1]

                            # solve for the wake zone diameter and its derivative w.r.t. the downwind
                            # location at the point of interest
                            wakeDiametersT_mat[turbI, turb, zone], dwakeDiametersT_dx[turbI, turb, zone] = Hermite_Spline\
                                (x, x0, x1, y0, dy0_dx, y1, dy1_dx)
                            # if turbI == 1 and turb == 0:
                            #     print dy1_dx
                            #     print dwakeDiametersT_dx[turbI, turb, zone]

                    elif zone == 1:
                        if turbineXw[turb]-p_far0*rotorDiameter[turb] < turbineXw[turbI]:
                            wakeDiametersT_mat[turbI, turb, zone] = wakeDiameter0 + 2*ke[turb]*me[zone]*deltax
                            dwakeDiametersT_dx[turbI, turb, zone] = 2*ke[turb]*me[zone]
                        else:
                            wakeDiametersT_mat[turbI, turb, zone] = wakeDiametersT_mat[turbI, turb, 0]

                    elif zone == 2:

                        if turbineXw[turb]+p_mix1*rotorDiameter[turb] < turbineXw[turbI]:
                            wakeDiametersT_mat[turbI, turb, zone] = wakeDiameter0+2*ke[turb]*me[zone]*deltax
                            dwakeDiametersT_dx[turbI, turb, zone] = 2*ke[turb]*me[zone]

                        elif turbineXw[turb]+p_mix1*rotorDiameter[turb] >= turbineXw[turbI] > \
                                        turbineXw[turb]-p_mix0*rotorDiameter[turb]:
                            x = turbineXw[turbI]                             # x position of interest
                            x0 = turbineXw[turb]-p_mix0*rotorDiameter[turb]  # point where all zones have equal diameter
                            x1 = turbineXw[turb]+p_mix1*rotorDiameter[turb]  # point downwind

                            # diameter at upwind point of spline
                            y0 = wakeDiameter0-2*ke[turb]*me[1]*p_mix0*rotorDiameter[turb]
                            # derivative of diameter at upwind point of spline w.r.t downwind position
                            dy0_dx = 2*ke[turb]*me[1]

                            # diameter at downwind point of spline
                            y1 = wakeDiameter0+2*ke[turb]*me[zone]*p_mix1*rotorDiameter[turb]
                            # derivative of diameter at downwind point of spline w.r.t. downwind position
                            dy1_dx = 2*ke[turb]*me[zone]

                            # solve for the wake zone diameter and its derivative w.r.t. the downwind
                            # location at the point of interest
                            wakeDiametersT_mat[turbI, turb, zone], dwakeDiametersT_dx[turbI, turb, zone] = Hermite_Spline\
                                (x, x0, x1, y0, dy0_dx, y1, dy1_dx)

                        else:
                            wakeDiametersT_mat[turbI, turb, zone] = wakeDiametersT_mat[turbI, turb, 0]

        wakeDiametersT_vec = np.zeros(3*nTurbines*nTurbines)
        # print 'shape of woTRel is %s' %wakeOverlapTRel_vec.shape
        for i in range(0, nTurbines):
            wakeDiametersT_vec[(3*nTurbines*i):(3*nTurbines*i+nTurbines)] = wakeDiametersT_mat[i, :, 0]
            wakeDiametersT_vec[3*nTurbines*i+nTurbines:3*nTurbines*i+2*nTurbines] = wakeDiametersT_mat[i, :, 1]
            wakeDiametersT_vec[3*nTurbines*i+2*nTurbines:3*nTurbines*i+3*nTurbines] = wakeDiametersT_mat[i, :, 2]

        wakeCentersYT_vec = np.zeros(nTurbines*nTurbines)
        # print 'shape of woTRel is %s' %wakeOverlapTRel_vec.shape
        for i in range(0, nTurbines):
            wakeCentersYT_vec[(nTurbines*i):(nTurbines*i+nTurbines)] = wakeCentersYT_mat[i, :]


        self.wakeCentersYT = wakeCentersYT_vec
        self.wakeDiametersT = wakeDiametersT_vec
        # self.wakeDiameters = wakeDiameters
        # self.wakeCentersY = wakeCentersY

        # print self.wakeCentersYT

        # testing
        # self.p_near0 = p_near0

    #
    # def list_deriv_vars(self):
    #
    #     return ('')
    #
    # def provideJ(self):
    #
    #
    #     return J


class floris_overlap(Component):
    """ Calculates the overlap between each turbine rotor and the existing turbine wakes """
    turbineXw = Array(iotype='in', units='m', desc='X positions of turbines wrt the wind direction')
    turbineYw = Array(iotype='in', units='m', desc='Y positions of turbines wrt the wind direction')
    rotorDiameter = Array(iotype='in', units='m', desc='diameters of all turbine rotors')
    wakeDiametersT = Array(iotype='in', units='m', desc='diameters of all turbines wake zones')
    wakeCentersYT = Array(iotype='in', units='m', desc='Y positions of all wakes at each turbine')
    rotorArea = Array(iotype='in', units='m*m', desc='Area of each turbine rotor')

    wakeOverlapTRel = Array(iotype='out', desc='relative wake zone overlap to rotor area')

    # p_near0 = Float(iotype='in', desc='upwind location of diameter spline in rotor diameters')

    def execute(self):
        print 'entering overlap'
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

        # p_near0 = self.p_near0
        p_near0 = 1
        # calculate overlap areas at rotors
        # wakeOverlapT(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake
        # of turbine TURB with rotor of turbine TURBI
        # wakeOverlapT = calcOverlapAreas(self.turbineXw, self.turbineYw, self.rotorDiameter, self.wakeDiametersT, self.wakeCentersYT, p_near0)

        # make overlap relative to rotor area (maximum value should be 1)
        # wakeOverlapTRel = wakeOverlapT
        # for turb in range(0, nTurbines): # Jared: I think it would make more sense to use turbI for consistency
        #     wakeOverlapTRel[turb] = wakeOverlapTRel[turb]/self.rotorArea[turb]

        # wakeOverlapTRel = calcOverlapAreas(self.turbineXw, self.turbineYw, self.rotorDiameter, self.wakeDiametersT, self.wakeCentersYT, p_near0)
        wakeOverlapTRel_mat = calcOverlapAreas(self.turbineXw, self.turbineYw, self.rotorDiameter, wakeDiametersT_mat, wakeCentersYT_mat, p_near0)

        # self.wakeOverlapTRel = wakeOverlapTRel
        # print self.wakeOverlapTRel
        # print '_'

        # convert matrix format to vector format (all are of type ndarray)
        wakeOverlapTRel_vec = np.zeros(3*nTurbines*nTurbines)
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
    velocitiesTurbines = Array(iotype='out', units='m/s')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    # input variables added so I don't have to use WISDEM while developing gradients
    rotorDiameter = Array(dtype='float', iotype='in', units='m', desc='rotor diameters of all turbine')
    rotorArea = Array(iotype='in', dtype='float', units='m*m', desc='rotor area of all turbines')
    axialInduction = Array(iotype='in', dtype='float', desc='axial induction of all turbines')
    Ct = Array(iotype='in', dtype='float', desc='Thrust coefficient for all turbines')
    Cp = Array(iotype='in', dtype='float', desc='power coefficient for all turbines')
    generator_efficiency = Array(iotype='in', dtype='float', desc='generator efficiency of all turbines')
    yaw = Array(iotype='in', desc='yaw of each turbine')
    turbineXw = Array(iotype='in', dtype='float', units='m', desc='X positions of turbines in the wind direction reference frame')
    wakeCentersYT = Array(iotype='in', units='m', desc='centers of the wakes at each turbine')
    wakeDiametersT = Array(iotype='in', units='m', desc='diameters of each of the wake zones for each of the wakes at each turbine')
    wakeOverlapTRel = Array(iotype='in', units='m', desc='ratios of wake overlap area per zone to rotor area')
    wsw_position = Array(iotype='in', units='m', desc='positions where measurements are desired in the windframe')
    # wakeDiameters = Array(iotype='in', units='m', desc='diameter of wake zones at measurement points')
    # wakeCentersY = Array(iotype='in', units='m', desc='Y positions of wakes at measurement points')

    # Flow property variables
    wind_speed = Float(iotype='in', units='m/s', desc='free stream wind velocity')
    air_density = Float(iotype='in', units='kg/(m*m*m)', desc='air density in free stream')

    # output variables added so I don't have to use WISDEM while developing gradients
    wt_power = Array(iotype='out', units='kW')
    power = Float(iotype='out', units='kW', desc='total power output of the wind farm')
    # ws_array = Array(iotype='out', units='m/s', desc='wind speed at measurement locations')

    def execute(self):
        print 'entering power'
        turbineXw = self.turbineXw
        #print 'turbineXw = %s' %turbineXw
        nTurbines = turbineXw.size
        #print 'number of turbines is: ', nTurbines

        wakeOverlapTRel = np.zeros((nTurbines, nTurbines, 3))

        # print 'wakeOverlapTRel = %s' %wakeOverlapTRel

        # convert the input vector to the array used for calculations
        for i in range(0, nTurbines):
                wakeOverlapTRel[i, :, 0] = self.wakeOverlapTRel[3*nTurbines*i:3*nTurbines*i+nTurbines]
                wakeOverlapTRel[i, :, 1] = self.wakeOverlapTRel[3*nTurbines*i+nTurbines:3*nTurbines*i+2*nTurbines]
                wakeOverlapTRel[i, :, 2] = self.wakeOverlapTRel[3*nTurbines*i+2*nTurbines:3*nTurbines*i+3*nTurbines]

        # print wakeOverlapTRel[0,0,0]

        # wakeOverlapTRel = self.wakeOverlapTRel
        ke = self.parameters.ke
        #print 'ke is: ', ke
        keCorrArray = self.parameters.keCorrArray
        keCorrCT = self.parameters.keCorrCT
        Region2CT = self.parameters.Region2CT
        CTcorrected = self.parameters.CTcorrected
        Ct = self.Ct
        Vinf = self.wind_speed
        turbineXw = self.turbineXw
        axialInduction = self.axialInduction
        rotorDiameter = self.rotorDiameter
        rotorArea = self.rotorArea
        rho = self.air_density
        generator_efficiency = self.generator_efficiency
        yaw = self.yaw
        Cp = self.Cp
        MU = self.parameters.MU

        # velX = self.wsw_position[0][:]
        # velY = self.wsw_position[1][:]
        # print 'wsw_position in power is: ', self.wsw_position
        # nLocations = np.size(velX)
        #print 'nlocations is', nLocations
        #print 'nturbines is', nTurbines
        # wakeCentersY = self.wakeCentersY
        # wakeDiameters = self.wakeDiameters

        axialIndProvided = self.parameters.axialIndProvided



        # how far upwind (in rotor diameters) to calculate power (must correspond to the value in wake overlap calculations)
        p_near0 = 1
        #print 'Vinf = %s' %Vinf

        # if CTcorrected == False:
        #     Ct = Ct * (np.cos(yaw)**2)
        #
        if axialIndProvided:
            axialInd = axialInduction
        else:
            axialInd = np.array([CTtoAxialInd(ct) for ct in Ct])



        # adjust k_e to C_T, adjusted to yaw
        ke = ke + keCorrCT*(Ct-Region2CT) # FT = Ct*0.5*rho*A*(U*cos(yaw))^2, hence, thrust decreases with cos^2
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
        # self.ws_array = np.tile(Vinf, nLocations)
        #print 'nLocations in power is: ', nLocations
        # for turb in range(0, nTurbines):
            #mU = MU/np.cos(aU*np.pi/180+bU*yaw[turb]) // CHANGE: ke now only corrected with CT, which is already corrected with yaw
            # for loc in range(0, nLocations):
            #     deltax = velX[loc] - turbineXw[turb]
            #     radiusLoc = abs(velY[loc]-wakeCentersY[loc, turb])
            #     axialIndAndNearRotor = 2*axialInd[turb]
            #
            #     if deltax > 0 and radiusLoc < wakeDiameters[loc, turb, 0]/2.0:    # check if in zone 1
            #         reductionFactor = axialIndAndNearRotor*\
            #                           np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*keArray[turb]*(MU[0])*np.maximum(0, deltax))), 2)
            #     elif deltax > 0 and radiusLoc < wakeDiameters[loc, turb, 1]/2.0:    # check if in zone 2
            #         reductionFactor = axialIndAndNearRotor*\
            #                           np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*keArray[turb]*(MU[1])*np.maximum(0, deltax))), 2)
            #     elif deltax > 0 and radiusLoc < wakeDiameters[loc, turb, 2]/2.0:    # check if in zone 3
            #         reductionFactor = axialIndAndNearRotor*\
            #                           np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*keArray[turb]*(MU[2])*np.maximum(0, deltax))), 2)
            #     elif deltax <= 0 and radiusLoc < rotorDiameter[turb]/2.0:     # check if axial induction zone in front of rotor
            #         reductionFactor = axialIndAndNearRotor*(0.5+np.arctan(2.0*np.minimum(0, deltax)/(rotorDiameter[turb]))/np.pi)
            #     else:
            #         reductionFactor = 0
            #     self.ws_array[loc] *= (1-reductionFactor)
        # print 'ws_array in floris_power is: ', self.ws_array
        # find effective wind speeds at downstream turbines, then predict power downstream turbine
        self.velocitiesTurbines = np.tile(Vinf, nTurbines)
        #print 'vel = %s' %self.velocitiesTurbines
        for turbI in range(0, nTurbines):

            # find overlap-area weighted effect of each wake zone
            wakeEffCoeff = 0
            for turb in range(0, nTurbines):

                wakeEffCoeffPerZone = 0
                deltax = turbineXw[turbI] - turbineXw[turb]

                if deltax > -1*p_near0*rotorDiameter[turb] and turbI != turb:
                # if deltax > 0:
                    #mU = MU / np.cos(aU*np.pi/180 + bU*yaw[turb]) // CHANGE: ke now only corrected with CT, which is already corrected with yaw
                    for zone in range(0, 3):
                        #wakeEffCoeffPerZone = wakeEffCoeffPerZone + np.power((rotorDiameter[turb])/(rotorDiameter[turb]+2*ke[turb]*mU[zone]*deltax),2.0) * wakeOverlapTRel[turbI,turb,zone]
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + np.power((rotorDiameter[turb])/(rotorDiameter[turb]+2*ke[turb]*MU[zone]*deltax), 2.0) * wakeOverlapTRel[turbI, turb, zone]

                    wakeEffCoeff = wakeEffCoeff + np.power(axialInd[turb]*wakeEffCoeffPerZone, 2.0)

            wakeEffCoeff = (1 - 2 * np.sqrt(wakeEffCoeff))
            #print wakeEffCoeff

            # multiply the inflow speed with the wake coefficients to find effective wind speed at turbine
            self.velocitiesTurbines[turbI] *= wakeEffCoeff

        if self.verbose:
            print "wind speed at turbines %s [m/s]" % self.velocitiesTurbines
            print "rotor area %s" % rotorArea
            print "rho %s" % rho
            print "generator_efficiency %s" % generator_efficiency

        # find turbine powers
        #print self.velocitiesTurbines
        self.wt_power = np.power(self.velocitiesTurbines, 3.0) * (0.5*rho*rotorArea*Cp) * generator_efficiency

        # # set outputs on turbine level
        # for turbI in range(0, nTurbines):
        #     turbineName = self.wt_layout.wt_names[turbI]
        #     getattr(self.wt_layout, turbineName).power = self.wt_power[turbI] # in W
        #     getattr(self.wt_layout, turbineName).wind_speed_eff = self.velocitiesTurbines[turbI]

        self.wt_power /= 1000  # in kW

        if self.verbose:
            print "powers turbines %s [kW]" % self.wt_power

        self.power = np.sum(self.wt_power)
        print self.power


def Hermite_Spline(x, x0, x1, y0, dy0, y1, dy1):
    """
    This function produces the y and dy values for a hermite cubic spline
    interpolating between two end points with known slopes

    :param x: x position of output y
    :param x0: x position of upwind endpoint of spline
    :param x1: x position of downwind endpoint of spline
    :param y0: y position of upwind endpoint of spline
    :param dy0: slope at upwind endpoint of spline
    :param y1: y position of downwind endpoint of spline
    :param dy1: slope at downwind endpoint of spline

    :return: [y: y value of spline at location x, dy: slope of spline at location x]

    """

    # initialize coefficients for parametric cubic spline
    # c3 = (2*(y1))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) - (2*(y0))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) + (dy0)/(x0**2 - 2*x0*x1 + x1**2) + (dy1)/(x0**2 - 2*x0*x1 + x1**2)
    # c2 = (3*(y0)*(x0 + x1))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) - ((dy1)*(2*x0 + x1))/(x0**2 - 2*x0*x1 + x1**2) - ((dy0)*(x0 + 2*x1))/(x0**2 - 2*x0*x1 + x1**2) - (3*(y1)*(x0 + x1))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3)
    # c1 = ((dy0)*(x1**2 + 2*x0*x1))/(x0**2 - 2*x0*x1 + x1**2) + ((dy1)*(x0**2 + 2*x1*x0))/(x0**2 - 2*x0*x1 + x1**2) - (6*x0*x1*(y0))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) + (6*x0*x1*(y1))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3)
    # c0 = ((y0)*(- x1**3 + 3*x0*x1**2))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) - ((y1)*(- x0**3 + 3*x1*x0**2))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) - (x0*x1**2*(dy0))/(x0**2 - 2*x0*x1 + x1**2) - (x0**2*x1*(dy1))/(x0**2 - 2*x0*x1 + x1**2)
    #
    # # Solve for y and dy values at the given point
    # y = c3*x**3 + c2*x**2 + c1*x + c0
    # dy = c3*3*x**3 + c2*2*x**2 + c1*x
    # print dy
    # print 'y = %s' %y

    # initialize coefficients for parametric cubic spline
    c3 = (2*(y1))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) - (2*(y0))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) + (dy0)/(x0**2 - 2*x0*x1 + x1**2) + (dy1)/(x0**2 - 2*x0*x1 + x1**2)
    c2 = (3*(y0)*(x0 + x1))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) - ((dy1)*(2*x0 + x1))/(x0**2 - 2*x0*x1 + x1**2) - ((dy0)*(x0 + 2*x1))/(x0**2 - 2*x0*x1 + x1**2) - (3*(y1)*(x0 + x1))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3)
    c1 = ((dy0)*(x1**2 + 2*x0*x1))/(x0**2 - 2*x0*x1 + x1**2) + ((dy1)*(x0**2 + 2*x1*x0))/(x0**2 - 2*x0*x1 + x1**2) - (6*x0*x1*(y0))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) + (6*x0*x1*(y1))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3)
    c0 = ((y0)*(- x1**3 + 3*x0*x1**2))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) - ((y1)*(- x0**3 + 3*x1*x0**2))/(x0**3 - 3*x0**2*x1 + 3*x0*x1**2 - x1**3) - (x0*x1**2*(dy0))/(x0**2 - 2*x0*x1 + x1**2) - (x0**2*x1*(dy1))/(x0**2 - 2*x0*x1 + x1**2)

    # print 'c3 = ', c3
    # print 'c2 = ', c2
    # print 'c1 = ', c1
    # print 'c0 = ', c0

    # Solve for y and dy values at the given point
    y = c3*x**3 + c2*x**2 + c1*x + c0
    dy_dx = c3*3*x**2 + c2*2*x + c1


    return y, dy_dx


def CTtoAxialInd(CT):
    if CT > 0.96: # Glauert condition
        axial_induction = 0.143+np.sqrt(0.0203-0.6427*(0.889-CT))
    else:
        axial_induction = 0.5*(1-np.sqrt(1-CT))
    return axial_induction


def calcOverlapAreas(turbineX,turbineY,rotorDiameter,wakeDiameters,wakeCenters,p_near0):
    """calculate overlap of rotors and wake zones (wake zone location defined by wake center and wake diameter)
    turbineX,turbineY is x,y-location of center of rotor

    wakeOverlap(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake of turbine TURB with rotor of downstream turbine
    TURBI"""

    nTurbines = turbineY.size

    wakeOverlap = np.zeros((nTurbines, nTurbines, 3))

    for turb in range(0, nTurbines):
        for turbI in range(0, nTurbines):
            if turbineX[turbI] > turbineX[turb] - p_near0*rotorDiameter[turb]:
                OVdYd = wakeCenters[turbI, turb]-turbineY[turbI]    # distance between wake center and rotor center
                OVr = rotorDiameter[turbI]/2    # rotor diameter
                for zone in range(0, 3):
                    OVR = wakeDiameters[turbI, turb, zone]/2    # wake diameter
                    OVdYd = abs(OVdYd)
                    if OVdYd != 0:
                        # calculate the distance from the wake center to the vertical line between
                        # the two circle intersection points
                        OVL = (-np.power(OVr, 2.0)+np.power(OVR, 2.0)+np.power(OVdYd, 2.0))/(2.0*OVdYd)
                    else:
                        OVL = 0

                    OVz = np.power(OVR, 2.0)-np.power(OVL, 2.0)

                    # Finish calculating the distance from the intersection line to the outer edge of the wake zone
                    if OVz > 0:
                        OVz = np.sqrt(OVz)
                    else:
                        OVz = 0

                    if OVdYd < (OVr+OVR): # if the rotor overlaps the wake zone
                        # if
                        if OVL < OVR and (OVdYd-OVL) < OVr:
                            wakeOverlap[turbI, turb, zone] = np.power(OVR, 2.0)*np.arccos(OVL/OVR) + np.power(OVr, 2.0)*np.arccos((OVdYd-OVL)/OVr) - OVdYd*OVz
                        elif OVR > OVr:
                            wakeOverlap[turbI, turb, zone] = np.pi*np.power(OVr, 2.0)
                        else:
                            wakeOverlap[turbI, turb, zone] = np.pi*np.power(OVR, 2.0)
                    else:
                        wakeOverlap[turbI, turb, zone] = 0


    for turb in range(0, nTurbines):
        for turbI in range(0, nTurbines):
            wakeOverlap[turbI, turb, 2] = wakeOverlap[turbI, turb, 2]-wakeOverlap[turbI, turb, 1]
            wakeOverlap[turbI, turb, 1] = wakeOverlap[turbI, turb, 1]-wakeOverlap[turbI, turb, 0]

    wakeOverlapTRel = wakeOverlap

    for turbI in range(0, nTurbines): # Jared: I think it would make more sense to use turbI for consistency
            wakeOverlapTRel[turbI] = wakeOverlapTRel[turbI]/((np.pi*rotorDiameter[turbI]**2)/4)


    return wakeOverlapTRel