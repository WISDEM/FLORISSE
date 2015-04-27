from fusedwind.plant_flow.comp import GenericWindFarm, GenericFlowModel
from openmdao.lib.datatypes.api import Array, Float, VarTree
from openmdao.main.api import VariableTree

import numpy as np

class FLORISParameters(VariableTree):
    """Container of FLORIS wake parameters"""

    pP = Float(1.88, iotype='in')
    ke = Float(0.065, iotype='in')
    keCorrDA = Float(0, iotype='in')
    kd = Float(0.15, iotype='in')
    me = Array(np.array([-0.5,0.22,1.0]), iotype='in')
    MU = Array(np.array([0.5,1.0,5.5]), iotype='in')
    ad = Float(-4.5, iotype='in')
    bd = Float(-0.01, iotype='in')
    aU = Float(5.0, iotype='in', units='deg')
    bU = Float(1.66, iotype='in')

class FLORIS(GenericWindFarm, GenericFlowModel):

    """Evaluates the FLORIS model and gives the FLORIS-predicted powers of the turbines at locations turbineX, turbineY,
    and, optionally, the FLORIS-predicted velocities at locations (velX,velY)"""

    parameters = VarTree(FLORISParameters(), iotype='in')
    velocitiesTurbines = Array(iotype='out', units='m/s')

    def execute(self):

        # rename inputs and outputs
        Vinf = self.wind_speed                               # from standard FUSED-WIND GenericWindFarm component
        windDirection = self.wind_direction*np.pi/180        # from standard FUSED-WIND GenericWindFarm component
        rho = self.air_density

        rotorDiameter = np.hstack(self.wt_layout.wt_array(attr='rotor_diameter')) # from standard FUSED-WIND GenericWindFarm component
        rotorArea = np.hstack(self.wt_layout.wt_array(attr='rotor_area'))         # from standard FUSED-WIND GenericWindFarm component
        yaw = np.hstack(self.wt_layout.wt_array(attr='yaw'))from fusedwind.plant_flow.comp import GenericWindFarm, GenericFlowModel
from openmdao.lib.datatypes.api import Array, Float, VarTree, Bool
from openmdao.main.api import VariableTree

import numpy as np

class FLORISParameters(VariableTree):
    """Container of FLORIS wake parameters"""

    pP = Float(1.88, iotype='in')
    ke = Float(0.065, iotype='in')
    keCorrDA = Float(0, iotype='in')
    kd = Float(0.15, iotype='in')
    me = Array(np.array([-0.5,0.22,1.0]), iotype='in')
    MU = Array(np.array([0.5,1.0,5.5]), iotype='in')
    ad = Float(-4.5, iotype='in')
    bd = Float(-0.01, iotype='in')
    aU = Float(5.0, iotype='in', units='deg')
    bU = Float(1.66, iotype='in')

class FLORIS(GenericWindFarm, GenericFlowModel):

    """Evaluates the FLORIS model and gives the FLORIS-predicted powers of the turbines at locations turbineX, turbineY,
    and, optionally, the FLORIS-predicted velocities at locations (velX,velY)"""

    parameters = VarTree(FLORISParameters(), iotype='in')
    velocitiesTurbines = Array(iotype='out', units='m/s')
    verbose = Bool(False, iotype='in', desc='verbosity of FLORIS, False is no output')

    def execute(self):

        # rename inputs and outputs
        Vinf = self.wind_speed                               # from standard FUSED-WIND GenericWindFarm component
        windDirection = self.wind_direction*np.pi/180.0        # from standard FUSED-WIND GenericWindFarm component
        rho = self.air_density

        rotorDiameter = np.hstack(self.wt_layout.wt_array(attr='rotor_diameter')) # from standard FUSED-WIND GenericWindFarm component
        rotorArea = np.hstack(self.wt_layout.wt_array(attr='rotor_area'))         # from standard FUSED-WIND GenericWindFarm component
        yaw = np.hstack(self.wt_layout.wt_array(attr='yaw'))*np.pi/180.0
        Cp = np.hstack(self.wt_layout.wt_array(attr='CP'))
        axialInd = np.hstack(self.wt_layout.wt_array(attr='axial_induction'))

        if self.verbose:
            np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print "wind direction %s deg" % [windDirection*180.0/np.pi]
            print "free-stream wind speed %s" % Vinf
            print "axial induction turbines %s" % axialInd
            print "Cp turbines %s" % Cp

        positions = self.wt_layout.wt_array(attr='position')
        turbineX = positions[:, 0]
        turbineY = positions[:, 1]

        if self.ws_positions.any():
            velX = self.ws_positions[:, 0]
            velY = self.ws_positions[:, 1]
        else:
            velX = np.zeros([0,0])
            velY = np.zeros([0,0])

        pP = self.parameters.pP
        ke = self.parameters.ke
        keCorrDA = self.parameters.keCorrDA
        kd = self.parameters.kd
        me = self.parameters.me
        MU = self.parameters.MU
        ad = self.parameters.ad
        bd = self.parameters.bd
        aU = self.parameters.aU
        bU = self.parameters.bU

        # find size of arrays
        nTurbines = turbineX.size
        nLocations = velX.size

        # convert to downwind-crosswind coordinates
        rotationMatrix = np.array([(np.cos(-windDirection), -np.sin(-windDirection)),
                                   (np.sin(-windDirection), np.cos(-windDirection))])
        turbineLocations = np.dot(rotationMatrix,np.array([turbineX,turbineY]))
        turbineX = turbineLocations[0]
        turbineY = turbineLocations[1]

        if velX.size>0:
            locations = np.dot(rotationMatrix,np.array([velX,velY]))
            velX = locations[0]
            velY = locations[1]

        # calculate y-location of wake centers
        wakeCentersY = np.zeros((nLocations,nTurbines))
        wakeCentersYT = np.zeros((nTurbines,nTurbines))
        for turb in range(0,nTurbines):
            wakeAngleInit = 0.5*np.power(np.cos(yaw[turb]),2)*np.sin(yaw[turb])*4*axialInd[turb]*(1-axialInd[turb])
            for loc in range(0,nLocations):  # at velX-locations
                deltax = np.maximum(velX[loc]-turbineX[turb],0)
                factor = (2*kd*deltax/rotorDiameter[turb])+1
                wakeCentersY[loc,turb] = turbineY[turb]
                wakeCentersY[loc,turb] = wakeCentersY[loc,turb] + ad+bd*deltax  # rotation-induced deflection
                wakeCentersY[loc,turb] = wakeCentersY[loc,turb] + \
                                         (wakeAngleInit*(15*np.power(factor,4)+np.power(wakeAngleInit,2))/((30*kd*np.power(factor,5))/rotorDiameter[turb]))-(wakeAngleInit*rotorDiameter[turb]*(15+np.power(wakeAngleInit,4))/(30*kd)) # yaw-induced deflection
            for turbI in range(0,nTurbines):  # at turbineX-locations
                deltax = np.maximum(turbineX[turbI]-turbineX[turb],0)
                factor = (2*kd*deltax/rotorDiameter[turb])+1
                wakeCentersYT[turbI,turb] = turbineY[turb]
                wakeCentersYT[turbI,turb] = wakeCentersYT[turbI,turb] + ad+bd*deltax  # rotation-induced deflection
                wakeCentersYT[turbI,turb] = wakeCentersYT[turbI,turb] + \
                                            (wakeAngleInit*(15*np.power(factor,4)+np.power(wakeAngleInit,4))/((30*kd*np.power(factor,5))/rotorDiameter[turb]))-(wakeAngleInit*rotorDiameter[turb]*(15+np.power(wakeAngleInit,4))/(30*kd)) # yaw-induced deflection

        # calculate wake zone diameters at velX-locations
        wakeDiameters = np.zeros((nLocations,nTurbines,3))
        wakeDiametersT = np.zeros((nTurbines,nTurbines,3))
        for turb in range(0,nTurbines):
            for loc in range(0,nLocations):  # at velX-locations
                deltax = velX[loc]-turbineX[turb]
                for zone in range(0,3):
                    wakeDiameters[loc,turb,zone] = rotorDiameter[turb]+2*ke*me[zone]*np.maximum(deltax,0)
            for turbI in range(0,nTurbines):  # at turbineX-locations
                deltax = turbineX[turbI]-turbineX[turb]
                for zone in range(0,3):
                    wakeDiametersT[turbI,turb,zone] = np.maximum(rotorDiameter[turb]+2*ke*me[zone]*deltax,0)

        # calculate overlap areas at rotors
        # wakeOverlapT(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake
        # of turbine TURB with rotor of turbine TURBI
        wakeOverlapT = calcOverlapAreas(turbineX,turbineY,rotorDiameter,wakeDiametersT,wakeCentersYT)

        # make overlap relative to rotor area (maximum value should be 1)
        wakeOverlapTRel = wakeOverlapT
        for turb in range(0,nTurbines):
            wakeOverlapTRel[turb] = wakeOverlapTRel[turb]/rotorArea[turb]

        # array effects with full or partial wake overlap:
        # use overlap area of zone 1+2 of upstream turbines to correct ke
        # Note: array effects only taken into account in calculating
        # velocity deficits, in order not to over-complicate code
        # (avoid loops in calculating overlaps)

        keUncorrected = ke
        ke = np.zeros(nTurbines)
        for turb in range(0,nTurbines):
            s = np.sum(wakeOverlapTRel[turb,:,0]+wakeOverlapTRel[turb,:,1])
            ke[turb] = keUncorrected*(1+s*keCorrDA)

        # calculate velocities in full flow field (optional)
        self.ws_array = np.tile(Vinf,nLocations)
        for turb in range(0,nTurbines):
            mU = MU/np.cos(aU*np.pi/180+bU*yaw[turb])
            for loc in range(0,nLocations):
                deltax = velX[loc] - turbineX[turb]
                radiusLoc = abs(velY[loc]-wakeCentersY[loc,turb])
                axialIndAndNearRotor = 2*axialInd[turb]

                if deltax > 0 and radiusLoc < wakeDiameters[loc,turb,0]/2.0:    # check if in zone 1
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*ke[turb]*(mU[0])*np.maximum(0,deltax))),2)
                elif deltax > 0 and radiusLoc < wakeDiameters[loc,turb,1]/2.0:    # check if in zone 2
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*ke[turb]*(mU[1])*np.maximum(0,deltax))),2)
                elif deltax > 0 and radiusLoc < wakeDiameters[loc,turb,2]/2.0:    # check if in zone 3
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*ke[turb]*(mU[2])*np.maximum(0,deltax))),2)
                elif deltax <= 0 and radiusLoc < rotorDiameter[turb]/2.0:     # check if axial induction zone in front of rotor
                    reductionFactor = axialIndAndNearRotor*(0.5+np.arctan(2.0*np.minimum(0,deltax)/(rotorDiameter[turb]))/np.pi)
                else:
                    reductionFactor = 0
                self.ws_array[loc] *= (1-reductionFactor)

        # find effective wind speeds at downstream turbines, then predict power downstream turbine
        self.velocitiesTurbines = np.tile(Vinf,nTurbines)

        for turbI in range(0,nTurbines):

            # find overlap-area weighted effect of each wake zone
            wakeEffCoeff = 0
            for turb in range(0, nTurbines):

                wakeEffCoeffPerZone = 0
                deltax = turbineX[turbI] - turbineX[turb]

                if deltax > 0:
                    mU = MU / np.cos(aU*np.pi/180 + bU*yaw[turb])
                    for zone in range(0,3):
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + np.power((rotorDiameter[turb])/(rotorDiameter[turb]+2*ke[turb]*mU[zone]*deltax),2.0) * wakeOverlapTRel[turbI,turb,zone]

                    wakeEffCoeff = wakeEffCoeff + np.power(axialInd[turb]*wakeEffCoeffPerZone,2.0)

            wakeEffCoeff = (1 - 2 * np.sqrt(wakeEffCoeff))

            # multiply the inflow speed with the wake coefficients to find effective wind speed at turbine
            self.velocitiesTurbines[turbI] *= wakeEffCoeff
        if self.verbose:
            print "wind speed at turbines %s [m/s]" % self.velocitiesTurbines

        # find turbine powers
        self.wt_power = np.power(self.velocitiesTurbines,3.0) * (0.5*rho*rotorArea*Cp*np.power(np.cos(yaw),pP))
        if self.verbose:
            print "powers turbines %s [kW]" % self.wt_power

        # set outputs on turbine level
        for turbI in range(0, nTurbines):
            turbineName = self.wt_layout.wt_names[turbI]
            getattr(self.wt_layout, turbineName).power = self.wt_power[turbI] # in W
            getattr(self.wt_layout, turbineName).wind_speed_eff = self.velocitiesTurbines[turbI]

        self.wt_power /= 1000  # in kW
        self.power = np.sum(self.wt_power)


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
        Cp = np.hstack(self.wt_layout.wt_array(attr='CP'))
        axialInd = np.hstack(self.wt_layout.wt_array(attr='axial_induction'))

        positions = self.wt_layout.wt_array(attr='position')
        turbineX = positions[:, 0]
        turbineY = positions[:, 1]

        if self.ws_positions.any():
            velX = self.ws_positions[:, 0]
            velY = self.ws_positions[:, 1]
        else:
            velX = np.zeros([0,0])
            velY = np.zeros([0,0])

        pP = self.parameters.pP
        ke = self.parameters.ke
        keCorrDA = self.parameters.keCorrDA
        kd = self.parameters.kd
        me = self.parameters.me
        MU = self.parameters.MU
        ad = self.parameters.ad
        bd = self.parameters.bd
        aU = self.parameters.aU
        bU = self.parameters.bU

        # find size of arrays
        nTurbines = turbineX.size
        nLocations = velX.size

        # convert to downwind-crosswind coordinates
        rotationMatrix = np.array([(np.cos(-windDirection), -np.sin(-windDirection)),
                                   (np.sin(-windDirection), np.cos(-windDirection))])
        turbineLocations = np.dot(rotationMatrix,np.array([turbineX,turbineY]))
        turbineX = turbineLocations[0]
        turbineY = turbineLocations[1]

        if velX.size>0:
            locations = np.dot(rotationMatrix,np.array([velX,velY]))
            velX = locations[0]
            velY = locations[1]

        # calculate y-location of wake centers
        wakeCentersY = np.zeros((nLocations,nTurbines))
        wakeCentersYT = np.zeros((nTurbines,nTurbines))
        for turb in range(0,nTurbines):
            wakeAngleInit = 0.5*np.power(np.cos(yaw[turb]),2)*np.sin(yaw[turb])*4*axialInd[turb]*(1-axialInd[turb])
            for loc in range(0,nLocations):  # at velX-locations
                deltax = np.maximum(velX[loc]-turbineX[turb],0)
                factor = (2*kd*deltax/rotorDiameter[turb])+1
                wakeCentersY[loc,turb] = turbineY[turb]
                wakeCentersY[loc,turb] = wakeCentersY[loc,turb] + ad+bd*deltax  # rotation-induced deflection
                wakeCentersY[loc,turb] = wakeCentersY[loc,turb] + \
                                         (wakeAngleInit*(15*np.power(factor,4)+np.power(wakeAngleInit,2))/((30*kd*np.power(factor,5))/rotorDiameter[turb]))-(wakeAngleInit*rotorDiameter[turb]*(15+np.power(wakeAngleInit,4))/(30*kd)) # yaw-induced deflection
            for turbI in range(0,nTurbines):  # at turbineX-locations
                deltax = np.maximum(turbineX[turbI]-turbineX[turb],0)
                factor = (2*kd*deltax/rotorDiameter[turb])+1
                wakeCentersYT[turbI,turb] = turbineY[turb]
                wakeCentersYT[turbI,turb] = wakeCentersYT[turbI,turb] + ad+bd*deltax  # rotation-induced deflection
                wakeCentersYT[turbI,turb] = wakeCentersYT[turbI,turb] + \
                                            (wakeAngleInit*(15*np.power(factor,4)+np.power(wakeAngleInit,4))/((30*kd*np.power(factor,5))/rotorDiameter[turb]))-(wakeAngleInit*rotorDiameter[turb]*(15+np.power(wakeAngleInit,4))/(30*kd)) # yaw-induced deflection

        # calculate wake zone diameters at velX-locations
        wakeDiameters = np.zeros((nLocations,nTurbines,3))
        wakeDiametersT = np.zeros((nTurbines,nTurbines,3))
        for turb in range(0,nTurbines):
            for loc in range(0,nLocations):  # at velX-locations
                deltax = velX[loc]-turbineX[turb]
                for zone in range(0,3):
                    wakeDiameters[loc,turb,zone] = rotorDiameter[turb]+2*ke*me[zone]*np.maximum(deltax,0)
            for turbI in range(0,nTurbines):  # at turbineX-locations
                deltax = turbineX[turbI]-turbineX[turb]
                for zone in range(0,3):
                    wakeDiametersT[turbI,turb,zone] = np.maximum(rotorDiameter[turb]+2*ke*me[zone]*deltax,0)

        # calculate overlap areas at rotors
        # wakeOverlapT(TURBI,TURB,ZONEI) = overlap area of zone ZONEI of wake
        # of turbine TURB with rotor of turbine TURBI
        wakeOverlapT = calcOverlapAreas(turbineX,turbineY,rotorDiameter,wakeDiametersT,wakeCentersYT)

        # make overlap relative to rotor area (maximum value should be 1)
        wakeOverlapTRel = wakeOverlapT
        for turb in range(0,nTurbines):
            wakeOverlapTRel[turb] = wakeOverlapTRel[turb]/rotorArea[turb]

        # array effects with full or partial wake overlap:
        # use overlap area of zone 1+2 of upstream turbines to correct ke
        # Note: array effects only taken into account in calculating
        # velocity deficits, in order not to over-complicate code
        # (avoid loops in calculating overlaps)

        keUncorrected = ke
        ke = np.zeros(nTurbines)
        for turb in range(0,nTurbines):
            s = np.sum(wakeOverlapTRel[turb,:,0]+wakeOverlapTRel[turb,:,1])
            ke[turb] = keUncorrected*(1+s*keCorrDA)

        # calculate velocities in full flow field (optional)
        self.ws_array = np.tile(Vinf,nLocations)
        for turb in range(0,nTurbines):
            mU = MU/np.cos(aU*np.pi/180+bU*yaw[turb])
            for loc in range(0,nLocations):
                deltax = velX[loc] - turbineX[turb]
                radiusLoc = abs(velY[loc]-wakeCentersY[loc,turb])
                axialIndAndNearRotor = 2*axialInd[turb]

                if deltax > 0 and radiusLoc < wakeDiameters[loc,turb,0]/2.0:    # check if in zone 1
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*ke[turb]*(mU[0])*np.maximum(0,deltax))),2)
                elif deltax > 0 and radiusLoc < wakeDiameters[loc,turb,1]/2.0:    # check if in zone 2
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*ke[turb]*(mU[1])*np.maximum(0,deltax))),2)
                elif deltax > 0 and radiusLoc < wakeDiameters[loc,turb,2]/2.0:    # check if in zone 3
                    reductionFactor = axialIndAndNearRotor*\
                                      np.power((rotorDiameter[turb]/(rotorDiameter[turb]+2*ke[turb]*(mU[2])*np.maximum(0,deltax))),2)
                elif deltax <= 0 and radiusLoc < rotorDiameter[turb]/2.0:     # check if axial induction zone in front of rotor
                    reductionFactor = axialIndAndNearRotor*(0.5+np.arctan(2.0*np.minimum(0,deltax)/(rotorDiameter[turb]))/np.pi)
                else:
                    reductionFactor = 0
                self.ws_array[loc] *= (1-reductionFactor)

        # find effective wind speeds at downstream turbines, then predict power downstream turbine
        self.velocitiesTurbines = np.tile(Vinf,nTurbines)

        for turbI in range(0,nTurbines):

            # find overlap-area weighted effect of each wake zone
            wakeEffCoeff = 0
            for turb in range(0, nTurbines):

                wakeEffCoeffPerZone = 0
                deltax = turbineX[turbI] - turbineX[turb]

                if deltax > 0:
                    mU = MU / np.cos(aU*np.pi/180 + bU*yaw[turb])
                    for zone in range(0,3):
                        wakeEffCoeffPerZone = wakeEffCoeffPerZone + np.power((rotorDiameter[turb])/(rotorDiameter[turb]+2*ke[turb]*mU[zone]*deltax),2.0) * wakeOverlapTRel[turbI,turb,zone]

                    wakeEffCoeff = wakeEffCoeff + np.power(axialInd[turb]*wakeEffCoeffPerZone,2.0)

            wakeEffCoeff = (1 - 2 * np.sqrt(wakeEffCoeff))

            # multiply the inflow speed with the wake coefficients to find effective wind speed at turbine
            self.velocitiesTurbines[turbI] *= wakeEffCoeff

        # find turbine powers
        self.wt_power = np.power(self.velocitiesTurbines,3.0) * (0.5*rho*rotorArea*Cp*np.power(np.cos(yaw),pP))

        # set outputs on turbine level
        for turbI in range(0, nTurbines):
            turbineName = self.wt_layout.wt_names[turbI]
            getattr(self.wt_layout, turbineName).power = self.wt_power[turbI] # in W
            getattr(self.wt_layout, turbineName).wind_speed_eff = self.velocitiesTurbines[turbI]

        self.wt_power /= 1000  # in kW
        self.power = np.sum(self.wt_power)


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