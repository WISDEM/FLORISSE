from openmdao.main.api import Assembly
from floris import FLORIS
from florisRotorSE import CCBladeCoefficients, AeroelasticHAWTVT_floris, Controller, windSpeedDistributor
from openmdao.lib.drivers.api import FixedPointIterator
import os
import numpy as np

# Define flow properties
freestream_wind_speed = 12.0  # m/s
air_density = 1.1716  # kg/m^3
wind_direction = 30  # deg
viscosity = 1.81206e-5

# define turbine positions
turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1])
turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7])

# create FLORIS model of wind plant consisting of AeroelasticHAWTVT_floris turbines
florisWindPlant = FLORIS()
for turbineI in range(0,turbineX.size):
    wt = AeroelasticHAWTVT_floris()

    wt.name = 'turbine%s' % turbineI
    wt.position = np.array([turbineX[turbineI], turbineY[turbineI]])

    # AeroelasticHAWTVT inputs
    wt.turbine_name = 'NREL 5-MW baseline turbine'

    wt.rotor_diameter = 126.4
    wt.rotor_area = 0.25*np.pi*wt.rotor_diameter*wt.rotor_diameter
    wt.yaw = 0.0

    florisWindPlant.wt_layout.add_wt(wt)

# Set up assembly and FixedPointIterator solver for combined model FLORISplusCCBlade
FLORISplusCCBlade = Assembly()
FLORISplusCCBlade.add('driver', FixedPointIterator()) # Note that FixedPointIterator is top-level (driver)

# add FLORIS model to FLORISplusCCBlade Assembly
FLORISplusCCBlade.add('floris',florisWindPlant)

# add CCBlade model for each turbine to FLORISplusCCBlade Assembly
for turbineI in range(0,turbineX.size):
    ccbladeName = 'CCblade%s' % turbineI
    FLORISplusCCBlade.add(ccbladeName,CCBladeCoefficients())

# add controller for each turbine to FLORISplusCCBlade Assembly
for turbineI in range(0,turbineX.size):
    controllerName = 'controller%s' % turbineI
    FLORISplusCCBlade.add(controllerName,Controller())

# add windSpeedDistributor for each turbine to FLORISplusCCBlade Assembly
for turbineI in range(0,turbineX.size):
    distName = 'wsDist%s' % turbineI
    FLORISplusCCBlade.add(distName, windSpeedDistributor())

# define workflow for FixedPointIterator
workflow = []
for turbineI in range(0,turbineX.size):
    distName = 'wsDist%s' % turbineI
    workflow.append(distName)
for turbineI in range(0,turbineX.size):
    controllerName = 'controller%s' % turbineI
    workflow.append(controllerName)
for turbineI in range(0,turbineX.size):
    ccbladeName = 'CCblade%s' % turbineI
    workflow.append(ccbladeName)
workflow.append('floris')
print workflow

FLORISplusCCBlade.driver.workflow.add(workflow)

# connect CCblade-predicted axial induction and CP to FLORIS
for turbineI in range(0, turbineX.size):
    FLORISplusCCBlade.connect("CCblade%s.CP" % turbineI, "floris.wt_layout.turbine%s.CP" % turbineI)
    FLORISplusCCBlade.connect("CCblade%s.axial_induction" % turbineI, "floris.wt_layout.turbine%s.axial_induction" % turbineI)

# connect windSpeedDistributor wind_speed to controller and CCBlade
for turbineI in range(0, turbineX.size):
     FLORISplusCCBlade.connect("wsDist%s.wind_speed_turbine" % turbineI, "CCblade%s.turbine.wind_speed_hub" % turbineI)
     FLORISplusCCBlade.connect("wsDist%s.wind_speed_controller" % turbineI, "controller%s.turbine.wind_speed_hub" % turbineI)

# connect controller rotor_speed and pitch to CCBlade
for turbineI in range(0, turbineX.size):
    FLORISplusCCBlade.connect("controller%s.turbine.rotor_speed" % turbineI, "CCblade%s.turbine.rotor_speed" % turbineI)
    FLORISplusCCBlade.connect("controller%s.turbine.pitch" % turbineI, "CCblade%s.turbine.pitch" % turbineI)

# connect FLORIS-predicted wind_speed_eff to windSpeedDistributor through FixedPointIterator
FLORISplusCCBlade.driver.tolerance = .00001
FLORISplusCCBlade.driver.max_iteration = 100
for turbineI in range(0, turbineX.size):
    FLORISplusCCBlade.driver.add_parameter("wsDist%s.wind_speed" % turbineI, low = FLORISplusCCBlade.driver.tolerance, high = freestream_wind_speed, start = freestream_wind_speed)
    FLORISplusCCBlade.driver.add_constraint("floris.wt_layout.turbine%s.wind_speed_eff = wsDist%s.wind_speed" % (turbineI, turbineI))

# define rotor properties for CCBlade
for turbineI in range(0,turbineX.size):
    wt = getattr(FLORISplusCCBlade,'CCblade%s' % turbineI).turbine

    wt.name = 'turbine%s' % turbineI
    wt.position = np.array([turbineX[turbineI], turbineY[turbineI]])
    wt.turbine_name = 'NREL 5-MW baseline turbine'
    wt.rotor_diameter = 126.4
    wt.rotor_area = 0.25*np.pi*wt.rotor_diameter*wt.rotor_diameter # TODO: remove double definitions
    wt.tip_radius = 63.2 # TODO: remove double definitions
    wt.hub_radius = 1.5
    wt.nblades = 3  # number of blades
    wt.hub_height = 90.0
    wt.tilt_angle = 5.0
    wt.cone_angle = 2.5

    # CCblade inputs
    wt.nSector = 8
    wt.radial_locations = np.array([2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500, 28.1500,
                                    32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500, 56.1667, 58.9000, 61.6333])
    wt.chord = np.array([3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
                         3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
    wt.theta = np.array([13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
                         6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106])
    # define airfoil
    basePath = os.path.join(os.getcwd(), '5MW_AFFiles')
    airfoil_types = ['Cylinder1.dat','Cylinder2.dat','DU40_A17.dat','DU35_A17.dat','DU30_A17.dat','DU25_A17.dat','DU21_A17.dat','NACA64_A17.dat']
    af_idx = [0, 0, 1, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7, 7, 7, 7, 7] # place at appropriate radial stations
    wt.airfoil_files = [basePath]*len(wt.radial_locations)
    for i in range(len(wt.radial_locations)):
        wt.airfoil_files[i] = os.path.join(basePath,airfoil_types[af_idx[i]])
    wt.air_density = air_density
    wt.viscosity = viscosity

    # control DOFs
    wt.yaw = 0.0

    getattr(FLORISplusCCBlade,'CCblade%s' % turbineI).turbine = wt

# augment turbine with controller settings for Controller()

attributes = ("nblades", "hub_height", "tilt_angle", "cone_angle", "hub_radius", "tip_radius", "rotor_area",
              "air_density","dynamic_viscosity","shear_exponent","radial_locations","chord","theta","precurve",
              "precurveTip","airfoil_files","nSector","tiploss","hubloss","wakerotation","usecd","yaw")
for turbineI in range(0,turbineX.size):
    # copy CCBlade properties to turbine definition for controller
    wt = getattr(FLORISplusCCBlade, 'CCblade%s' % turbineI).turbine
    for attribute in attributes:
        setattr(getattr(FLORISplusCCBlade, 'controller%s' % turbineI).turbine, attribute, getattr(wt,attribute))
    # add control-relevant attributes
    getattr(FLORISplusCCBlade, 'controller%s' % turbineI).turbine.rated_generator_speed = 1173.7
    getattr(FLORISplusCCBlade, 'controller%s' % turbineI).turbine.gearbox_ratio = 97.0
    getattr(FLORISplusCCBlade, 'controller%s' % turbineI).turbine.rated_power = 5.0e6
    getattr(FLORISplusCCBlade, 'controller%s' % turbineI).turbine.generator_efficiency = 0.944
    getattr(FLORISplusCCBlade, 'controller%s' % turbineI).turbine.opt_tsr = 7.55
    getattr(FLORISplusCCBlade, 'controller%s' % turbineI).verbose = True
    getattr(FLORISplusCCBlade, 'controller%s' % turbineI).turbine.turbineName = wt.name

# Define flow properties in FLORIS
FLORISplusCCBlade.floris.wind_speed = freestream_wind_speed  # m/s
FLORISplusCCBlade.floris.air_density = air_density  # kg/m^3
FLORISplusCCBlade.floris.wind_direction = wind_direction  # deg

FLORISplusCCBlade.run()

# print the CP and axial induction of each turbine, and the rotor speed and pitch
for turbineI in range(0, turbineX.size):
    florisTurbine = "turbine%s" % turbineI
    print getattr(FLORISplusCCBlade,"controller%s" % turbineI).turbine.turbineName
    print " axial induction %s" % getattr(FLORISplusCCBlade.floris.wt_layout,florisTurbine).axial_induction
    print " C_P %s" % getattr(FLORISplusCCBlade.floris.wt_layout,florisTurbine).CP
    print " power %s W" % getattr(FLORISplusCCBlade.floris.wt_layout,florisTurbine).power
    wind_speed_floris = getattr(FLORISplusCCBlade.floris.wt_layout,florisTurbine).wind_speed_eff
    print " wind speed floris side %s m/s" % wind_speed_floris
    print " error with wind speed turbine side %s m/s" % (wind_speed_floris - getattr(FLORISplusCCBlade, 'wsDist%s' % turbineI).wind_speed)
    rotor_speed = getattr(FLORISplusCCBlade,"CCblade%s" % turbineI).turbine.rotor_speed
    tip_radius = getattr(FLORISplusCCBlade,"CCblade%s" % turbineI).turbine.tip_radius
    print " rotor speed %s RPM" % rotor_speed
    print " pitch %s deg" % getattr(FLORISplusCCBlade,"CCblade%s" % turbineI).turbine.pitch
    print " tip-speed ratio %s" % (tip_radius * np.pi/30.0 * rotor_speed / wind_speed_floris)

print "number of open-loop model executions: %s" % FLORISplusCCBlade.floris.exec_count