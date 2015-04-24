import numpy as np
from florisRotorSE import AeroelasticHAWTVT_CCBlade_floris, CoupledFLORIS_CCBlade_control
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout
from fusedwind.plant_flow.generate_fake_vt import generate_random_GenericWindRoseVT
from fusedwind.plant_flow.asym import AEPSingleWindRose
import os
from copy import copy
import matplotlib.pyplot as plt

wind_rose = generate_random_GenericWindRoseVT()
air_density = 1.1716  # kg/m^3
viscosity = 1.81206e-5




# Define the turbine

wt = AeroelasticHAWTVT_CCBlade_floris()

wt.turbine_name = 'NREL 5-MW baseline turbine'
wt.rotor_diameter = 126.4
wt.rotor_area = 0.25*np.pi*wt.rotor_diameter*wt.rotor_diameter # TODO: remove double definitions
wt.tip_radius = 63.2 # TODO: remove double definitions
wt.hub_radius = 1.5
wt.nblades = 3  # number of blades
wt.hub_height = 90.0
wt.tilt_angle = 5.0
wt.cone_angle = 2.5

wt.yaw = 0.0

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

# control properties
wt.rated_generator_speed = 1173.7
wt.gearbox_ratio = 97.0
wt.rated_power = 5.0e6
wt.generator_efficiency = 0.944
wt.opt_tsr = 7.55
wt.cut_in_wind_speed = 3.0
wt.cut_out_wind_speed = 25.0
wt.verbose = True
wt.turbineName = 'NREL 5MW'  # double names needed, that is to be fixed! Related to name being overwritten by signal name (e.g. turbineIn)
wt.minimum_generator_speed = 670.0 #rpm
wt.transitional_generator_speed =  871.0 #rpm - transitional rotor speed between region 1.5 and 2

# control DOFs
wt.yaw = 0.0

# pre-calculate controller
print "============== PRE-CALCULATING CONTROLLER ==================="
wt.preCalculateController(visual=True)
print "============================================================="






# Add turbines to wind plant

# locations of turbines in Princess Amalia Wind Park, North Sea, The Netherlands
# turbineX = np.array([972.91, 1552.6, 2157.5, 720.86,1290.5,1895.4,2515.5,3135.5,468.81,1033.4,1623.2,2238.2,2853.2,3483.3,216.76,776.31,1356,1955.9,2570.9,3185.9,3800.9,519.22,1093.9,1683.7,2283.6,2893.5,3503.5,257.09,821.68,1406.4,1996.2,2601.2,3201,3780.8,0,559.55,1129.2,1713.9,2308.8,2903.6,3468.2,287.34,851.93,1431.6,2006.3,2596.1,3155.7,3710.2,20.164,574.67,1139.3,1719,2283.6,2843.1,297.42,851.93,1426.6,1986.2,1124.1,1683.7])
# turbineY = np.array([4879.69, 4889.77, 4869.61, 4380.63,4395.75,4380.63,4350.38,4310.06,3886.61,3911.82,3901.73,3866.45,3831.16,3785.79,3392.59,3417.80,3412.76,3392.59,3357.31,3316.98,3256.49,2928.82,2928.82,2903.62,2883.45,2843.12,2792.71,2429.76,2439.84,2429.76,2404.56,2369.27,2323.90,2278.53,1940.78,1955.91,1945.83,1935.74,1900.46,1865.17,1819.80,1466.93,1461.89,1451.81,1431.64,1396.36,1356.03,1320.74,977.95,977.95,977.95,962.83,927.54,897.30,499.06,504.10,488.98,468.81,20.16,0.00])

turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1])
turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7])

numberOfTurbines = turbineX.size

wt_layout = GenericWindFarmTurbineLayout()
for turbineI in range(0, numberOfTurbines):
   wt.name = 'turbine%s' % turbineI
   wt.turbineName = wt.name
   wt.position = np.array([turbineX[turbineI], turbineY[turbineI]])
   wt_layout.add_wt(copy(wt))


# TODO: prevent that we have to run before adding to AEP calculation
FLORIS_CCBlade_control = CoupledFLORIS_CCBlade_control()
FLORIS_CCBlade_control.coupledModel_verbosity = True
FLORIS_CCBlade_control.wt_layout = wt_layout
FLORIS_CCBlade_control.wind_speed = 30.0  # m/s
FLORIS_CCBlade_control.air_density = air_density  # kg/m^3
FLORIS_CCBlade_control.wind_direction = 0.0  # deg
FLORIS_CCBlade_control.controller_verbosity = True
FLORIS_CCBlade_control.run()






# AEP calculation
FLORIS_CCBlade_control.controller_verbosity = False
FLORIS_CCBlade_control.floris_verbosity = False
FLORIS_CCBlade_control.coupledModel_verbosity = True

aep = AEPSingleWindRose()
aep.add('wf', FLORIS_CCBlade_control)
aep.wind_speeds = wind_rose.wind_speeds
aep.wind_directions = wind_rose.wind_directions
aep.wind_rose = wind_rose.frequency_array
aep.configure()
aep.run()

print 'The energy production per sector [nWD, nWS]:'
print aep.array_aep

print 'Net Annual Energy Production after availability and loss impacts'
print aep.net_aep

# # TODO: get meaningful number out of capacity factor
# #print 'Capacity factor'
# #print aep.capacity_factor









