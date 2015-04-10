import numpy as np
from florisRotorSE2 import AeroelasticHAWTVT_CCBlade_floris, constructCoupledFLORIS_CCBlade_control
from fusedwind.plant_flow.vt import GenericWindFarmTurbineLayout
import os
import matplotlib.pyplot as plt

freestream_wind_speed = 12.0  # m/s
air_density = 1.1716  # kg/m^3
wind_direction = 30  # deg
viscosity = 1.81206e-5

# locations of turbines in Princess Amalia Wind Park, North Sea, The Netherlands
turbineX = np.array([972.91, 1552.6, 2157.5, 720.86,1290.5,1895.4,2515.5,3135.5,468.81,1033.4,1623.2,2238.2,2853.2,3483.3,216.76,776.31,1356,1955.9,2570.9,3185.9,3800.9,519.22,1093.9,1683.7,2283.6,2893.5,3503.5,257.09,821.68,1406.4,1996.2,2601.2,3201,3780.8,0,559.55,1129.2,1713.9,2308.8,2903.6,3468.2,287.34,851.93,1431.6,2006.3,2596.1,3155.7,3710.2,20.164,574.67,1139.3,1719,2283.6,2843.1,297.42,851.93,1426.6,1986.2,1124.1,1683.7])
turbineY = np.array([4879.69, 4889.77, 4869.61, 4380.63,4395.75,4380.63,4350.38,4310.06,3886.61,3911.82,3901.73,3866.45,3831.16,3785.79,3392.59,3417.80,3412.76,3392.59,3357.31,3316.98,3256.49,2928.82,2928.82,2903.62,2883.45,2843.12,2792.71,2429.76,2439.84,2429.76,2404.56,2369.27,2323.90,2278.53,1940.78,1955.91,1945.83,1935.74,1900.46,1865.17,1819.80,1466.93,1461.89,1451.81,1431.64,1396.36,1356.03,1320.74,977.95,977.95,977.95,962.83,927.54,897.30,499.06,504.10,488.98,468.81,20.16,0.00])

#turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1])
#turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7])

wt_layout = GenericWindFarmTurbineLayout()

for turbineI in range(0,turbineX.size):
    wt = AeroelasticHAWTVT_CCBlade_floris()

    wt.name = 'turbine%s' % turbineI
    wt.position = np.array([turbineX[turbineI], turbineY[turbineI]])
    wt.turbine_name = 'NREL 5-MW baseline turbine'
    wt.wind_speed_hub = freestream_wind_speed
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
    wt.verbose = True
    wt.turbineName = wt.name  # double names needed, that is to be fixed! Related to name being overwritten by signal name (e.g. turbineIn)

    # control DOFs
    wt.yaw = 0.0
    wt_layout.add_wt(wt)

FLORIS_CCBlade_control = constructCoupledFLORIS_CCBlade_control(len(wt_layout.wt_list), startpoint_wind_speed=freestream_wind_speed, tolerance=0.00001, max_iteration=100)

FLORIS_CCBlade_control.control.listOfTurbinesIn = wt_layout.wt_list
FLORIS_CCBlade_control.floris.wind_speed = freestream_wind_speed  # m/s
FLORIS_CCBlade_control.floris.air_density = air_density  # kg/m^3
FLORIS_CCBlade_control.floris.wind_direction = wind_direction  # deg
FLORIS_CCBlade_control.control.verbose = True

FLORIS_CCBlade_control.run()

print "number of iterations: %s" % FLORIS_CCBlade_control.floris.exec_count

power = np.zeros(len(wt_layout.wt_list))

for turbineI in range(0, turbineX.size):
    florisTurbine = "turbine%s" % turbineI
    print FLORIS_CCBlade_control.ccblade.listOfTurbinesOut[turbineI].turbineName
    print " axial induction %s" % FLORIS_CCBlade_control.floris.florisWindPlant.wt_layout.wt_list[turbineI].axial_induction
    print " C_P %s" % FLORIS_CCBlade_control.floris.florisWindPlant.wt_layout.wt_list[turbineI].CP
    power[turbineI] = FLORIS_CCBlade_control.floris.florisWindPlant.wt_layout.wt_list[turbineI].power*1e-6*FLORIS_CCBlade_control.control.listOfTurbinesIn[turbineI].generator_efficiency
    print " power %s MW" % power[turbineI]
    wind_speed_floris = FLORIS_CCBlade_control.floris.florisWindPlant.wt_layout.wt_list[turbineI].wind_speed_eff
    print " wind speed floris side %s m/s" % wind_speed_floris
    print " error with wind speed turbine side %s m/s" % (wind_speed_floris - FLORIS_CCBlade_control.windSpeedDistributor.wind_speed_in[turbineI])
    rotor_speed = FLORIS_CCBlade_control.ccblade.listOfTurbinesOut[turbineI].rotor_speed
    tip_radius = FLORIS_CCBlade_control.ccblade.listOfTurbinesOut[turbineI].tip_radius
    print " rotor speed %s RPM" % rotor_speed
    print " pitch %s deg" % FLORIS_CCBlade_control.ccblade.listOfTurbinesOut[turbineI].pitch
    print " tip-speed ratio %s" % (tip_radius * np.pi/30.0 * rotor_speed / wind_speed_floris)

plt.scatter(turbineX,turbineY,c=power,s=50,cmap='bwr')
plt.axis('equal')
cb = plt.colorbar()
cb.set_label('turbine generator power (MW)')
plt.xlabel('distance (m)')
plt.ylabel('distance (m)')
for turbineI in range(0, turbineX.size):
    plt.text(turbineX[turbineI],turbineY[turbineI]+20,"%.2f MW" % power[turbineI], fontsize=8, horizontalalignment='center')

plt.show()

