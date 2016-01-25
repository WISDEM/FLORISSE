from openmdao.api import Problem, Group, IndepVarComp
from floris_openmdao1 import DirectionGroupFLORIS
from GeneralWindfarmComponents import AdjustCtCpYaw

import time
import numpy as np


if __name__ == "__main__":

    # define turbine locations in global reference frame
    turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1])
    turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7])

    # initialize input variable arrays
    nTurbs = turbineX.size
    rotorDiameter = np.zeros(nTurbs)
    axialInduction = np.zeros(nTurbs)
    Ct = np.zeros(nTurbs)
    Cp = np.zeros(nTurbs)
    generator_efficiency = np.zeros(nTurbs)
    yaw = np.zeros(nTurbs)

    # define initial values
    for turbI in range(0, nTurbs):
        rotorDiameter[turbI] = 126.4            # m
        axialInduction[turbI] = 1.0/3.0
        Ct[turbI] = 4.0*axialInduction[turbI]*(1.0-axialInduction[turbI])
        Cp[turbI] = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
        generator_efficiency[turbI] = 0.944
        yaw[turbI] = 0.     # deg.

    # Define flow properties
    wind_speed = 8.0        # m/s
    air_density = 1.1716    # kg/m^3
    wind_direction = 240    # deg (N = 0 deg., using direction FROM, as in met-mast data)

    # set up problem
    prob = Problem(root=Group())
    prob.root.add('FLORIS', DirectionGroupFLORIS(nTurbs, resolution=0), promotes=['*'])

    # initialize problem
    prob.setup()

    # assign values to turbine states
    prob['turbineX'] = turbineX
    prob['turbineY'] = turbineY
    prob['yaw'] = yaw

    # assign values to constant inputs (not design variables)
    prob['rotorDiameter'] = rotorDiameter
    prob['axialInduction'] = axialInduction
    prob['generator_efficiency'] = generator_efficiency
    prob['wind_speed'] = wind_speed
    prob['air_density'] = air_density
    prob['wind_direction'] = wind_direction
    prob['Ct_in'] = Ct
    prob['Cp_in'] = Cp
    prob['floris_params:FLORISoriginal'] = False

    # run the problem
    print 'start FLORIS run'
    tic = time.time()
    prob.run()
    toc = time.time()

    # print the results
    print 'FLORIS calculation took %.06f sec.' % (toc-tic)
    print 'turbine X positions in wind frame (m): %s' % prob['turbineX']
    print 'turbine Y positions in wind frame (m): %s' % prob['turbineY']
    print 'yaw (deg) = ', prob['yaw']
    print 'Effective hub velocities (m/s) = ', prob['velocitiesTurbines0']
    print 'Turbine powers (kW) = ', prob['wt_power0']
    print 'wind farm power (kW): %s' % prob['power0']