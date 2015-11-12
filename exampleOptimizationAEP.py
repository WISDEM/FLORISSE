from openmdao.api import Problem, pyOptSparseDriver
from OptimizationGroups import OptAEP

import time
import numpy as np


if __name__ == "__main__":

    # define turbine locations in global reference frame
    turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1])   # m
    turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7])   # m

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
        rotorDiameter[turbI] = 126.4      # m
        axialInduction[turbI] = 1.0/3.0
        Ct[turbI] = 4.0*axialInduction[turbI]*(1.0-axialInduction[turbI])
        Cp[turbI] = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
        generator_efficiency[turbI] = 0.944
        yaw[turbI] = 0.     # deg.

    # Define flow properties
    wind_speed = 8.0        # m/s
    air_density = 1.1716    # kg/m^3
    wind_direction = 240    # deg (N = 0 deg., using direction FROM, as in met-mast data)
    windDirections = np.array([wind_direction, wind_direction-90, wind_direction+90])

    # initialize problem
    prob = Problem(root=OptAEP(nTurbines=nTurbs, nDirections=windDirections.size, resolution=0))

    # set up optimizer
    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    prob.driver.add_objective('obj')

    # select design variables
    # prob.driver.add_desvar('turbineX', low=np.ones(nTurbs)*min(turbineX), high=np.ones(nTurbs)*max(turbineX))
    # prob.driver.add_desvar('turbineY', low=np.ones(nTurbs)*min(turbineY), high=np.ones(nTurbs)*max(turbineY))
    for i in range(0, windDirections.size):
        prob.driver.add_desvar('yaw%i' % i, low=-30.0, high=30.0)


    # initialize problem
    prob.setup()

    # assign initial values to design variables
    prob['turbineX'] = turbineX
    prob['turbineY'] = turbineY
    for i in range(0, windDirections.size):
        prob['yaw%i' % i] = yaw

    # assign values to constant inputs (not design variables)
    prob['rotorDiameter'] = rotorDiameter
    prob['axialInduction'] = axialInduction
    prob['generator_efficiency'] = generator_efficiency
    prob['wind_speed'] = wind_speed
    prob['air_density'] = air_density
    prob['windDirections'] = windDirections
    prob['Ct_in'] = Ct
    prob['Cp_in'] = Cp
    prob['floris_params:FLORISoriginal'] = True
    prob['floris_params:CPcorrected'] = False
    prob['floris_params:CTcorrected'] = False
    prob['windrose_frequencies'] = np.array([0.25, 0.5, 0.25])

    # run the problem
    print 'start FLORIS run'
    tic = time.time()
    prob.run()
    toc = time.time()

    # print the results
    # print the results
    print('FLORIS Opt. calculation took %.03f sec.' % (toc-tic))

    for i in range(0, windDirections.size):
        print 'yaw%i (deg) = ' % i, prob.root.unknowns['yaw%i' % i]
        exec("print 'velocities in dir%i: ', prob.root.AEPgroup.dir%i.unknowns['velocitiesTurbines']" % (i, i))

    print 'turbine X positions in wind frame (m): %s' % prob.root.unknowns['turbineX']
    print 'turbine Y positions in wind frame (m): %s' % prob.root.unknowns['turbineY']
    print 'power in each direction (kW): %s' % prob.root.AEPgroup.AEPcomp.params['power_directions']
    print 'AEP (kWh): %s' % prob.root.unknowns['AEP']