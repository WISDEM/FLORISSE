from openmdao.api import Problem, Group, IndepVarComp, ScipyOptimizer
from floris_openmdao1 import FLORIS, adjustCtCp_yaw

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

    # set up problem
    prob = Problem(root=Group())
    prob.root.add('myFloris', FLORIS(nTurbs, resolution=0))
    prob.root.add('CtCp', adjustCtCp_yaw(nTurbs))

    # initialize design variables for optimization
    prob.root.add('p1', IndepVarComp('turbineX', turbineX))
    prob.root.add('p2', IndepVarComp('turbineY', turbineY))

    # connect design variables to corresponding components
    prob.root.connect('p1.turbineX', 'myFloris.turbineX')
    prob.root.connect('p2.turbineY', 'myFloris.turbineY')

    # set up optimizer
    prob.driver = ScipyOptimizer()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.driver.options['tol'] = 1.0E-8
    prob.driver.add_desvar('p1.turbineX', low=np.ones(nTurbs)*min(turbineX), high=np.ones(nTurbs)*max(turbineX))
    prob.driver.add_desvar('p2.turbineY', low=np.ones(nTurbs)*min(turbineY), high=np.ones(nTurbs)*max(turbineY))
    prob.driver.add_objective('myFloris.power')

    # connect components within the problem
    prob.root.connect('CtCp.Ct_out', 'myFloris.Ct')
    prob.root.connect('CtCp.Cp_out', 'myFloris.Cp')
    prob.root.connect('CtCp.yaw', 'myFloris.yaw')

    # initialize problem
    prob.setup()

    # assign values to constant inputs (not design variables)
    prob['myFloris.rotorDiameter'] = rotorDiameter
    prob['myFloris.axialInduction'] = axialInduction
    prob['myFloris.generator_efficiency'] = generator_efficiency
    prob['myFloris.wind_speed'] = wind_speed
    prob['myFloris.air_density'] = air_density
    prob['myFloris.wind_direction'] = wind_direction
    prob['CtCp.yaw'] = yaw
    prob['CtCp.Ct_in'] = Ct
    prob['CtCp.Cp_in'] = Cp
    prob['myFloris.floris_params:FLORISoriginal'] = True
    prob['CtCp.floris_params:FLORISoriginal'] = True

    # run the problem
    print 'start FLORIS run'
    tic = time.time()
    prob.run()
    toc = time.time()

    # print the results
    print('FLORIS Opt. calculation took %.03f sec.' % (toc-tic))
    print 'turbine powers (kW): %s' % prob.root.myFloris.unknowns['wt_power']
    print 'turbine X positions in wind frame (m): %s' % prob.root.myFloris.f_1.params['turbineX']
    print 'turbine Y positions in wind frame (m): %s' % prob.root.myFloris.f_1.params['turbineY']
    print 'yaw (deg) = ', prob.root.myFloris.f_4.params['yaw']
    print 'effective wind speeds (m/s): %s' % prob.root.myFloris.f_4.unknowns['velocitiesTurbines']
    print 'wind farm power (kW): %s' % prob.root.myFloris.f_4.unknowns['power']

