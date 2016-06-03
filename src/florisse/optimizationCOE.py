from __future__ import print_function

from openmdao.api import Problem, ExecComp, pyOptSparseDriver, ScipyGMRES, IndepVarComp
#from florisse.OptimizationGroups import OptAEP
from florisse import config
from florisse.GeneralWindFarmComponents import SpacingComp
from COE import *

import time
import numpy as np
import pylab as plt

import cProfile


import sys

if __name__ == "__main__":
    """
    #Load in Amalia locations
    filename = "layout_amalia.txt"

    amalia = open(filename)
    x_y = np.loadtxt(amalia)
    turbineX = x_y[:,0]
    turbineY = x_y[:,1]
    
    #Amalia Data
    ##############################################################################
    windSpeeds = np.array([6.53163342, 6.11908394, 6.13415514, 6.0614625,  6.21344602,
                                5.87000793, 5.62161519, 5.96779107, 6.33589422, 6.4668016,
                                7.9854581,  7.6894432,  7.5089221,  7.48638098, 7.65764618,
                                6.82414044, 6.36728201, 5.95982999, 6.05942132, 6.1176321,
                                5.50987893, 4.18461796, 4.82863115, 0.,         0.,         0.,
                                5.94115843, 5.94914252, 5.59386528, 6.42332524, 7.67904937,
                                7.89618066, 8.84560463, 8.51601497, 8.40826823, 7.89479475,
                                7.86194762, 7.9242645,  8.56269962, 8.94563889, 9.82636368,
                               10.11153102, 9.71402212, 9.95233636,  10.35446959, 9.67156182,
                                9.62462527, 8.83545158, 8.18011771, 7.9372492,  7.68726143,
                                7.88134508, 7.31394723, 7.01839896, 6.82858346, 7.06213432,
                                7.01949894, 7.00575122, 7.78735165, 7.52836352, 7.21392201,
                                7.4356621,  7.54099962, 7.61335262, 7.90293531, 7.16021596,
                                7.19617087, 7.5593657,  7.03278586, 6.76105501, 6.48004694,
                                6.94716392])

    windFrequencies = np.array([1.17812570e-02, 1.09958570e-02, 9.60626600e-03, 1.21236860e-02,
                               1.04722450e-02, 1.00695140e-02, 9.68687400e-03, 1.00090550e-02,
                               1.03715390e-02, 1.12172280e-02, 1.52249700e-02, 1.56279300e-02,
                               1.57488780e-02, 1.70577560e-02, 1.93535770e-02, 1.41980570e-02,
                               1.20632100e-02, 1.20229000e-02, 1.32111160e-02, 1.74605400e-02,
                               1.72994400e-02, 1.43993790e-02, 7.87436000e-03, 0.00000000e+00,
                               2.01390000e-05, 0.00000000e+00, 3.42360000e-04, 3.56458900e-03,
                               7.18957000e-03, 8.80068000e-03, 1.13583200e-02, 1.41576700e-02,
                               1.66951900e-02, 1.63125500e-02, 1.31709000e-02, 1.09153300e-02,
                               9.48553000e-03, 1.01097900e-02, 1.18819700e-02, 1.26069900e-02,
                               1.58895900e-02, 1.77021600e-02, 2.04208100e-02, 2.27972500e-02,
                               2.95438600e-02, 3.02891700e-02, 2.69861000e-02, 2.21527500e-02,
                               2.12465500e-02, 1.82861400e-02, 1.66147400e-02, 1.90111800e-02,
                               1.90514500e-02, 1.63932050e-02, 1.76215200e-02, 1.65341460e-02,
                               1.44597600e-02, 1.40370300e-02, 1.65745000e-02, 1.56278200e-02,
                               1.53459200e-02, 1.75210100e-02, 1.59702700e-02, 1.51041500e-02,
                               1.45201100e-02, 1.34527800e-02, 1.47819600e-02, 1.33923300e-02,
                               1.10562900e-02, 1.04521380e-02, 1.16201970e-02, 1.10562700e-02])

    index = np.where(windSpeeds==0.0)
    windSpeeds = np.delete(windSpeeds, index[0])
    windFrequencies = np.delete(windFrequencies, index[0])
    print(len(windSpeeds))
    print(len(windFrequencies))

    nDirections = len(windSpeeds)

    ############################################################################################
    """
    config.floris_single_component = True

    
    ######################### for MPI functionality #########################
    from openmdao.core.mpi_wrap import MPI

    if MPI: # pragma: no cover
        # if you called this script with 'mpirun', then use the petsc data passing
        from openmdao.core.petsc_impl import PetscImpl as impl

    else:
        # if you didn't use 'mpirun', then use the numpy data passing
        from openmdao.api import BasicImpl as impl

    def mpi_print(prob, *args): 
        """helper function to only print on rank 0 """
        if prob.root.comm.rank == 0:
            print(*args)

    prob = Problem(impl=impl)

    size = 8 # number of processors (and number of wind directions to run)

    #########################################################################
    # define turbine size
    rotor_diameter = 126.4  # (m)

    # define turbine locations in global reference frame
    # original example case
    # turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1])   # m
    # turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7])   # m

    # Scaling grid case
    # nRows = int(sys.argv[1])     # number of rows and columns in grid
    nRows = 5
    spacing = 3    # turbine grid spacing in diameters

    # Set up position arrays
    points = np.linspace(start=spacing*rotor_diameter, stop=nRows*spacing*rotor_diameter, num=nRows)
    xpoints, ypoints = np.meshgrid(points, points)
    """
    turbineX = np.ndarray.flatten(xpoints)
    turbineY = np.ndarray.flatten(ypoints)
    """
    """
    turbineX1 = np.array([379.2,1137.6,1896,758.4,1516.8,379.2,1137.6,1896,758.4,1516.8,379.2,1137.6,1896])
    turbineX2 = np.array([758.4,1516.8,379.2,1137.6,1896,758.4,1516.8,379.2,1137.6,1896,758.4,1516.8])
    turbineX = np.hstack([turbineX1,turbineX2])
    turbineY1 = np.array([379.2,379.2,379.2,758.4,758.4,1137.6,1137.6,1137.6,1516.8,1516.8,1896,1896,1896])
    turbineY2 = np.array([379.2,379.2,758.4,758.4,758.4,1137.6,1137.6,1516.8,1516.8,1516.8,1896,1896])
    turbineY = np.hstack([turbineY1,turbineY2])
    """

    
    inverted = "no"
    

    nRows = 5
    spacing = 3   # turbine grid spacing in diameters
    nTurbs = nRows*nRows
    # Set up position arrays
    rotor_diameter = 126.4 
    points = np.linspace(start=spacing*rotor_diameter, stop=nRows*spacing*rotor_diameter, num=nRows)
    xpoints, ypoints = np.meshgrid(points, points)
    
    turbineX = np.ndarray.flatten(xpoints)
    turbineY = np.ndarray.flatten(ypoints)
    turbineX1 = np.array([])
    turbineX2 = np.array([])
    turbineY1 = np.array([])
    turbineY2 = np.array([])

    for i in range(nTurbs):
        if i%2 == 0:
            turbineX1 = np.append(turbineX1, turbineX[i])
            turbineY1 = np.append(turbineY1, turbineY[i])
        else:
            turbineX2 = np.append(turbineX2, turbineX[i])
            turbineY2 = np.append(turbineY2, turbineY[i])

    if inverted == "yes":
        turbineX = np.hstack([turbineX2,turbineX1])
        turbineY = np.hstack([turbineY2,yturbineY1])
    else:
        turbineX = np.hstack([turbineX1,turbineX2])
        turbineY = np.hstack([turbineY1,turbineY2])



    print('Turbine X', turbineX)
    print('Turbine Y', turbineY)
    # initialize input variable arrays
    nTurbs = turbineX.size
    #87.6m is the height of the 5MW NREL turbines
    turbineMin = 70
    turbineH1 = 130
    turbineH2 = 110
    #turbineH2 = 87.6
    #nTurbsH1 = 5
    nTurbsH1 = len(turbineX1)
    nTurbsH2 = nTurbs-nTurbsH1
    rotorDiameter = np.zeros(nTurbs)
    axialInduction = np.zeros(nTurbs)
    Ct = np.zeros(nTurbs)
    Cp = np.zeros(nTurbs)
    generatorEfficiency = np.zeros(nTurbs)
    yaw = np.zeros(nTurbs)
    minSpacing = 2.                         # number of rotor diameters

    # define initial values
    for turbI in range(0, nTurbs):
        rotorDiameter[turbI] = rotor_diameter      # m
        axialInduction[turbI] = 1.0/3.0
        Ct[turbI] = 4.0*axialInduction[turbI]*(1.0-axialInduction[turbI])
        Cp[turbI] = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
        generatorEfficiency[turbI] = 0.944
        yaw[turbI] = 0.     # deg.

    # Define flow properties
    wind_speed = 8.0        # m/s
    air_density = 1.1716    # kg/m^3
    nDirections = size
    windDirections = np.linspace(0, 360.-360/nDirections, nDirections)
    windFrequencies = np.ones_like(windDirections)*1.0/size

    # initialize problem

    # set up problem
    # prob = Problem(impl=impl, root=COE(nTurbs))
    prob = Problem()
    root = prob.root = Group()
    
    #TODO How to I pass TurbineZ down to calculate AEP? Basically, How do I get AEP?
    root.add('dv5', IndepVarComp('turbineH1', 0., units='m'), promotes=['*'])
    root.add('dv6', IndepVarComp('turbineH2', 0., units='m'), promotes=['*'])
    root.add('getTurbineZ', getTurbineZ(nTurbs), promotes=['*'])
    root.add('AEPGroup', AEPGroup(nTurbs, nDirections=size,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['*'])
    root.add('COEComponent', COEComponent(nTurbs), promotes=['*'])
    
    root.add('spacing_comp', SpacingComp(nTurbines=nTurbs), promotes=['*'])

    # add constraint definitions
    root.add('spacing_con', ExecComp('sc = wtSeparationSquared-(minSpacing*rotorDiameter[0])**2',
                                 minSpacing=minSpacing, rotorDiameter=np.zeros(nTurbs),
                                 sc=np.zeros(((nTurbs-1.)*nTurbs/2.)),
                                 wtSeparationSquared=np.zeros(((nTurbs-1.)*nTurbs/2.))),
         promotes=['*'])
    root.fd_options['force_fd'] = True
    #root.ln_solver = ScipyGMRES()
    

    # set up optimizer
    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    prob.driver.add_objective('COE', scaler=1E-2) #TODO???? COE is the objective ya?

    # set optimizer options
    prob.driver.opt_settings['Verify level'] = 3
    prob.driver.opt_settings['Print file'] = 'SNOPT_print_exampleOptAEP.out'
    prob.driver.opt_settings['Summary file'] = 'SNOPT_summary_exampleOptAEP.out'
    prob.driver.opt_settings['Major iterations limit'] = 1000
     

    # select design variables
    print('nTurbs: ', nTurbs)
    print('nTurbsH1: ', nTurbsH1)
    print('nturbsH2: ', nTurbsH2)
    prob.driver.add_desvar('turbineX', lower=np.ones(nTurbs)*min(turbineX), upper=np.ones(nTurbs)*max(turbineX), scaler=1E-2)
    prob.driver.add_desvar('turbineY', lower=np.ones(nTurbs)*min(turbineY), upper=np.ones(nTurbs)*max(turbineY), scaler=1E-2)
    prob.driver.add_desvar('turbineH1', turbineMin, upper=175, scaler=1E-2)
    prob.driver.add_desvar('turbineH2', turbineMin, upper=175, scaler=1E-2)
    # prob.driver.add_desvar('turbineZ', lower=np.ones(nTurbs)*100., upper=np.ones(nTurbs)*500., scaler=1E-2)
    """for direction_id in range(0, windDirections.size):
        prob.driver.add_desvar('yaw%i' % direction_id, lower=-30.0, upper=30.0, scaler=1E-1)"""
    

    # add constraints
    prob.driver.add_constraint('sc', lower=np.zeros(((nTurbs-1.)*nTurbs/2.)), scaler=1.0/rotor_diameter)

    tic = time.time()
    prob.setup(check=True)
    toc = time.time()

    # print the results
    mpi_print(prob, ('FLORIS setup took %.03f sec.' % (toc-tic)))

    # time.sleep(10)
    # assign initial values to design variables
    prob['turbineX'] = turbineX
    prob['turbineY'] = turbineY
    prob['turbineH1'] = turbineH1
    prob['turbineH2'] = turbineH2
    prob['nTurbsH1'] = nTurbsH1
    prob['nTurbsH2'] = nTurbsH2
    for direction_id in range(0, windDirections.size):
        prob['yaw%i' % direction_id] = yaw

    # assign values to constant inputs (not design variables)
    prob['rotorDiameter'] = rotorDiameter
    prob['axialInduction'] = axialInduction
    prob['generatorEfficiency'] = generatorEfficiency
    prob['windSpeeds'] = np.ones(nDirections)*wind_speed
    prob['air_density'] = air_density
    prob['windDirections'] = windDirections
    prob['windFrequencies'] = windFrequencies
    prob['Ct_in'] = Ct
    prob['Cp_in'] = Cp
    prob['floris_params:cos_spread'] = 1E12

    # set options
    # prob['floris_params:FLORISoriginal'] = True
    # prob['floris_params:CPcorrected'] = False
    # prob['floris_params:CTcorrected'] = False

    # run the problem
    mpi_print(prob, 'start FLORIS run')
    tic = time.time()
    # cProfile.run('prob.run()')
    prob.run()
    toc = time.time()
    #self.J = prob.check_total_derivatives()

    # print the results
    mpi_print(prob, ('FLORIS Opt. calculation took %.03f sec.' % (toc-tic)))

    """for direction_id in range(0, windDirections.size):
        mpi_print(prob,  'yaw%i (deg) = ' % direction_id, prob['yaw%i' % direction_id])"""
    # for direction_id in range(0, windDirections.size):
        # mpi_print(prob,  'velocitiesTurbines%i (m/s) = ' % direction_id, prob['velocitiesTurbines%i' % direction_id])
    # for direction_id in range(0, windDirections.size):
    #     mpi_print(prob,  'wt_power%i (kW) = ' % direction_id, prob['wt_power%i' % direction_id])

    mpi_print(prob,  'turbine X positions in wind frame (m): %s' % prob['turbineX'])
    mpi_print(prob,  'turbine Y positions in wind frame (m): %s' % prob['turbineY'])
    mpi_print(prob,  'turbine Height 1 (m): %s' % prob['turbineH1'])
    mpi_print(prob,  'turbine Height 2 (m): %s' % prob['turbineH2'])
    #mpi_print(prob,  'wind farm power in each direction (kW): %s' % prob['dirPowers'])
    #mpi_print(prob,  'AEP: %s' % prob['AEP'])
    mpi_print(prob,  'COE: %s' % prob['COE'])
    mpi_print(prob,  'AEP: %s' % prob['AEP'])

    optCOE = prob['COE']
    optAEP = prob['AEP']

    np.savetxt('XYZ5test.txt',np.c_[prob['turbineX'],prob['turbineY'],prob['turbineZ']])

    xbounds = [min(turbineX), min(turbineX), max(turbineX), max(turbineX), min(turbineX)]
    ybounds = [min(turbineY), max(turbineY), max(turbineY), min(turbineY), min(turbineX)]

    plt.figure(1)
    #plt.plot(turbineX, turbineY, 'ok', label='Original')
    plt.plot(turbineX1, turbineY1, '.k', label='Short Original')
    for i in range(0, nTurbsH1):
        plt.plot(prob['turbineX'][i],prob['turbineY'][i], 'ok')
    plt.plot(turbineX2, turbineY2, '.g', label='Tall Original')
    for i in range(nTurbsH1, nTurbsH1+nTurbsH2):
        plt.plot(prob['turbineX'][i],prob['turbineY'][i], 'og')
    #plt.plot(prob['turbineX'], prob['turbineY'], 'og', label='Optimized')
    plt.plot(xbounds, ybounds, ':k')
    for i in range(0, nTurbs):
        plt.plot([turbineX[i], prob['turbineX'][i]], [turbineY[i], prob['turbineY'][i]], '--k')
    plt.legend()
    plt.xlabel('Turbine X Position (m)')
    plt.ylabel('Turbine Y Position (m)')
    plt.legend(bbox_to_anchor=(1.14, 1.14))
    plt.show() 

    ##########################################################################################
    ##########################################################################################
    ##########################################################################################
    #Only modify x and y


    turbineH1 = turbineMin
    turbineH2 = turbineMin

    # set up problem
    prob = Problem()
    root = prob.root = Group()
    
    #TODO How to I pass TurbineZ down to calculate AEP? Basically, How do I get AEP?
    root.add('dv5', IndepVarComp('turbineH1', 0., units='m'), promotes=['*'])
    root.add('dv6', IndepVarComp('turbineH2', 0., units='m'), promotes=['*'])
    root.add('getTurbineZ', getTurbineZ(nTurbs), promotes=['*'])
    root.add('AEPGroup', AEPGroup(nTurbs, nDirections=size,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['*'])
    root.add('AEPobj', AEPobj(), promotes=['*'])
    root.add('COEComponent', COEComponent(nTurbs), promotes=['*'])
    root.add('spacing_comp', SpacingComp(nTurbines=nTurbs), promotes=['*'])

    # add constraint definitions
    root.add('spacing_con', ExecComp('sc = wtSeparationSquared-(minSpacing*rotorDiameter[0])**2',
                                 minSpacing=minSpacing, rotorDiameter=np.zeros(nTurbs),
                                 sc=np.zeros(((nTurbs-1.)*nTurbs/2.)),
                                 wtSeparationSquared=np.zeros(((nTurbs-1.)*nTurbs/2.))),
         promotes=['*'])
    root.fd_options['force_fd'] = True    

    # set up optimizer
    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'SNOPT'
    prob.driver.add_objective('maxAEP', scaler=1E-6)

    # set optimizer options
    prob.driver.opt_settings['Verify level'] = 3
    prob.driver.opt_settings['Print file'] = 'SNOPT_print_exampleOptAEP.out'
    prob.driver.opt_settings['Summary file'] = 'SNOPT_summary_exampleOptAEP.out'
    prob.driver.opt_settings['Major iterations limit'] = 1000
     

    # select design variables
    prob.driver.add_desvar('turbineX', lower=np.ones(nTurbs)*min(turbineX), upper=np.ones(nTurbs)*max(turbineX), scaler=1E-2)
    prob.driver.add_desvar('turbineY', lower=np.ones(nTurbs)*min(turbineY), upper=np.ones(nTurbs)*max(turbineY), scaler=1E-2)
    #prob.driver.add_desvar('turbineH1', 87.6, upper=None, scaler=1E-2)
    #prob.driver.add_desvar('turbineH2', 87.6, upper=None, scaler=1E-2)
    # prob.driver.add_desvar('turbineZ', lower=np.ones(nTurbs)*100., upper=np.ones(nTurbs)*500., scaler=1E-2)
    """for direction_id in range(0, windDirections.size):
        prob.driver.add_desvar('yaw%i' % direction_id, lower=-30.0, upper=30.0, scaler=1E-1)"""
    

    # add constraints
    prob.driver.add_constraint('sc', lower=np.zeros(((nTurbs-1.)*nTurbs/2.)), scaler=1.0/rotor_diameter)

    tic = time.time()
    prob.setup(check=True)
    toc = time.time()

    # print the results
    mpi_print(prob, ('FLORIS setup took %.03f sec.' % (toc-tic)))

    # time.sleep(10)
    # assign initial values to design variables
    print('Turbine X', turbineX)
    print('Turbine Y', turbineY)
    prob['turbineX'] = turbineX
    prob['turbineY'] = turbineY
    prob['turbineH1'] = turbineH1
    prob['turbineH2'] = turbineH2
    prob['nTurbsH1'] = nTurbsH1
    prob['nTurbsH2'] = nTurbsH2
    for direction_id in range(0, windDirections.size):
        prob['yaw%i' % direction_id] = yaw

    # assign values to constant inputs (not design variables)
    prob['rotorDiameter'] = rotorDiameter
    prob['axialInduction'] = axialInduction
    prob['generatorEfficiency'] = generatorEfficiency
    prob['windSpeeds'] = np.ones(nDirections)*wind_speed
    prob['air_density'] = air_density
    prob['windDirections'] = windDirections
    prob['windFrequencies'] = windFrequencies
    prob['Ct_in'] = Ct
    prob['Cp_in'] = Cp
    prob['floris_params:cos_spread'] = 1E12

    # set options
    # prob['floris_params:FLORISoriginal'] = True
    # prob['floris_params:CPcorrected'] = False
    # prob['floris_params:CTcorrected'] = False

    # run the problem
    mpi_print(prob, 'start FLORIS run')
    tic = time.time()
    # cProfile.run('prob.run()')
    prob.run()
    toc = time.time()
    #self.J = prob.check_total_derivatives()

    # print the results
    mpi_print(prob, ('FLORIS Opt. calculation took %.03f sec.' % (toc-tic)))

    """for direction_id in range(0, windDirections.size):
        mpi_print(prob,  'yaw%i (deg) = ' % direction_id, prob['yaw%i' % direction_id])"""
    # for direction_id in range(0, windDirections.size):
        # mpi_print(prob,  'velocitiesTurbines%i (m/s) = ' % direction_id, prob['velocitiesTurbines%i' % direction_id])
    # for direction_id in range(0, windDirections.size):
    #     mpi_print(prob,  'wt_power%i (kW) = ' % direction_id, prob['wt_power%i' % direction_id])

    mpi_print(prob,  'turbine X positions in wind frame (m): %s' % prob['turbineX'])
    mpi_print(prob,  'turbine Y positions in wind frame (m): %s' % prob['turbineY'])
    mpi_print(prob,  'turbine Height 1 (m): %s' % prob['turbineH1'])
    mpi_print(prob,  'turbine Height 2 (m): %s' % prob['turbineH2'])
    #mpi_print(prob,  'wind farm power in each direction (kW): %s' % prob['dirPowers'])
    #mpi_print(prob,  'AEP: %s' % prob['AEP'])
    mpi_print(prob,  'COE: %s' % prob['COE'])
    mpi_print(prob,  'AEP: %s' % prob['AEP'])

    tdCOE = prob['COE']
    tdAEP = prob['AEP']
    np.savetxt('XY5test.txt',np.c_[prob['turbineX'],prob['turbineY'],prob['turbineZ']])

    xbounds = [min(turbineX), min(turbineX), max(turbineX), max(turbineX), min(turbineX)]
    ybounds = [min(turbineY), max(turbineY), max(turbineY), min(turbineY), min(turbineX)]

    plt.figure(2)
    plt.plot(turbineX, turbineY, 'ok', label='Original')
    plt.plot(prob['turbineX'], prob['turbineY'], 'og', label='Optimized')
    plt.plot(xbounds, ybounds, ':k')
    for i in range(0, nTurbs):
        plt.plot([turbineX[i], prob['turbineX'][i]], [turbineY[i], prob['turbineY'][i]], '--k')
    plt.legend()
    plt.xlabel('Turbine X Position (m)')
    plt.ylabel('Turbine Y Position (m)')
    plt.legend(bbox_to_anchor=(1.14, 1.14))
    perDec = (tdCOE-optCOE)/tdCOE*100
    print('Percent Decrease: ', perDec)
    plt.show() 

