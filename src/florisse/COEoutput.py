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

    size = 36 # number of processors (and number of wind directions to run)

    #########################################################################
    # define turbine size
    rotor_diameter = 126.4  # (m)
    
    
    numRows = 6

    filename = "XYZ%stest.txt"%(numRows)

    file = open(filename)
    xin = np.loadtxt(file)
    n = len(xin)
    turbineX = np.zeros(n)
    turbineY = np.zeros(n)
    turbineZ = np.zeros(n)
    for i in range(n):
        turbineX[i] = xin[i][0]
        turbineY[i] = xin[i][1]
        turbineZ[i] = xin[i][2]


    print('Turbine X', turbineX)
    print('Turbine Y', turbineY)
    # initialize input variable arrays
    nTurbs = turbineX.size
    #87.6m is the height of the 5MW NREL turbines
    turbineMin = rotor_diameter/2+10
    turbineH1 = 180
    turbineH2 = 73.2
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
    prob = Problem()
    root = prob.root = Group()
    
    root.add('AEPGroup', AEPGroup(nTurbs, nDirections=size,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['*'])
    root.add('COEComponent', COEComponent(nTurbs), promotes=['*'])
    
    tic = time.time()
    prob.setup(check=True)
    toc = time.time()

    # print the results
    mpi_print(prob, ('FLORIS setup took %.03f sec.' % (toc-tic)))

    # assign initial values to design variables
    prob['turbineX'] = turbineX
    prob['turbineY'] = turbineY
    prob['turbineZ'] = turbineZ
    #prob['turbineH1'] = turbineH1
    #prob['turbineH2'] = turbineH2
    #prob['nTurbsH1'] = nTurbsH1
    #prob['nTurbsH2'] = nTurbsH2
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
    
    mpi_print(prob,  'turbine X positions in wind frame (m): %s' % prob['turbineX'])
    mpi_print(prob,  'turbine Y positions in wind frame (m): %s' % prob['turbineY'])
    #mpi_print(prob,  'wind farm power in each direction (kW): %s' % prob['dirPowers'])
    #mpi_print(prob,  'AEP: %s' % prob['AEP'])
    mpi_print(prob,  '3D COE: %s' % prob['COE'])
    mpi_print(prob,  '3D AEP: %s' % prob['AEP'])

    optCOE = prob['COE']
    zTurbs = prob['turbineZ']







    filename = "XY%stest.txt"%(numRows)

    file = open(filename)
    xin = np.loadtxt(file)
    n = len(xin)
    turbineX = np.zeros(n)
    turbineY = np.zeros(n)
    turbineZ = np.zeros(n)
    for i in range(n):
        turbineX[i] = xin[i][0]
        turbineY[i] = xin[i][1]
        turbineZ[i] = xin[i][2]


    prob = Problem()
    root = prob.root = Group()
    
    root.add('AEPGroup', AEPGroup(nTurbs, nDirections=size,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['*'])
    root.add('COEComponent', COEComponent(nTurbs), promotes=['*'])
    
    tic = time.time()
    prob.setup(check=True)
    toc = time.time()

    # print the results
    mpi_print(prob, ('FLORIS setup took %.03f sec.' % (toc-tic)))

    # assign initial values to design variables
    prob['turbineX'] = turbineX
    prob['turbineY'] = turbineY
    prob['turbineZ'] = turbineZ
    #prob['turbineH1'] = turbineH1
    #prob['turbineH2'] = turbineH2
    #prob['nTurbsH1'] = nTurbsH1
    #prob['nTurbsH2'] = nTurbsH2
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
    
    mpi_print(prob,  'turbine X positions in wind frame (m): %s' % prob['turbineX'])
    mpi_print(prob,  'turbine Y positions in wind frame (m): %s' % prob['turbineY'])
    #mpi_print(prob,  'wind farm power in each direction (kW): %s' % prob['dirPowers'])
    #mpi_print(prob,  'AEP: %s' % prob['AEP'])
    mpi_print(prob,  'COE: %s' % prob['COE'])
    mpi_print(prob,  'AEP: %s' % prob['AEP'])
    
    tdCOE = prob['COE']

    perDec = (tdCOE-optCOE)/tdCOE*100
    print('Percent Decrease: ', perDec)
    print('zTurbine: ', zTurbs)
