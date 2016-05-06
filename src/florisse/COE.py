import numpy as np
from math import pi
import time
from openmdao.api import Group, Component, Problem, ScipyGMRES
from florisse.floris import AEPGroup


class COEComponent(Component):
    """
    Componenet to calculate the cost of energy (COE)
    """

    def __init__(self, nTurbines):

        super(COEComponent, self).__init__()

        self.fd_options['form'] = 'forward'
        self.fd_options['step_size'] = 1.0e-6
        self.fd_options['step_type'] = 'relative'
        self.fd_options['force_fd'] = True
        

        self.add_param('turbineZ', np.zeros(nTurbines), units='m',
                       desc='z coordinates of turbines')
        self.add_param('AEP', 0.0, desc='AEP of the wind farm')
        self.add_param('rotorDiameter', np.zeros(nTurbines) + 126.4, units='m', desc='rotor diameter of each turbine')

        self.add_output('COE', 0.0, desc='Cost of Energy for the wind farm')

    def solve_nonlinear(self, params, unknowns, resids):

        turbineZ = params['turbineZ']
        AEP = params['AEP']
        RotorDiam = params['rotorDiameter']

        nTurbines = len(turbineZ)

        # Local Variables
        fixed_charge_rate = 0.095
        tax_rate = 0.4
        ppi_mat   = 1.0465528035
        slope   = 13.0
        intercept     = 5813.9
        bos = 559. * 5e3
        array_losses = 0.059
        other_losses = 0.0
        availability = 0.94
        losses = availability * (1-array_losses) * (1-other_losses)
        assemblyCostMultiplier = 0.30
        profitMultiplier = 0.20
        overheadCostMultiplier = 0.0
        transportMultiplier = 0.0

        rotor_cost = 1505102.53 
        nacelle_cost = 3000270.

        windpactMassSlope = 0.397251147546925
        windpactMassInt   = -1414.381881
        
        twrCostEscalator  = 1.5944
        twrCostCoeff      = 1.5 # $/kg   

        tower_cost = np.zeros(nTurbines)
        for i in range(nTurbines):
            mass = windpactMassSlope * pi * (RotorDiam[i]/2.)**2 * turbineZ[i] + windpactMassInt
            tower_cost[i] = mass*twrCostEscalator*twrCostCoeff
            # tower_cost = 1390588.80 # to change

        parts_cost_farm = nTurbines*(rotor_cost + nacelle_cost) + np.sum(tower_cost) #parts cost for the entire wind farm
        turbine_multiplier = (1 + transportMultiplier + profitMultiplier) * (1+overheadCostMultiplier+assemblyCostMultiplier)
        turbine_cost = turbine_multiplier * parts_cost_farm

        unknowns['COE'] = fixed_charge_rate*(turbine_cost+bos)/AEP + 0.0122*(1-tax_rate)





class getTurbineZ(Component):

    def __init__(self, nTurbines):

        super(getTurbineZ, self).__init__()

        self.add_param('turbineH1', 0.0, units='m', desc='Turbine height 1')
        self.add_param('turbineH2', 0.0, units='m', desc='Turbine height 2')
        self.add_param('nTurbsH1', 1, desc='The number of turbines of height 1')
        self.add_param('nTurbsH2', 1, desc='The number of turbines of height 2')

        self.add_output('turbineZ', np.zeros(nTurbines), units='m', desc='The array of turbine heights')


    def solve_nonlinear(self, params, unknowns, resids):
        turbineH1 = params['turbineH1']
        turbineH2 = params['turbineH2']
        nTurbsH1 = params['nTurbsH1']
        nTurbsH2 = params['nTurbsH2']

        unknowns['turbineZ'] = np.hstack([np.ones(nTurbsH1)*turbineH1, np.ones(nTurbsH2)*turbineH2])

    def linearize(self, params, unknowns, resids):
        turbineH1 = params['turbineH1']
        turbineH2 = params['turbineH2']
        nTurbsH1 = params['nTurbsH1']
        nTurbsH2 = params['nTurbsH2']

        J = {}

        J['turbineZ', 'turbineH1'] = np.hstack([np.ones(nTurbsH1), np.zeros(nTurbsH2)])
        J['turbineZ', 'turbineH2'] = np.hstack([np.zeros(nTurbsH1), np.ones(nTurbsH2)])
        return J




if __name__=="__main__":
    """
    This is just to test during development
    """
    rotor_diameter = 126.4
    nRows = 3
    spacing = 5     # turbine grid spacing in diameters

    # Set up position arrays
    points = np.linspace(start=spacing*rotor_diameter, stop=nRows*spacing*rotor_diameter, num=nRows)
    xpoints, ypoints = np.meshgrid(points, points)
    turbineX = np.ndarray.flatten(xpoints)
    turbineY = np.ndarray.flatten(ypoints)
    

    # initialize input variable arrays
    nTurbs = turbineX.size
    #turbineZ = np.array([100,100,100,100,100,200,200,200,200])
    rotorDiameter = np.zeros(nTurbs)
    axialInduction = np.zeros(nTurbs)
    Ct = np.zeros(nTurbs)
    Cp = np.zeros(nTurbs)
    generatorEfficiency = np.zeros(nTurbs)
    yaw = np.zeros(nTurbs)

    # define initial values
    for turbI in range(0, nTurbs):
        rotorDiameter[turbI] = 126.4            # m
        axialInduction[turbI] = 1.0/3.0
        Ct[turbI] = 4.0*axialInduction[turbI]*(1.0-axialInduction[turbI])
        # Cp[turbI] = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
        Cp[turbI] = 0.7737 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
        generatorEfficiency[turbI] = 1.0#0.944
        yaw[turbI] = 0.     # deg.

    # Define flow properties
    wind_speed = 8.0        # m/s
    air_density = 1.1716    # kg/m^3
    # wind_direction = 240    # deg (N = 0 deg., using direction FROM, as in met-mast data)
    
    turbineH1 = 90.
    turbineH2 = 100.
    nTurbsH1 = 3
    nTurbsH2 = nTurbs-nTurbsH1
    rotorDiameter = np.ones(nTurbs)*126.4 
    nDirections = 50
    wind_frequency = 1./nDirections    # probability of wind in this direction at this speed

    # set up problem
    prob = Problem()
    root = prob.root = Group()
    

    root.add('getTurbineZ', getTurbineZ(nTurbs), promotes=['*'])
    root.add('COEComponent', COEComponent(nTurbs), promotes=['*'])
    root.add('AEPGroup', AEPGroup(nTurbs, nDirections=nDirections,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['*'])
    #root.ln_solver = ScipyGMRES()

    # initialize problem
    prob.setup()

    prob['turbineH1'] = turbineH1
    prob['turbineH2'] = turbineH2
    prob['nTurbsH1'] = nTurbsH1
    prob['nTurbsH2'] = nTurbsH2

    prob['turbineX'] = turbineX
    prob['turbineY'] = turbineY
    # prob['turbineZ'] = turbineZ
    prob['yaw0'] = yaw

    # assign values to constant inputs (not design variables)
    prob['rotorDiameter'] = rotorDiameter
    prob['axialInduction'] = axialInduction
    prob['generatorEfficiency'] = generatorEfficiency
    prob['windSpeeds'] = np.array([wind_speed])
    prob['air_density'] = air_density
    #prob['windDirections'] = np.array([wind_direction])
    prob['windFrequencies'] = np.ones([nDirections])*wind_frequency
    prob['Ct_in'] = Ct
    prob['Cp_in'] = Cp
    prob['floris_params:cos_spread'] = 1E12         # turns off cosine spread (just needs to be very large)

    # run the problem
    print 'start run'
    tic = time.time()
    prob.run()
    toc = time.time()

    print 'turbineZ: ', prob['turbineZ']
    print 'AEP: ', prob['AEP']
    print 'COE: ', prob['COE']
    
