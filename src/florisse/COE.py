import numpy as np
from math import pi
import time
from openmdao.api import Group, Component, Problem
from florisse.floris import AEPGroup


class COE(Component):
    """
    Componenet to calculate the cost of energy (COE)
    """

    def __init__(self, nTurbines):

        super(COE, self).__init__()

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

    def __init__(self):

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

    # define turbine locations in global reference frame
    turbineH1 = 90.
    turbineH2 = 150.
    nTurbsH1 = 3
    nTurbsH2 = 5
    nTurbines = nTurbsH1+nTurbsH2
    rotorDiameter = np.ones(nTurbines)*126.4 
    nDirections = 72

    # set up problem
    prob = Problem()
    root = prob.root = Group()
    
    #TODO How to I pass TurbineZ down to calculate AEP? Basically, How do I get AEP?

    root.add('getTurbineZ', getTurbineZ(), promotes=['*'])
    root.add('COE', COE(nTurbines), promotes=['*'])
    root.add('AEP', AEPGroup(nTurbines, nDirections=nDirections,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['AEP'])

    # initialize problem
    prob.setup()

    prob['turbineH1'] = turbineH1
    prob['turbineH2'] = turbineH2
    prob['nTurbsH1'] = nTurbsH1
    prob['nTurbsH2'] = nTurbsH2
    prob['rotorDiameter'] = rotorDiameter

    # run the problem
    print 'start run'
    tic = time.time()
    prob.run()
    toc = time.time()

    print 'turbineZ: ', prob['turbineZ']
    print 'COE: ', prob['COE']
    
