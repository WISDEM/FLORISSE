import numpy as np
from math import pi, log
import time
from datetime import datetime
from openmdao.api import Group, Component, Problem, ScipyGMRES
from florisse.floris import AEPGroup
from commonse.environment import PowerWind, LogWind

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

        #windpactMassSlope = 0.397251147546925
        #windpactMassInt   = -1414.381881

        tower_mass_coeff = 19.828
        tower_mass_exp = 2.0282


        #twrCostEscalator  = 1.5944
        #twrCostCoeff      = 1.5 # $/kg

        tower_mass_cost_coefficient = 3.08 #$/kg

        tower_cost = np.zeros(nTurbines)
        for i in range(nTurbines):
            #mass = windpactMassSlope * pi * (RotorDiam[i]/2.)**2 * turbineZ[i] + windpactMassInt
            mass = tower_mass_coeff*turbineZ[i]**tower_mass_exp #new mass from Katherine
            #tower_cost[i] = mass*twrCostEscalator*twrCostCoeff
            # tower_cost = 1390588.80 # to change
            tower_cost[i] = tower_mass_cost_coefficient*mass #new cost from Katherine


        parts_cost_farm = nTurbines*(rotor_cost + nacelle_cost) + np.sum(tower_cost) #parts cost for the entire wind farm
        turbine_multiplier = (1 + transportMultiplier + profitMultiplier) * (1+overheadCostMultiplier+assemblyCostMultiplier)
        turbine_cost = turbine_multiplier * parts_cost_farm

        unknowns['COE'] = (fixed_charge_rate*(turbine_cost+bos)/AEP + 0.0122*(1-tax_rate))*1E8


class getTurbineZ(Component):

    def __init__(self, nTurbines):

        super(getTurbineZ, self).__init__()

        self.add_param('turbineH1', 0.0, units='m', desc='Turbine height 1')
        self.add_param('turbineH2', 0.0, units='m', desc='Turbine height 2')
        #self.add_param('nTurbsH1', 1, desc='The number of turbines of height 1')
        #self.add_param('nTurbsH2', 1, desc='The number of turbines of height 2')
        self.add_param('H1_H2', np.zeros(nTurbines), desc='An array indicating which turbines are of each height: 0 indicates H1, 1 indicates H2')

        self.add_output('turbineZ', np.zeros(nTurbines), units='m', desc='The array of turbine heights')


    def solve_nonlinear(self, params, unknowns, resids):
        turbineH1 = params['turbineH1']
        turbineH2 = params['turbineH2']
        H1_H2 = params['H1_H2']
        nTurbines = len(H1_H2)
        #nTurbsH1 = params['nTurbsH1']
        #nTurbsH2 = params['nTurbsH2']

        turbineZ = np.array([])
        for i in range(nTurbines):
            if H1_H2[i] == 0:
                turbineZ = np.append(turbineZ, turbineH1)
            elif H1_H2[i] == 1:
                turbineZ = np.append(turbineZ, turbineH2)
        #unknowns['turbineZ'] = np.hstack([np.ones(nTurbsH1)*turbineH1, np.ones(nTurbsH2)*turbineH2])
        unknowns['turbineZ'] = turbineZ


    def linearize(self, params, unknowns, resids):
        turbineH1 = params['turbineH1']
        turbineH2 = params['turbineH2']
        H1_H2 = params['H1_H2']
        nTurbs = len(H1_H2)


        J = {}

        J['turbineZ', 'turbineH1'] = np.array([])
        for i in range(nTurbs):
            if H1_H2[i] == 0:
                J['turbineZ', 'turbineH1'] = np.append(J['turbineZ', 'turbineH1'], 1)
            else:
                J['turbineZ', 'turbineH1'] = np.append(J['turbineZ', 'turbineH1'], 0)

        J['turbineZ', 'turbineH2'] = np.array([])
        for i in range(nTurbs):
            if H1_H2[i] == 0:
                J['turbineZ', 'turbineH2'] = np.append(J['turbineZ', 'turbineH2'], 0)
            else:
                J['turbineZ', 'turbineH2'] = np.append(J['turbineZ', 'turbineH2'], 1)

        return J


class AEPobj(Component):
    """
    Objective to maximize AEP
    """

    def __init__(self):

        super(AEPobj, self).__init__()

        self.fd_options['form'] = 'forward'
        self.fd_options['step_size'] = 1.0e-6
        self.fd_options['step_type'] = 'relative'
        self.fd_options['force_fd'] = True

        self.add_param('AEP', 0.0, desc='AEP of the wind farm')

        self.add_output('maxAEP', 0.0, desc='negative AEP')

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['maxAEP'] = -1*params['AEP']


class getUeffintegrateU1_U2(Component):
    """
    Integrate across the turbine to get effective wind speed
    """
    def __init__(self, nDirections):

        super(getUeffintegrate, self).__init__()

        self.fd_options['form'] = 'forward'
        self.fd_options['step_size'] = 1.0e-6
        self.fd_options['step_type'] = 'relative'
        self.fd_options['force_fd'] = True

        # inputs
        self.add_param('nIntegrationPoints', 5, desc='number of integration points')
        self.add_param('rotorDiameter', 126.4, units='m', desc='rotor diameter of each turbine')
        self.add_param('turbineH1', 87.6, units='m', desc='height of turbine1')
        self.add_param('turbineH2', 87.6, units='m', desc='height of turbine2')
        self.add_param('wind', 'PowerWind', desc='Wind shear calculation method')
        self.add_param('Uref', np.zeros(nDirections), units='m/s', desc='refenence wind speed for each direction')
        self.add_param('zref', 90, units='m', desc='height at which Uref was measured')
        self.add_param('z_roughness', 0.01, units='m', desc='ground roughness height')
        self.add_param('z0', 0, units='m', desc='height of ground')
        self.add_param('shearExp', 0.2, desc='PowerWind exponent')

        # outputs
        self.add_output('UeffH1', np.zeros(nDirections), units='m/s', desc='The effective wind speed on turbine with height H1 from each direction')
        self.add_output('UeffH2', np.zeros(nDirections), units='m/s', desc='The effective wind speed on turbine with height H2 from each direction')



    def solve_nonlinear(self, params, unknowns, resids):

        D = params['rotorDiameter']
        r = D/2.

        nPoints = params['nIntegrationPoints']
        wind = params['wind']
        turbineH1 = params['turbineH1']
        turbineH2 = params['turbineH2']
        Uref = params['Uref']
        nDirections = len(Uref)
        zref = params['zref']
        z_roughness = params['z_roughness']
        z0 = params['z0']
        shearExp = params['shearExp']
        UeffH1 = np.zeros(nDirections)
        UeffH2 = np.zeros(nDirections)

        #start at bottom work the way up to the top
        for j in range(nDirections):
            z1 = turbineH1-r
            z2 = turbineH2-r
            Usum1 = 0.
            Usum2 = 0.
            Asum = 0
            for i in range(nPoints+1):
                dz = D/nPoints
                if i == 0 or i == nPoints:
                    dz = dz/2.
                if z1 < turbineH1:
                    a1 = 2*np.sqrt(r**2-(turbineH1-z1)**2)
                elif z1 == turbineH1:
                    a1 = D
                else:
                    a1 = 2*np.sqrt(r**2-(z1-turbineH1)**2)
                if z1+dz < turbineH1:
                    a2 = 2*np.sqrt(r**2-(turbineH1-(z1+dz))**2)
                elif z1+dz == turbineH1:
                    a2 = D
                else:
                    a2 = 2*np.sqrt(r**2-(z1+dz-turbineH1)**2)
                if wind == 'PowerWind':
                    U1b = PowWind(Uref[j], z1, zref, z0, shearExp)
                    U1t = PowWind(Uref[j], z1+dz, zref, z0, shearExp)
                    U2b = PowWind(Uref[j], z2, zref, z0, shearExp)
                    U2t = PowWind(Uref[j], z2+dz, zref, z0, shearExp)
                if wind == 'LogWind':
                    U1b = LnWind(Uref[j], z, z0, z_roughness, zref)
                    U1t = LnWind(Uref[j], z, z0, z_roughness, zref)
                    U2b = LnWind(Uref[j], z, z0, z_roughness, zref)
                    U2t = LnWind(Uref[j], z, z0, z_roughness, zref)

                Usum1 += dz/2.*(a1*U1b+a2*U1t)
                Usum2 += dz/2.*(a1*U2b+a2*U2t)
                Asum += dz/2.*(a1+a2)
                z1 += dz
                z2 += dz
            UeffH1[j] = Usum1/Asum
            UeffH2[j] = Usum2/Asum
        unknowns['UeffH1'] = UeffH1
        unknowns['UeffH2'] = UeffH2


if __name__=="__main__":
    """
    This is just to test during development
    """
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

    turbineH1 = 100
    turbineH2 = 1000
    nTurbsH1 = nTurbs/2
    nTurbsH2 = nTurbs-nTurbsH1
    H1_H2 = np.array([0,1,0,1,0,1,0,1,0], dtype=int)
    print H1_H2
    rotorDiameter = np.ones(nTurbs)*126.4
    nDirections = 50
    wind_frequency = 1./nDirections    # probability of wind in this direction at this speed

    # set up problem
    prob = Problem()
    root = prob.root = Group()


    root.add('getTurbineZ', getTurbineZ(nTurbs), promotes=['*'])
    root.add('AEPGroup', AEPGroup(nTurbs, nDirections=nDirections,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['*'])
    root.add('COEComponent', COEComponent(nTurbs), promotes=['*'])
    root.add('getUeff', getUeffintegrate(10,11))

    #root.ln_solver = ScipyGMRES()

    # initialize problem
    prob.setup()

    prob['turbineH1'] = turbineH1
    prob['turbineH2'] = turbineH2
    #prob['nTurbsH1'] = nTurbsH1
    #prob['nTurbsH2'] = nTurbsH2
    prob['H1_H2'] = H1_H2

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

    COE5 = prob['COE']
    AEP5 = prob['AEP']

    print 'turbineZ5: ', prob['turbineZ']
    print 'AEP5: ', prob['AEP']
    print 'COE5: ', prob['COE']
    """
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
    nDirections = len(windSpeeds)
    nIntegrationPoints = 100
    turbineH1 = 75
    turbineH2 = 90
    nTurbines = 15
    turbineZ = np.ones(nTurbines)*75.
    rotorDiameter = np.ones(nTurbines)*126.4
    wind = 'PowerWind'

    prob = Problem()
    root = prob.root = Group()

    root.add('p', getUeffintegrate(nDirections, nTurbines))
    prob.setup()

    prob['p.rotorDiameter'] = rotorDiameter
    prob['p.turbineZ'] = turbineZ
    prob['p.nIntegrationPoints'] = 100
    prob['p.Uref'] = windSpeeds
    prob['p.wind'] = wind

    start = time.time()
    prob.run()
    print 'Time to run: ', time.time() - start

    #print 'Uref: ', windSpeeds
    #print 'WindSpeeds: ', prob['p.Ueff']
    #print np.shape(prob['p.Ueff'])


    """
    ref1 = prob['p.UeffH1']
    ref2 = prob['p.UeffH2']

    error1 = np.zeros(100)
    error2 = np.zeros(100)

    import matplotlib.pyplot as plt

    x = np.linspace(1,nDirections,nDirections)

    plt.figure(1)
    plt.plot(x, windSpeeds, label='Reference')
    plt.plot(x, prob['p.UeffH1'], label='Ueff1')
    plt.plot(x, prob['p.UeffH2'], label='Ueff2')

    plt.legend()

    for i in range(0, 100):
        prob = Problem()
        root = prob.root = Group()

        root.add('p', getUeffintegrate(nDirections))

        prob.setup()

        prob['p.turbineH1'] = turbineH1
        prob['p.turbineH2'] = turbineH2
        prob['p.nIntegrationPoints'] = i+1
        prob['p.Uref'] = windSpeeds
        prob['p.wind'] = wind

        prob.run()

        print i

        error1[i] = (np.sum(ref1)-np.sum(prob['p.UeffH1']))/np.sum(ref1)*100
        error2[i] = (np.sum(ref2)-np.sum(prob['p.UeffH2']))/np.sum(ref2)*100

    plt.figure(2)
    y = np.linspace(1,100,100)
    plt.plot(y, error1, label='error1')
    plt.plot(y, error2, label='error2')
    plt.legend(loc=4)
    plt.show()

    times = np.zeros(100)

    #Testing time
    for i in range(100):
        sumtime = 0
        for j in range(25):
            start = time.time()
            prob = Problem()
            root = prob.root = Group()

            root.add('p', getUeffintegrate(nDirections))

            prob.setup()

            prob['p.turbineH1'] = turbineH1
            prob['p.turbineH2'] = turbineH2
            prob['p.nIntegrationPoints'] = i+1
            prob['p.Uref'] = windSpeeds
            prob['p.wind'] = wind

            prob.run()

            sumtime += time.time()-start
            print i
        times[i] = sumtime/25

    plt.figure(3)
    plt.plot(y,times)
    plt.show()
    """
