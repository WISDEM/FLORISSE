from florisse.COE import *
import numpy as np
import matplotlib.pyplot as plt
from commonse.environment import PowerWind, LogWind

if __name__=="__main__":

    # initialize input variable arrays
    nRows = 3
    nTurbs = nRows**2
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
    # wind_speed = 8.0        # m/s
    air_density = 1.1716    # kg/m^3
    # wind_direction = 240    # deg (N = 0 deg., using direction FROM, as in met-mast data)

    turbineH1 = 75
    turbineH2 = 150
    nTurbsH1 = nTurbs/2
    nTurbsH2 = nTurbs-nTurbsH1
    H1_H2 = np.array([])
    for i in range(nTurbs/2):
        H1_H2 = np.append(H1_H2, 0)
        H1_H2 = np.append(H1_H2, 1)
    if len(H1_H2) < nTurbs:
        H1_H2 = np.append(H1_H2, 0)
    print H1_H2
    rotorDiameter = np.ones(nTurbs)*126.4
    # nDirections = 50
    # wind_frequency = 1./nDirections    # probability of wind in this direction at this speed

    rotor_diameter = 126.4

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

    nDirections = len(windSpeeds)
    windDirections = np.linspace(0,360-360/nDirections, nDirections)
    print(len(windSpeeds))
    print(len(windFrequencies))



    spacing = 5     # turbine grid spacing in diameters

    # Set up position arrays
    points = np.linspace(start=spacing*rotor_diameter, stop=nRows*spacing*rotor_diameter, num=nRows)
    xpoints, ypoints = np.meshgrid(points, points)
    turbineX = np.ndarray.flatten(xpoints)
    turbineY = np.ndarray.flatten(ypoints)
    nPoints = 100

    # set up problem
    prob = Problem()
    root = prob.root = Group()


    root.add('getTurbineZ', getTurbineZ(nTurbs), promotes=['*'])
    for i in range(nDirections):
        root.add('H1PowerWind_%s'%i, PowerWind(nPoints))
    for i in range(nDirections):
        root.add('H2PowerWind_%s'%i, PowerWind(nPoints))
    root.add('getUeff', getUeffintegrate(nPoints,nDirections))
    root.add('AEPGroup', AEPGroup(nTurbs, nDirections=nDirections,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['*'])
    root.add('COEComponent', COEComponent(nTurbs), promotes=['*'])

    for i in range(nDirections):
        root.connect('H1PowerWind_%s.U'%i, 'getUeff.windSpeedsH1_%s'%i)
        root.connect('H2PowerWind_%s.U'%i, 'getUeff.windSpeedsH2_%s'%i)

    root.connect('windSpeeds', )
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
    prob['windSpeeds'] = windSpeeds
    prob['air_density'] = air_density
    prob['windDirections'] = windDirections
    prob['windFrequencies'] = windFrequencies
    prob['Ct_in'] = Ct
    prob['Cp_in'] = Cp
    prob['floris_params:cos_spread'] = 1E12         # turns off cosine spread (just needs to be very large)

    # Wind Data
    for i in range(nDirections):
        prob['H1PowerWind_%s.Uref'%i] = windSpeeds[i]
        prob['H1PowerWind_%s.zref'%i] = 87.6
        prob['H1PowerWind_%s.z'%i] = np.linspace(prob['turbineH1']-rotorDiameter[0]/2, prob['turbineH1']+rotorDiameter[0]/2, nPoints)

        prob['H2PowerWind_%s.Uref'%i] = windSpeeds[i]
        prob['H2PowerWind_%s.zref'%i] = 87.6
        prob['H2PowerWind_%s.z'%i] = np.linspace(prob['turbineH2']-rotorDiameter[0]/2, prob['turbineH2']+rotorDiameter[0]/2, nPoints)

    prob['getUeff.rotorDiameter'] = 126.4


    # run the problem
    print 'start run'
    tic = time.time()
    prob.run()
    toc = time.time()

    print 'turbineZ: ', prob['turbineZ']
    print 'AEP: ', prob['AEP']
    print 'COE: ', prob['COE']
    print 'Uref: ', windSpeeds
    print 'UeffH1: ', prob['getUeff.UeffH1']
    print 'UeffH2: ', prob['getUeff.UeffH2']
