from florisse.COE import *
from towerse.tower import TowerSE
import numpy as np
import matplotlib.pyplot as plt
from commonse.environment import PowerWind, LogWind
import matplotlib.pyplot as plt
from openmdao.api import Group

if __name__=="__main__":

    # --- tower setup ------

    # --- geometry ----
    z_paramH1 = np.array([0.0, 43.8, 87.6])
    d_paramH1 = np.array([6.0, 4.935, 3.87])
    t_paramH1 = [0.027*1.3, 0.023*1.3, 0.019*1.3]
    n = 15
    z_fullH1 = np.linspace(0.0, 87.6, n)
    L_reinforced = 30.0*np.ones(n)  # [m] buckling length
    yawLoads = 0.0

    # --- material props ---
    E = 210e9*np.ones(n)
    G = 80.8e9*np.ones(n)
    rho = 8500.0*np.ones(n)
    sigma_y = 450.0e6*np.ones(n)

    # --- spring reaction data.  Use float('inf') for rigid constraints. ---
    kidx = np.array([0], dtype=int)  # applied at base
    kx = np.array([float('inf')])
    ky = np.array([float('inf')])
    kz = np.array([float('inf')])
    ktx = np.array([float('inf')])
    kty = np.array([float('inf')])
    ktz = np.array([float('inf')])
    nK = len(kidx)

    # --- extra mass ----
    midx = np.array([n-1], dtype=int)  # RNA mass at top
    m = np.array([285598.8])
    mIxx = np.array([1.14930678e+08])
    mIyy = np.array([2.20354030e+07])
    mIzz = np.array([1.87597425e+07])
    mIxy = np.array([0.00000000e+00])
    mIxz = np.array([5.03710467e+05])
    mIyz = np.array([0.00000000e+00])
    mrhox = np.array([-1.13197635])
    mrhoy = np.array([0.])
    mrhoz = np.array([0.50875268])
    nMass = len(midx)
    addGravityLoadForExtraMass = True
    # -----------

    # --- wind ---
    wind_zref = 90.0
    wind_z0 = 0.0
    shearExp = 0.2
    # ---------------

    # if addGravityLoadForExtraMass=True be sure not to double count by adding those force here also
    # # --- loading case 1: max Thrust ---
    wind_Uref1 = 11.73732
    plidx1 = np.array([n-1], dtype=int)  # at  top
    Fx1 = np.array([1284744.19620519])
    Fy1 = np.array([0.])
    Fz1 = np.array([-2914124.84400512])
    Mxx1 = np.array([3963732.76208099])
    Myy1 = np.array([-2275104.79420872])
    Mzz1 = np.array([-346781.68192839])
    nPL = len(plidx1)
    # # ---------------

    # # --- loading case 2: max wind speed ---
    wind_Uref2 = 70.0
    plidx2 = np.array([n-1], dtype=int)  # at  top
    Fx2 = np.array([930198.60063279])
    Fy2 = np.array([0.])
    Fz2 = np.array([-2883106.12368949])
    Mxx2 = np.array([-1683669.22411597])
    Myy2 = np.array([-2522475.34625363])
    Mzz2 = np.array([147301.97023764])
    # # ---------------

    # --- safety factors ---
    gamma_f = 1.35
    gamma_m = 1.3
    gamma_n = 1.0
    gamma_b = 1.1
    # ---------------

    # --- fatigue ---
    z_DELH1 = np.array([0.000, 1.327, 3.982, 6.636, 9.291, 11.945, 14.600, 17.255, 19.909, 22.564, 25.218, 27.873, 30.527, 33.182, 35.836, 38.491, 41.145, 43.800, 46.455, 49.109, 51.764, 54.418, 57.073, 59.727, 62.382, 65.036, 67.691, 70.345, 73.000, 75.655, 78.309, 80.964, 83.618, 86.273, 87.600])
    M_DEL = 1e3*np.array([8.2940E+003, 8.1518E+003, 7.8831E+003, 7.6099E+003, 7.3359E+003, 7.0577E+003, 6.7821E+003, 6.5119E+003, 6.2391E+003, 5.9707E+003, 5.7070E+003, 5.4500E+003, 5.2015E+003, 4.9588E+003, 4.7202E+003, 4.4884E+003, 4.2577E+003, 4.0246E+003, 3.7942E+003, 3.5664E+003, 3.3406E+003, 3.1184E+003, 2.8977E+003, 2.6811E+003, 2.4719E+003, 2.2663E+003, 2.0673E+003, 1.8769E+003, 1.7017E+003, 1.5479E+003, 1.4207E+003, 1.3304E+003, 1.2780E+003, 1.2673E+003, 1.2761E+003])
    nDEL = len(z_DELH1)
    gamma_fatigue = 1.35*1.3*1.0
    life = 20.0
    m_SN = 4
    # ---------------


    # --- constraints ---
    min_d_to_t = 120.0
    min_taper = 0.4
    # ---------------

    # # V_max = 80.0  # tip speed
    # # D = 126.0
    # # .freq1p = V_max / (D/2) / (2*pi)  # convert to Hz

    nPoints = len(z_paramH1)
    nFull = len(z_fullH1)
    wind = 'PowerWind'

    """

    z = prob['TowerSE_H1.z_full']

    print 'mass (kg) =', prob['TowerSE_H1.tower1.mass']
    print 'f1 (Hz) =', prob['TowerSE_H1.tower1.f1']
    print 'f2 (Hz) =', prob['TowerSE_H1.tower2.f2']
    print 'top_deflection1 (m) =', prob['TowerSE_H1.tower1.top_deflection']
    print 'top_deflection2 (m) =', prob['TowerSE_H1.tower2.top_deflection']
    print 'weldability =', prob['TowerSE_H1.gc.weldability']
    print 'manufacturability =', prob['TowerSE_H1.gc.manufacturability']
    print 'stress1 =', prob['TowerSE_H1.tower1.stress']
    print 'stress2 =', prob['TowerSE_H1.tower2.stress']
    print 'zs=', z
    print 'ds=', prob['TowerSE_H1.d_full']
    print 'ts=', prob['TowerSE_H1.t_full']
    print 'GL buckling =', prob['TowerSE_H1.tower1.global_buckling']
    print 'GL buckling =', prob['TowerSE_H1.tower2.global_buckling']
    print 'Shell buckling =', prob['TowerSE_H1.tower1.shell_buckling']
    print 'Shell buckling =', prob['TowerSE_H1.tower2.shell_buckling']
    print 'damage =', prob['TowerSE_H1.tower1.damage']


    print 'wind1: ', prob['TowerSE_H1.wind1.Uref']
    print 'wind2: ', prob['TowerSE_H1.wind2.Uref']


    import matplotlib.pyplot as plt
    plt.figure(figsize=(5.0, 3.5))
    plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=3)
    plt.plot(prob['TowerSE_H1.tower1.stress'], z, label='stress1')
    plt.plot(prob['TowerSE_H1.tower2.stress'], z, label='stress2')
    plt.plot(prob['TowerSE_H1.tower1.shell_buckling'], z, label='shell buckling 1')
    plt.plot(prob['TowerSE_H1.tower2.shell_buckling'], z, label='shell buckling 2')
    plt.plot(prob['TowerSE_H1.tower1.global_buckling'], z, label='global buckling 1')
    plt.plot(prob['TowerSE_H1.tower2.global_buckling'], z, label='global buckling 2')
    plt.plot(prob['TowerSE_H1.tower1.damage'], z, label='damage')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc=2)
    plt.xlabel('utilization')
    plt.ylabel('height along tower (m)')
    plt.show()


    """

    # initialize input variable arrays
    nRows = 2
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
    turbineH2 = 125
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
    nIntegrationPoints = 100
    # set up problem
    prob = Problem()
    root = prob.root = Group()


    root.add('getTurbineZ', getTurbineZ(nTurbs), promotes=['*'])
    root.add('getUeff', getUeffintegrate(nIntegrationPoints,nDirections))
    root.add('AEPGroup', AEPGroup(nTurbs, nDirections=nDirections,
                use_rotor_components=False, datasize=0, differentiable=True,
                optimizingLayout=False, nSamples=0), promotes=['*'])
    root.add('COEComponent', COEComponent(nTurbs), promotes=['*'])
    root.add('TowerSE_H1', TowerSE(nPoints, nFull, nK, nMass, nPL, nDEL, wind=wind))

    root.connect('getUeff.H1','turbineH1')
    root.connect('getUeff.H2','turbineH2')
    root.connect('getUeff.windSpeeds', 'windSpeeds')

    #root.connect('windSpeeds', )
    # initialize problem
    prob.setup()

    prob['turbineH1'] = turbineH1
    prob['turbineH2'] = turbineH2
    prob['H1_H2'] = H1_H2

    prob['turbineX'] = turbineX
    prob['turbineY'] = turbineY
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
    prob['getUeff.rotorDiameter'] = 126.4
    prob['getUeff.nPoints'] = nIntegrationPoints
    prob['getUeff.zref'] = 90
    prob['getUeff.z0'] = 0


    """**********************************************************************"""
    #TowerSE
    if wind=='PowerWind':
        prob['TowerSE_H1.wind1.shearExp'] = shearExp
        prob['TowerSE_H1.wind2.shearExp'] = shearExp

    # assign values to params

    # --- geometry ----
    prob['TowerSE_H1.z_param'] = z_paramH1
    prob['TowerSE_H1.d_param'] = d_paramH1
    prob['TowerSE_H1.t_param'] = t_paramH1
    prob['TowerSE_H1.z_full'] = z_fullH1
    prob['TowerSE_H1.tower1.L_reinforced'] = L_reinforced
    prob['TowerSE_H1.distLoads1.yaw'] = yawLoads

    # --- material props ---
    prob['TowerSE_H1.tower1.E'] = E
    prob['TowerSE_H1.tower1.G'] = G
    prob['TowerSE_H1.tower1.rho'] = rho
    prob['TowerSE_H1.tower1.sigma_y'] = sigma_y

    # --- spring reaction data.  Use float('inf') for rigid constraints. ---
    prob['TowerSE_H1.tower1.kidx'] = kidx
    prob['TowerSE_H1.tower1.kx'] = kx
    prob['TowerSE_H1.tower1.ky'] = ky
    prob['TowerSE_H1.tower1.kz'] = kz
    prob['TowerSE_H1.tower1.ktx'] = ktx
    prob['TowerSE_H1.tower1.kty'] = kty
    prob['TowerSE_H1.tower1.ktz'] = ktz

    # --- extra mass ----
    prob['TowerSE_H1.tower1.midx'] = midx
    prob['TowerSE_H1.tower1.m'] = m
    prob['TowerSE_H1.tower1.mIxx'] = mIxx
    prob['TowerSE_H1.tower1.mIyy'] = mIyy
    prob['TowerSE_H1.tower1.mIzz'] = mIzz
    prob['TowerSE_H1.tower1.mIxy'] = mIxy
    prob['TowerSE_H1.tower1.mIxz'] = mIxz
    prob['TowerSE_H1.tower1.mIyz'] = mIyz
    prob['TowerSE_H1.tower1.mrhox'] = mrhox
    prob['TowerSE_H1.tower1.mrhoy'] = mrhoy
    prob['TowerSE_H1.tower1.mrhoz'] = mrhoz
    prob['TowerSE_H1.tower1.addGravityLoadForExtraMass'] = addGravityLoadForExtraMass
    # -----------

    # --- wind ---
    prob['TowerSE_H1.wind1.zref'] = wind_zref
    prob['TowerSE_H1.wind1.z0'] = wind_z0
    # ---------------

    # # --- loading case 1: max Thrust ---
    prob['TowerSE_H1.wind1.Uref'] = wind_Uref1
    prob['TowerSE_H1.tower1.plidx'] = plidx1
    prob['TowerSE_H1.tower1.Fx'] = Fx1
    prob['TowerSE_H1.tower1.Fy'] = Fy1
    prob['TowerSE_H1.tower1.Fz'] = Fz1
    prob['TowerSE_H1.tower1.Mxx'] = Mxx1
    prob['TowerSE_H1.tower1.Myy'] = Myy1
    prob['TowerSE_H1.tower1.Mzz'] = Mzz1
    # # ---------------

    # # --- loading case 2: max Wind Speed ---
    prob['TowerSE_H1.wind2.Uref'] = wind_Uref2
    prob['TowerSE_H1.tower2.plidx'] = plidx2
    prob['TowerSE_H1.tower2.Fx'] = Fx2
    prob['TowerSE_H1.tower2.Fy'] = Fy2
    prob['TowerSE_H1.tower2.Fz'] = Fz2
    prob['TowerSE_H1.tower2.Mxx'] = Mxx2
    prob['TowerSE_H1.tower2.Myy'] = Myy2
    prob['TowerSE_H1.tower2.Mzz'] = Mzz2
    # # ---------------

    # --- safety factors ---
    prob['TowerSE_H1.tower1.gamma_f'] = gamma_f
    prob['TowerSE_H1.tower1.gamma_m'] = gamma_m
    prob['TowerSE_H1.tower1.gamma_n'] = gamma_n
    prob['TowerSE_H1.tower1.gamma_b'] = gamma_b
    # ---------------

    # --- fatigue ---
    prob['TowerSE_H1.tower1.z_DEL'] = z_DELH1
    prob['TowerSE_H1.tower1.M_DEL'] = M_DEL
    prob['TowerSE_H1.tower1.gamma_fatigue'] = gamma_fatigue
    prob['TowerSE_H1.tower1.life'] = life
    prob['TowerSE_H1.tower1.m_SN'] = m_SN
    # ---------------

    # --- constraints ---
    prob['TowerSE_H1.gc.min_d_to_t'] = min_d_to_t
    prob['TowerSE_H1.gc.min_taper'] = min_taper
    # ---------------
    """**********************************************************************"""


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

    x = np.linspace(0,nDirections,nDirections)
    plt.figure()
    plt.plot(x,windSpeeds, label='Measured wind speed at %s m'%prob['H1PowerWind_0.zref'])
    plt.plot(x,prob['getUeff.UeffH1'], label='H1: %s m'%turbineH1)
    plt.plot(x,prob['getUeff.UeffH2'], label='H2: %s m'%turbineH2)
    plt.legend(loc=4)
    plt.xlabel('Direction #')
    plt.ylabel('Effective Wind Speed (m/s)')
    plt.show()
