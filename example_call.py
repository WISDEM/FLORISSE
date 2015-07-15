# from openmdao.main.api import Assembly
from example_assembly import floris_assembly_opt
from FLORIS_visualization import floris_assembly
# from FLORIS_components import floris_assembly as vis_assembly
import time
import numpy as np


if __name__ == "__main__":

    myFloris = floris_assembly_opt()

    x = 0.
    y = 250
    # define all turbine properties
    # turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1, 1164.7+x, 947.2+x,  1682.4+x, 1464.9+x, 1982.6+x, 2200.1+x])
    # turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7, 1024.7+y, 1335.3+y, 1387.2+y, 1697.8+y, 2060.3+y, 1749.7+y])
    turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1])#*0.05
    turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7])#*0.05

    nTurbs = turbineX.size
    position = np.zeros([nTurbs, 2])
    rotorDiameter = np.zeros(nTurbs)
    rotorArea = np.zeros(nTurbs)
    axialInduction = np.zeros(nTurbs)
    Ct = np.zeros(nTurbs)
    Cp = np.zeros(nTurbs)
    generator_efficiency = np.zeros(nTurbs)
    yaw = np.zeros(nTurbs)

    resolution = 200
    xlim = [0, 3000]
    ylim = [0, 3000]
    x = np.linspace(xlim[0], xlim[1], resolution)
    y = np.linspace(ylim[0], ylim[1], resolution)
    x, y = np.meshgrid(x, y)
    myFloris.ws_position = np.array([x.flatten(), y.flatten()]).transpose()

    for turbI in range(0, nTurbs):

        # position[turbI, 0] = turbineX[turbI]
        # position[turbI, 1] = turbineY[turbI]

        rotorDiameter[turbI] = 126.4
        rotorArea[turbI] = 0.25*np.pi*rotorDiameter[turbI]*rotorDiameter[turbI]
        axialInduction[turbI] = 1.0/3.0
        Ct[turbI] = 4.0*axialInduction[turbI]*(1.0-axialInduction[turbI])
        Cp[turbI] = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
        generator_efficiency[turbI] = 0.944
        yaw[turbI] = 0.
        # yaw[turbI] = 0

    # myFloris.position = position
    myFloris.turbineX = turbineX
    myFloris.turbineY = turbineY
    # myFloris.floris_windframe.turbineX = turbineX
    # myFloris.floris_windframe.turbineY = turbineY
    myFloris.rotorDiameter = rotorDiameter
    myFloris.rotorArea = rotorArea
    myFloris.axialInduction = axialInduction
    myFloris.Ct = Ct
    myFloris.Cp = Cp
    myFloris.generator_efficiency = generator_efficiency
    # myFloris.yaw = yaw
    myFloris.yaw = yaw

    # Define flow properties
    myFloris.wind_speed = 8.0  # m/s
    myFloris.air_density = 1.1716  # kg/m^3
    myFloris.wind_direction = 30  # deg
    myFloris.verbose = False

    myFloris.parameters.CPcorrected = False
    myFloris.parameters.CTcorrected = False

    # define sampling points for optimization run (GenericFlowModel)
    # resolution = 0.
    # xlim = [0, 3000]
    # ylim = [0, 3000]
    # x = np.linspace(xlim[0], xlim[1], resolution)
    # y = np.linspace(ylim[0], ylim[1], resolution)
    # x, y = np.meshgrid(x, y)
    # myFloris.ws_position = np.array([x.flatten(), y.flatten()]).transpose()
    #print myFloris.ws_position
    # run model
    print 'start FLORIS Opt. run'
    tic = time.time()
    myFloris.run()
    toc = time.time()
    print('FLORIS Opt. calculation took %.03f sec.' % (toc-tic))

    # Display returns
    # print 'turbine powers (kW): %s' % myFloris2.wt_power
    # print 'turbine powers (kW): %s' % myFloris2.wt_power
    print 'turbine powers (kW): %s' % myFloris.wt_power
    print 'turbine X positions in wind frame (m): %s' % myFloris.floris_windframe.turbineX
    print 'turbine Y positions in wind frame (m): %s' % myFloris.floris_windframe.turbineY
    # print 'yaw (deg) = ', myFloris.yaw
    print 'yaw (deg) = ', myFloris.floris_adjustCtCp.yaw
    # print 'yaw (deg) = ', myFloris.floris_adjustCtCp.yaw
    print 'effective wind speeds (m/s): %s' % myFloris.floris_power.velocitiesTurbines
    # print 'Wake center Y positions (m): %s' % myFloris2.wakeCentersYT
    # print 'Wake diameters (m): %s' % myFloris2.wakeDiametersT
    # print 'Relative wake overlap (m*m): %s' % myFloris2.wakeOverlapTRel
    print 'wind farm power (kW): %s' % myFloris.power

    # windDirection = myFloris.wind_direction
    # rotationMatrix = np.array([(np.cos(windDirection), -np.sin(windDirection)),
    #                                (np.sin(windDirection), np.cos(windDirection))])
    # # print 'rotation matrix = ', rotationMatrix
    # turbinepos_opt = np.dot(rotationMatrix, np.array([myFloris.turbineXw, myFloris.turbineYw]))
    # print turbinepos_opt
    #
    myFloris2 = floris_assembly()
    myFloris2.turbineX = myFloris.floris_windframe.turbineX
    myFloris2.turbineY = myFloris.floris_windframe.turbineY

    myFloris2.rotorDiameter = rotorDiameter
    myFloris2.rotorArea = rotorArea
    myFloris2.axialInduction = axialInduction
    myFloris2.Ct = Ct
    myFloris2.Cp = Cp
    myFloris2.generator_efficiency = generator_efficiency
    myFloris2.yaw = myFloris.floris_adjustCtCp.yaw

    # # Define flow properties
    myFloris2.wind_speed = 8.0  # m/s
    myFloris2.air_density = 1.1716  # kg/m^3
    myFloris2.wind_direction = 30  # deg
    # myFloris2.verbose = True

    myFloris2.parameters.CPcorrected = False
    myFloris2.parameters.CTcorrected = False

    # define sampling points for final visualization (GenericFlowModel)
    resolution = 200
    xlim = [0, 3000]
    ylim = [0, 3000]
    x = np.linspace(xlim[0], xlim[1], resolution)
    y = np.linspace(ylim[0], ylim[1], resolution)
    x, y = np.meshgrid(x, y)
    myFloris2.ws_position = np.array([x.flatten(), y.flatten()]).transpose()
    #print myFloris.ws_position
    # run model
    print 'start FLORIS Vis. run'
    tic = time.time()
    myFloris2.run()
    toc = time.time()
    print('FLORIS Vis. calculation took %.03f sec.' % (toc-tic))



    np.set_printoptions(linewidth=150, precision=4)

    # Display returns
    # print 'turbine powers (kW): %s' % myFloris2.wt_power
    # print 'turbine powers (kW): %s' % myFloris2.wt_power
    print 'turbine powers (kW): %s' % myFloris2.wt_power
    # print 'turbine X positions in wind frame (m): %s' % myFloris2.turbineXw
    # print 'turbine Y positions in wind frame (m): %s' % myFloris2.turbineYw
    # print 'Wake center Y positions (m): %s' % myFloris2.wakeCentersYT
    # print 'Wake diameters (m): %s' % myFloris2.wakeDiametersT
    # print 'Relative wake overlap (m*m): %s' % myFloris2.wakeOverlapTRel
    print 'wind farm power (kW): %s' % myFloris2.power
    # print myFloris2.ws_array
    velocities = myFloris2.ws_array.reshape(resolution, resolution)

    import matplotlib.pyplot as plt

    fig, (ax1) = plt.subplots(nrows=1)
    im = ax1.pcolormesh(x, y, velocities, cmap='coolwarm')
    plt.colorbar(im, orientation='vertical')
    ax1.set_aspect('equal')
    ax1.autoscale(tight=True)
    plt.show()