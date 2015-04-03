#!/usr/bin/env python
# encoding: utf-8
"""

- This script illustrates an example call to the FLORIS model with some
fusedwind classes used

Created by Pieter Gebraad, 2015
"""

from floris import FLORIS
import numpy as np
from florisRotorSE import AeroelasticHAWTVT_floris
import time

myFloris = FLORIS()

# define all turbine properties
turbineX = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1])
turbineY = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7])

for turbI in range(0,turbineX.size):
    wt = AeroelasticHAWTVT_floris()
    wt.turbine_name = 'NREL 5-MW baseline turbine'
    wt.name = 'turbine%s' % turbI
    wt.position = np.array([turbineX[turbI], turbineY[turbI]])
    wt.rotor_diameter = 126.4
    wt.rotor_area = 0.25*np.pi*wt.rotor_diameter*wt.rotor_diameter
    wt.axial_induction = 1.0/3.0
    wt.CP = 0.7737 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
    wt.yaw = 0
    myFloris.wt_layout.add_wt(wt)

# Define flow properties
myFloris.wind_speed = 8  # m/s
myFloris.air_density = 1.1716  # kg/m^3
myFloris.wind_direction = 30  # deg

# define sampling points (GenericFlowModel)
resolution = 10
x = np.linspace(0,3000,resolution)
y = np.linspace(0,3000,resolution)
x, y = np.meshgrid(x, y)
myFloris.ws_positions = np.array([x.flatten(),
                                  y.flatten()]).transpose()
# run model
print 'start FLORIS run'
tic = time.time()
myFloris.run()
toc = time.time()
print('FLORIS calculation took %.03f sec.' % (toc-tic))

np.set_printoptions(linewidth=150,precision=4)

# Display returns
print 'turbine powers (kW): %s' % myFloris.wt_power
print 'velocity sampling (m/s):\n %s' % myFloris.ws_array.reshape(resolution, resolution)