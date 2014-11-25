#!/usr/bin/env python
# encoding: utf-8
"""
example.py

Created by Andrew Ning on 2014-04-22.
Copyright (c) NREL. All rights reserved.

Some updates by Paul Fleming April 30, 2014
-Assigning effU
-Re-arranging

-This script illustrates an example call to the FLORIS model

"""


import sys
sys.path.insert(0,'src')
import _floris
import numpy as np
import math

Cp = 0.7737 * 4.0 * 1.0/3.0 * math.pow((1 - 1.0/3.0),2)
PI = 3.1415926535897932385;
yawVal = 0 # Now defined realtive to inflow30.0 * PI/180.0

# Defube turbine locations and orientation
X = np.array([1164.7, 947.2,  1682.4, 1464.9, 1982.6, 2200.1 ])
Y = np.array([1024.7, 1335.3, 1387.2, 1697.8, 2060.3, 1749.7 ])
yaw = np.array([yawVal, yawVal, yawVal, yawVal, yawVal, yawVal])

# Define turbine characteristics
axialInd = np.array([1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0])
rotorDiameter = np.array([126.4, 126.4, 126.4, 126.4, 126.4, 126.4])
Cp = np.array([Cp,Cp,Cp,Cp,Cp,Cp])

# Define turbine measurements
measPower = np.array([100.0, 100.0, 100.0, 100.0, 100.0, 100.0])
effU = np.array([8,8,8,8,8,8])
useMeasPower = False

# Define site measurements
effUdXY = 0.523599
rho = 1.1716

# Call FLORIS
power, effU = _floris.floris(X, Y, axialInd, yaw, rotorDiameter,
                             Cp, measPower, effU, effUdXY, rho, useMeasPower)

# Display retunrs
print power
print effU
