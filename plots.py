import numpy as np
from matplotlib import pyplot as plt
from scipy.io import loadmat
from openmdao.api import Problem

from floris_openmdao1 import AEPGroupFLORIS


def writetofile(towrite, filename):
    # print 'entered write func'
    #open file to write to
    txtfile = open(filename, 'w')

    # print the results one element at a time
    if np.size(towrite) == 1:
        txtfile.write(str(towrite))
    else:
        for i in range(np.size(towrite)):

            txtfile.write(str(towrite[i]))
            txtfile.write('\n')

    # close file
    txtfile.close()


def readfile(filename,n):
    # print 'entered write func'
    #open file to write to
    txtfile = open(filename, 'r')
    contents = np.zeros(n)
    for i in range(n):
        contents[i] = float(txtfile.readline())

    # close file
    txtfile.close()

    return contents

# NREL5MWCPCT = pickle.load(open('NREL5MWCPCT.p'))
# datasize = NREL5MWCPCT.CP.size
nTurbines = 2
myFloris = Problem(root=AEPGroupFLORIS(nTurbines=nTurbines, nDirections=1, resolution=0.0))
myFloris.setup()

rotorDiameter = 126.4
rotorArea = np.pi*rotorDiameter*rotorDiameter/4.0
axialInduction = 1.0/3.0
CP = 0.7737/0.944 * 4.0 * 1.0/3.0 * np.power((1 - 1.0/3.0), 2)
CT = 4.0*axialInduction*(1.0-axialInduction)
generator_efficiency = 0.944

myFloris['floris_params:CPcorrected'] = False
myFloris['floris_params:CTcorrected'] = False
myFloris['floris_params:FLORISoriginal'] = False

# Define turbine characteristics
myFloris['axialInduction'] = np.array([axialInduction, axialInduction])
myFloris['rotorDiameter'] = np.array([rotorDiameter, rotorDiameter])

myFloris['generator_efficiency'] = np.array([generator_efficiency, generator_efficiency])

# Define site measurements
myFloris['windDirections'] = np.array([240.])
myFloris['wind_speed'] = 8.    # m/s
myFloris['air_density'] = 1.1716

myFloris['Ct_in'] = np.array([CT, CT])
myFloris['Cp_in'] = np.array([CP, CP])

FLORISeffu = list()
res = 10000
y = np.linspace(-1.5*rotorDiameter, 1.5*rotorDiameter, res)
myFloris['rotorDiameter'] = np.array([rotorDiameter, 0.1*rotorDiameter])
for i in range(0, res):
    myFloris['windDirections'] = np.array([270])
    X = np.array([0, 10.*rotorDiameter])
    Y = np.array([0, y[i]])
    myFloris['turbineX'] = X
    myFloris['turbineY'] = Y
    myFloris['yaw0'] = np.array([0.0, 0.0])

    # Call FLORIS
    myFloris.run()

    FLORISeffu.append(list(myFloris.root.dir0.unknowns['velocitiesTurbines']))

FLORISeffu = np.array(FLORISeffu)
# writetofile(FLORISeffu[:, 1],'effu_discon_p1diam.txt')
FLORISeffu_disc = readfile('effu_discon_p1diam.txt', 10000)
plt.figure()
plt.plot(y/rotorDiameter, FLORISeffu[:, 1], label='Differentiable')
plt.plot(y/rotorDiameter, FLORISeffu_disc, label='Original')
plt.legend()
plt.xlabel('cross-wind rotor/wake center separation (rotor diams)')
plt.ylabel('effective hub velocity (m/s)')
plt.show()

#
# FLORISindiam = list()
# res = 100
# x = np.linspace(-0.25*rotorDiameter, 5.0*rotorDiameter, res)
# for i in range(0, res):
#     myFloris['windDirections'] = np.array([270])
#     X = np.array([0, x[i]])
#     Y = np.array([0, myFloris['wakeCentersYT'][2]])
#     print myFloris['wakeCentersYT'][2]
#     # Y = np.array([0, 0])
#     myFloris['turbineX'] = X
#     myFloris['turbineY'] = Y
#     myFloris['yaw0'] = np.array([0.0, 0.0])
#
#     # Call FLORIS
#     myFloris.run()
#
#     FLORISindiam.append(list(myFloris.root.dir0.unknowns['velocitiesTurbines']))
#
# FLORISindiam = np.array(FLORISindiam)
# plt.figure()
# plt.plot(x/rotorDiameter, FLORISindiam[:, 1])
# plt.show()




