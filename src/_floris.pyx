# distutils: language = c

import numpy as np
cimport numpy as np


cdef extern from "FLORISmodel.h":

    void FLORISmodel(double* power, double* effU, double* effU_In, double* X, double* Y, double* axialInd,
        double* yaw, double* rotorDiameter, double* Cp, double* measPower, int nTurb,
        double effUdXY, double rho,
        double kd, double ad, double bd, double ke, double me1, double me2, double me3,
        double pP, double MU1, double MU2, double MU3, double aU, double bU,
        int flagUseMeasPower)



def floris(X, Y, axialInd, yaw, rotorDiameter, Cp, measPower,
           effU_In, effUdXY, rho, useMeasPower=True,
           kd=0.15, ad=-4.5, bd=-0.01, ke=0.065, me1=-0.5, me2=0.22,
           me3=1.0, pP=1.88, MU1=0.5, MU2=1.0, MU3=5.5, aU=5.0, bU=1.66):

    nTurb = len(X)
    power = np.zeros(nTurb)
    effU = np.zeros(nTurb)

    cdef np.ndarray[np.double_t, ndim=1] power_c
    power_c = np.ascontiguousarray(power, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] effU_c
    effU_c = np.ascontiguousarray(effU, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] effU_In_c
    effU_In_c = np.ascontiguousarray(effU_In, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] X_c
    X_c = np.ascontiguousarray(X, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] Y_c
    Y_c = np.ascontiguousarray(Y, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] axialInd_c
    axialInd_c = np.ascontiguousarray(axialInd, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] yaw_c
    yaw_c = np.ascontiguousarray(yaw, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] rotorDiameter_c
    rotorDiameter_c = np.ascontiguousarray(rotorDiameter, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] Cp_c
    Cp_c = np.ascontiguousarray(Cp, dtype=np.double)

    cdef np.ndarray[np.double_t, ndim=1] measPower_c
    measPower_c = np.ascontiguousarray(measPower, dtype=np.double)

    FLORISmodel(&power_c[0], &effU_c[0], &effU_In_c[0], &X_c[0], &Y_c[0], &axialInd_c[0], &yaw_c[0],
        &rotorDiameter_c[0], &Cp_c[0], &measPower_c[0], nTurb, effUdXY, rho,
        kd, ad, bd, ke, me1, me2, me3, pP, MU1, MU2, MU3, aU, bU, useMeasPower)


    return power, effU