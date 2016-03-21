from openmdao.main.api import VariableTree
from openmdao.lib.datatypes.api import Array, Bool, Float
import numpy as np

class FLORISParameters(VariableTree):
    """Container of FLORIS wake parameters"""

    # original tuning parameters
    pP = Float(1.88, iotype='in') # yaw power correction parameter
    ke = Float(0.05, iotype='in') # wake expansion paramters
    keCorrArray = Float(0.0, iotype='in') # array-correction factor
    keCorrCT = Float(0.0, iotype='in') # CT-correction factor
    baselineCT = Float(4.0*(1.0/3.0)*(1.0-(1.0/3.0)), iotype='in') # Baseline CT for ke-correction
    me = Array(np.array([-0.5, 0.22, 1.0]), iotype='in') # relative expansion of wake zones

    kd = Float(0.17, iotype='in') # wake deflection recovery factor
    
    # define initial wake displacement and angle (not determined by yaw angle)
    useWakeAngle = Bool(True, iotype = 'in')
    initialWakeDisplacement = Float(-4.5, iotype='in')
    initialWakeAngle = Float(1.5, iotype='in')
    bd = Float(-0.01, iotype='in')

    # correction recovery coefficients with yaw
    useaUbU = Bool(True, iotype = 'in')
    aU = Float(12.0, iotype='in', units='deg')
    bU = Float(1.3, iotype='in')

    MU = Array(np.array([0.5, 1.0, 5.5]), iotype='in')
    CTcorrected = Bool(True, iotype='in', desc = 'CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
    CPcorrected = Bool(True, iotype='in', desc = 'CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)')
    axialIndProvided = Bool(False, iotype='in', desc = 'CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)') 

    # adjust initial wake diameter to yaw
    adjustInitialWakeDiamToYaw = Bool(False, iotype = 'in')

    # shear layer (only influences visualization)
    shearCoefficientAlpha = Float(0.10805, iotype='in')
    shearZh = Float(90., iotype='in')
