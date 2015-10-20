from openmdao.api import Component
import numpy as np

class FLORISParameters(Component):
    """Container of FLORIS wake parameters"""

    def __init_(self):
        # original tuning parameters
        self.add_param('floris_params:pP', 1.88)
        self.add_param('floris_params:ke', 0.065)
        self.add_param('floris_params:keCorrArray', 0.0)
        self.add_param('floris_params:keCorrCT', 0.0)
        self.add_param('floris_params:Region2CT', 4.0*(1.0/3.0)*(1.0-(1.0/3.0)))
        self.add_param('floris_params:kd', 0.15)
        self.add_param('floris_params:me', np.array([-0.5, 0.22, 1.0]))

        self.add_param('floris_params:initialWakeDisplacement', -4.5)
        self.add_param('floris_params:initialWakeAngle', 3.0)

        self.add_param('floris_params:baselineCT', 4./3.*(1.-1./3.))

        self.add_param('floris_params:keCorrTI', 0.0)
        self.add_param('floris_params:baselineTI', 0.045)

        self.add_param('floris_params:keCorrHR', 0.0) # neutral, with heating rate 0, is baseline

        self.add_param('floris_params:keCorrHRTI', 0.0)

        self.add_param('floris_params:keSaturation', 0.0)

        self.add_param('floris_params:kdCorrYawDirection', 0.0)


        self.add_param('floris_params:MU', np.array([0.5, 1.0, 10]))

        self.add_param('floris_params:CTcorrected', True, desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')
        self.add_param('floris_params:CPcorrected', True, desc = 'CP factor already corrected by CCBlade calculation (assumed with approximately factor cos(yaw)^3)')

        self.add_param('floris_params:axialIndProvided', False, desc='CT factor already corrected by CCBlade calculation (approximately factor cos(yaw)^2)')


        self.add_param('floris_params:useWakeAngle', True)

        self.add_param('floris_params:bd', -0.01)

        self.add_param('floris_params:useaUbU', False)
        self.add_param('floris_params:aU', 5.0, units='deg')
        self.add_param('floris_params:bU', 1.66)

        self.add_param('floris_params:adjustInitialWakeDiamToYaw', True)

        self.add_param('floris_params:FLORISoriginal', False, desc='override all parameters and use FLORIS as original in first Wind Energy paper')
