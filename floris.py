import numpy as np

from openmdao.api import Group, Component, Problem, IndepVarComp


class FLORIS(Group):
    stuff = 0


class floris_adjustCtCp(Component):

    """ Adjust Cp and Ct to yaw if they are not already adjusted """

    def __init__(self, nTurbines):

        # print 'entering adjustCtCp __init__ - analytic'
        super(floris_adjustCtCp, self).__init__()

        # Explicitly size input arrays
        self.add_param('Ct_in', val=np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_param('Cp_in', val=np.zeros(nTurbines), desc='power coefficient for all turbines')
        self.add_param('generator_efficiency', val=np.zeros(nTurbines), desc='generator efficiency of all turbines')
        self.add_param('yaw', val=np.zeros(nTurbines), desc='yaw of each turbine')

        # Explicitly size output arrays
        self.add_output('Ct_out', val=np.zeros(nTurbines), desc='Thrust coefficient for all turbines')
        self.add_output('Cp_out', val=np.zeros(nTurbines), desc='power coefficient for all turbines')

        # parameters since var trees are not supports
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


    def solve_nonlinear(self, params, unknowns, resids):

        # print 'entering adjustCtCP - analytic'

        # collect inputs
        Ct = params['Ct_in']
        Cp = params['Cp_in']
        yaw = params['yaw']*np.pi/180.

        # determine floris_parameter values
        if params['floris_params:FLORISoriginal']:
            pP = 1.88
        else:
            pP = params['floris_params:pP']

        CTcorrected = params['floris_params:CTcorrected']
        CPcorrected = params['floris_params:CPcorrected']

        # calculate new CT values, if desired
        if not CTcorrected:
            unknowns['Ct_out'] = Ct*np.cos(yaw)*np.cos(yaw)
        else:
            unknowns['Ct_out'] = Ct

        # calculate new CP values, if desired
        if not CPcorrected:
            unknowns['Cp_out'] = Cp*np.cos(yaw)**pP
        else:
            unknowns['Cp_out'] = Cp

    def jacobian(self, params, unknowns, resids):

        # collect inputs
        Ct = params['Ct_in']
        Cp = params['Cp_in']
        nTurbines = np.size(Ct)
        yaw = params['yaw']*np.pi/180.

        # determine floris_parameter values
        if self.parameters.FLORISoriginal:
            pP = 1.88
        else:
            pP = params['floris_params:pP']

        CTcorrected = params['floris_params:CTcorrected']
        CPcorrected = params['floris_params:CPcorrected']

        # calculate gradients
        J = {}

        if not CTcorrected:
            J[('Ct_out', 'Ct_in')] = np.eye(nTurbines)*np.cos(yaw)*np.cos(yaw)
            J[('Ct_out', 'Cp_in')] = np.zeros((nTurbines, nTurbines))
            J[('Ct_out', 'yaw')] = np.eye(nTurbines)*(-2.*Ct*np.sin(yaw)*np.cos(yaw))*np.pi/180.
        else:
            J[('Ct_out', 'Ct_in')] = np.eye(nTurbines, nTurbines)
            J[('Ct_out', 'Cp_in')] = np.zeros((nTurbines, nTurbines))
            J[('Ct_out', 'yaw')] = np.zeros((nTurbines, nTurbines))

        if not CPcorrected:
            J[('Cp_out', 'Cp_in')] = np.eye(nTurbines, nTurbines)*np.cos(yaw)**pP
            J[('Cp_out', 'Ct_in')] = np.zeros((nTurbines, nTurbines))
            J[('Cp_out', 'yaw')] = np.eye(nTurbines, nTurbines)*(-Cp*pP*np.sin(yaw)*np.cos(yaw)**(pP-1.0))*np.pi/180.
        else:
            J[('Cp_out', 'Cp_in')] = np.eye(nTurbines, nTurbines)
            J[('Cp_out', 'Ct_in')] = np.zeros((nTurbines, nTurbines))
            J[('Cp_out', 'yaw')] = np.zeros((nTurbines, nTurbines))

        return J

class floris_windframe(Component):
    stuff = 0

class floris_wcent_wdiam(Component):
    stuff = 0

class floris_overlap(Component):
    stuff = 0

class floris_power(Component):
    stuff = 0

if __name__ == "__main__":

    top = Problem()

    root = top.root = Group()

    root.add('p1', IndepVarComp('x', np.array([3.0])))
    root.add('p2', IndepVarComp('y', np.array([2.0])))
    root.add('p3', IndepVarComp('z', np.array([1.0])))
    root.add('p', floris_adjustCtCp(nTurbines=np.array([1])))

    root.connect('p1.x', 'p.Ct_in')
    root.connect('p2.y', 'p.Cp_in')
    root.connect('p3.z', 'p.yaw')

    top.setup()
    top.run()

    print(root.p.unknowns['Ct_out'])
    print(root.p.unknowns['Cp_out'])
