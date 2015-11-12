import numpy as np

from openmdao.api import Group, IndepVarComp, ExecComp

from floris_openmdao1 import FLORIS, DirectionGroupFLORIS, AEPGroupFLORIS
from GeneralWindfarmComponents import AdjustCtCpYaw, WindFarmAEP, MUX, DeMUX


class OptPowerOneDir(Group):
    """ Group connecting the floris model and adjustCtCp for optimization """

    def __init__(self, nTurbines, resolution=0):

        super(OptPowerOneDir, self).__init__()

        # add major components
        self.add('myFloris', DirectionGroupFLORIS(nTurbines, resolution=0), promotes=['*'])

        # add objective component
        self.add('obj_comp', ExecComp('obj = -1*power', power=0.0), promotes=['*'])
        #
        # # connect components within the problem
        # self.connect('Ct_out', 'Ct')
        # self.connect('Cp_out', 'Cp')

         # initialize design variables for optimization
        self.add('p1', IndepVarComp('turbineX', np.zeros(nTurbines)), promotes=['*'])
        self.add('p2', IndepVarComp('turbineY', np.zeros(nTurbines)), promotes=['*'])
        self.add('p3', IndepVarComp('yaw', np.zeros(nTurbines)), promotes=['*'])


class OptAEP(Group):
    """ Group connecting the floris model and adjustCtCp for optimization """

    def __init__(self, nTurbines, resolution=0, nDirections=1):

        super(OptAEP, self).__init__()

        # add major components
        self.add('AEPgroup', AEPGroupFLORIS(nTurbines=nTurbines, nDirections=nDirections), promotes=['*'])

        # add objective component
        self.add('obj_comp', ExecComp('obj = -1*AEP', AEP=0.0), promotes=['*'])


# class SetupAEP(Group):
#     """ Group connecting the floris model and adjustCtCp for optimization """
#
#     def __init__(self, nTurbines, resolution=0, nDirections=1):
#
#         super(OptAEP, self).__init__()
#
#         # add major components
#         for i in range(0, nDirections):
#             self.add('CtCp%i' % i, AdjustCtCpYaw(nTurbines), promotes=['Ct_in', 'Cp_in', 'params:*'])
#             self.add('myFloris%i' % i, FLORIS(nTurbines, resolution=resolution),
#                      promotes=['floris_params:*', 'wind_speed', 'air_density', 'generator_efficiency',
#                                'turbineX', 'turbineY', 'rotorDiameter'])
#         self.add('powerMUX', MUX(nDirections))
#         self.add('AEPcomp', WindFarmAEP(nDirections), promotes=['*'])
#
#         # connect components
#         for i in range(0, nDirections):
#             self.connect('CtCp%i.Ct_out' % i, 'myFloris%i.Ct' % i)
#             self.connect('CtCp%i.Cp_out' % i, 'myFloris%i.Cp' % i)
#             self.connect('myFloris%i.power' % i, 'powerMUX.input%i' % i)
#             self.connect('myFloris%i.wind_direction' % i, 'powerMUX.input%i' % i)
#         self.connect('powerMUX.Array', 'power_directions')
#         self.connect('floris_params:CTcorrected', 'params:CTcorrected')
#         self.connect('floris_params:CPcorrected', 'params:CPcorrected')
#
#         # add objective component
#         self.add('obj_comp', ExecComp('obj = -1*AEP', AEP=0.0), promotes=['*'])
#
#          # initialize design variables for optimization
#         self.add('p1', IndepVarComp('turbineX', np.zeros(nTurbines)), promotes=['*'])
#         self.add('p2', IndepVarComp('turbineY', np.zeros(nTurbines)), promotes=['*'])
#
#         for i in range(3, nDirections+3):
#             self.add('p%i' % i, IndepVarComp('yaw%i' % (i-3), np.zeros(nTurbines)), promotes=['*'])

