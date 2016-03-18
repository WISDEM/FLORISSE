#!/usr/bin/env python
# encoding: utf-8

"""
OptimizationGroups.py
Created by Jared J. Thomas, Nov. 2015.
Brigham Young University
"""

import numpy as np

from openmdao.api import Group, IndepVarComp, ExecComp

from floris import DirectionGroupFLORIS, AEPGroupFLORIS
from GeneralWindFarmComponents import SpacingComp
from Parameters import FLORISParameters


class OptPowerOneDir(Group):
    """ Group connecting the floris model for optimization with one wind direction"""

    def __init__(self, nTurbines, resolution=0, minSpacing=2., differentiable=True, use_rotor_components=True):

        super(OptPowerOneDir, self).__init__()

        # add major components
        self.add('dirComp', AEPGroupFLORIS(nTurbines, resolution=0, differentiable=differentiable,
                                           use_rotor_components=use_rotor_components), promotes=['*'])
        self.add('spacing_comp', SpacingComp(nTurbines=nTurbines), promotes=['*'])

        # add constraint definitions
        self.add('spacing_con', ExecComp('sc = separation_squared-(minSpacing*rotorDiameter[0])**2',
                                         minSpacing=minSpacing, rotorDiameter=np.zeros(nTurbines),
                                         sc=np.zeros(((nTurbines-1.)*nTurbines/2.)),
                                         separation_squared=np.zeros(((nTurbines-1.)*nTurbines/2.))),
                 promotes=['*'])

        # add objective component
        self.add('obj_comp', ExecComp('obj = -1.*power0', power0=0.0), promotes=['*'])

        # initialize design variables for optimization
        # self.add('p1', IndepVarComp('turbineX', np.zeros(nTurbines)), promotes=['*'])
        # self.add('p2', IndepVarComp('turbineY', np.zeros(nTurbines)), promotes=['*'])
        # self.add('p3', IndepVarComp('yaw', np.zeros(nTurbines)), promotes=['*'])


class OptAEP(Group):
    """
        Group adding optimization parameters to an AEPGroup


        ----------------
        Design Variables
        ----------------
        turbineX:   1D numpy array containing the x coordinates of each turbine in the global reference frame
        turbineY:   1D numpy array containing the x coordinates of each turbine in the global reference frame
        yaw_i:      1D numpy array containing the yaw angle of each turbine in the wind direction reference frame for
                    direction i

        ---------------
        Constant Inputs
        ---------------
        rotorDiameter:                          1D numpy array containing the rotor diameter of each turbine

        axialInduction:                         1D numpy array containing the axial induction of each turbine. These
                                                values are not actually used unless the appropriate floris_param is set.

        generator_efficiency:                   1D numpy array containing the efficiency of each turbine generator

        wind_speed:                             scalar containing a generally applied inflow wind speed

        air_density:                            scalar containing the inflow air density

        windDirections:                         1D numpy array containing the angle from N CW to the inflow direction

        windrose_frequencies:                   1D numpy array containing the probability of each wind direction

        Ct:                                     1D numpy array containing the thrust coefficient of each turbine

        Cp:                                     1D numpy array containing the power coefficient of each turbine

        floris_params:FLORISoriginal(False):    boolean specifying which formulation of the FLORIS model to use. (True
                                                specfies to use the model as originally formulated and published).

        floris_params:CPcorrected(True):        boolean specifying whether the Cp values provided have been adjusted
                                                for yaw

        floris_params:CTcorrected(True):        boolean specifying whether the Ct values provided have been adjusted
                                                for yaw

        -------
        Returns
        -------
        AEP:                scalar containing the final AEP for the wind farm

        power_directions:   1D numpy array containing the power production for each wind direction (unweighted)

        velocitiesTurbines: 1D numpy array of velocity at each turbine in each direction. Currently only accessible by
                            *.AEPgroup.dir%i.unknowns['velocitiesTurbines']

        wt_powers: 1D numpy array of power production at each turbine in each direction. Currently only accessible by
                            *.AEPgroup.dir%i.unknowns['velocitiesTurbines']

    """

    def __init__(self, nTurbines, resolution=0, nDirections=1, minSpacing=2., use_rotor_components=True,
                 datasize=0, differentiable=True, optimizingLayout=False):

        super(OptAEP, self).__init__()

        # add major components and groups
        self.add('AEPgroup', AEPGroupFLORIS(nTurbines=nTurbines, nDirections=nDirections,
                                            use_rotor_components=use_rotor_components,
                                            datasize=datasize, differentiable=differentiable,
                                            optimizingLayout=optimizingLayout),
                 promotes=['*'])

        self.add('spacing_comp', SpacingComp(nTurbines=nTurbines), promotes=['*'])

        # add constraint definitions
        self.add('spacing_con', ExecComp('sc = separation_squared-(minSpacing*rotorDiameter[0])**2',
                                         minSpacing=minSpacing, rotorDiameter=np.zeros(nTurbines),
                                         sc=np.zeros(((nTurbines-1.)*nTurbines/2.)),
                                         separation_squared=np.zeros(((nTurbines-1.)*nTurbines/2.))),
                 promotes=['*'])

        # add objective component
        self.add('obj_comp', ExecComp('obj = -1.*AEP', AEP=0.0), promotes=['*'])
