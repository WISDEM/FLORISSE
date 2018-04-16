import numpy as np
from openmdao.api import Group, Problem, IndepVarComp, pyOptSparseDriver, view_tree, profile
# from FLORISSE3D.floris import Floris as FLORIS3D
from florisse.floris import Floris as floris

if __name__ == '__main__':
    """setup the turbine locations"""

    rotor_diameter = 126.4
    turbineX = np.array([0.,200.,400.,600.])
    turbineY = np.array([0.,0.,0.,0.])
    turbineZ = np.array([100.,150.,45.,55.])

    nTurbs = len(turbineX)

    rotorDiameter = np.ones(nTurbs)*100.

    nDirections = 1

    shearExp = 0.15
    z0 = 0.
    zref = 50.
    Uref = 8.

    windArray = np.zeros(nTurbs)
    for turbine_id in range(nTurbs):
            windArray[turbine_id] = Uref*((turbineZ[turbine_id]-z0)/(zref-z0))**shearExp


    """OpenMDAO"""

    prob = Problem()
    root = prob.root = Group()

    root.add('floris', floris(nTurbs, differentiable=True, use_rotor_components=False, nSamples=0,
                 verbose=False),promotes=['turbineXw','turbineYw','hubHeight','rotorDiameter','wind_speed',
                 'floris_params:shearExp','floris_params:z_ref'])
    # root.add('FLORIS', FLORIS3D(nTurbs, differentiable=True, use_rotor_components=False, nSamples=0,
    #              verbose=False),promotes=['turbineXw','turbineYw','turbineZ','rotorDiameter'])


    prob.root.ln_solver.options['single_voi_relevance_reduction'] = True
    profile.setup(prob)
    profile.start()
    prob.setup(check=True)

    prob['turbineXw'] = turbineX
    prob['turbineYw'] = turbineY
    prob['hubHeight'] = turbineZ
    # prob['turbineZ'] = turbineZ
    # prob['FLORIS.wind_speed'] = windArray
    prob['floris_params:z_ref'] = zref
    prob['floris_params:shearExp'] = shearExp

    prob.run()
    print prob['floris.wtVelocity0']
    # print prob['FLORIS.wtVelocity0']
    #
    #
    # prob['turbineY'] = turbineY
    #
    # prob.run()
