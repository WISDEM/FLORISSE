FLORISSE ELLIPSE
========


Created by Pieter Gebraad.
Copyright (c) NREL. All rights reserved.

-This archive contains the OpenMDAO 0.13 implementation of FLORIS, with adjustments to incorporate veer effects
-This model also implements downstream corrections close to the rotor using an arctan profile, and a shear profile
-For overview of this FLORIS model refer to the pptx file

-The main examples to run are fitSkewedWakes.py and ICOWES3_veerProfiles_manualTune.py

fitSkewedWakes.py fits to the average SOWFA profiles generated for Matts poster:
"Large-Eddy Simulations of Wind Turbine Wakes  Subject to Different Atmospheric Stabilities"

ICOWES3_veerProfiles_manualTune.py fits to averaged SOWFA ICOWES3 profiles ( generated with ICOWES3_averageProfiles.py ) using the neutral and stable case, and performs power prediction for those two cases. Note that:
* a linear relationship is found between the slope of the 5D upstream cross-stream wind across the rotor (as a measure for veer), to the coefficient wakeSkewCoefficient  that defines the increase in in-plane skew with downstream distance (note that more cases with more or less veer are needed to confirm that first relationship is linear). The linear relationship between downstream distance and skew is built-in into FLORIS, although also a constant skew can be predefined (wakeSkew)
* Adjustments had to be made to CP and CT curve to make a good prediction of power in shear. Note that CCBlade should have some option to correct for veer (but assumes a power-law).
* Velocity shear profile in FLORIS can be predefined as a scaling factor that changes with height, normalized by hub-height (shearProfileZnorm, shearProfileUnorm), or assuming a power-law w
ith coefficient shearCoefficientAlpha
* We use downstream corrections close to the rotor using an arctan profile, k_nearwake could be used to adjust somewhat.



