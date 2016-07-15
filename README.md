FLORISSE ELLIPSE
========


Created by Pieter Gebraad.
Copyright (c) NREL. All rights reserved.

-This archive contains the OpenMDAO 0.13 implementation of FLORIS, with adjustments to incorporate veer effects, based on rotated elliptical shapes of the wake
-This model also implements downstream corrections close to the rotor using an arctan profile, and a shear profile
-For overview of this FLORIS model refer to the pptx file

- Uses the C-code EEOver for calculating overlaps of ellipses, https://github.com/chraibi/EEOver, I made a Python wrapper for it

-- INSTALLATION

- Git clone this branch
- Setup openmdao 0.13 and source it
- For Window only : EEOver needs Gnu Scientific Library, which is standard in Linux, but you need to download them for Windows:

  1) Download GSL Windows DLL and headers for either 32 and 64-bit from https://code.google.com/archive/p/oscats/downloads, For example, gsl-1.15-dev-win64.zip for the 64-bit.

  2) Place contents of zip-file in folder called GSL in the EEOver directory, such that EEOver/GSL/lib, EEOver/GSL/include, EEOver/GSL/bin are there.
  
- Go into EEOver and build this part of the code: run "python setup.py build_ext --inplace" for Linux or "python setupWindows.py build_ext --inplace" for Windows
- Test EEOver by calling "python callEE.py"

-- EXAMPLES

The main examples to run are fitSkewedWakes.py and ICOWES3_veerProfiles_manualTune.py

fitSkewedWakes.py fits to the average SOWFA profiles generated for Matts poster:
"Large-Eddy Simulations of Wind Turbine Wakes  Subject to Different Atmospheric Stabilities"

ICOWES3_veerProfiles_manualTune.py fits to averaged SOWFA ICOWES3 profiles (generated with ICOWES3_averageProfiles.py) using the neutral and stable case, and performs power prediction for those two cases. Note that:

* a linear relationship is found between the slope of the 5D upstream cross-stream wind across the rotor (as a measure for veer), to the coefficient wakeSkewCoefficient  that defines the increase in in-plane skew with downstream distance (note that more cases with more or less veer are needed to confirm that first relationship is linear). The linear relationship between downstream distance and skew is built-in into FLORIS, although also a constant skew can be predefined (wakeSkew)
* Adjustments had to be made to CP and CT curve to make a good prediction of power in shear. Note that CCBlade should have some option to correct for veer (but assumes a power-law).
* Velocity shear profile in FLORIS can be predefined as a scaling factor that changes with height, normalized by hub-height (shearProfileZnorm, shearProfileUnorm), or assuming a power-law w
ith coefficient shearCoefficientAlpha
* We use downstream corrections close to the rotor using an arctan profile, k_nearwake could be used to adjust somewhat.



