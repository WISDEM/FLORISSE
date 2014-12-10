FLORISSE
========


Created by Pieter Gebraad and Paul Fleming.
Copyright (c) NREL. All rights reserved.

-This archive contains the c implementation of FLORIS, with python wrapper through CYTHON

REQUIRED PYTHON LIBRARIES:
  -Cython

-For summary of the FLORIS model refer to:
P. M. O. Gebraad, F. W. Teeuwisse, J.-W. van Wingerden, P. A. Fleming, S. D. Ruben, J. R. Marden, and L. Pao, “A Data-Driven Model for Wind Plant Power Optimization by Yaw Control,” in Proceedings of the American Control Conference, 2014, pp. 3128–3134.

-For full details refer to:
Data-driven wind plant control, PhD Thesis, Pieter Gebraad 2014 [see: Chapter 4]
http://dx.doi.org/10.4233/uuid:5c37b2d7-c2da-4457-bff9-f6fd27fe8767

FILES:
setup.py: Setup function to build the python-wrapped FLORIS function
  Usage: "python setup.py build_ext --inplace"

example.py: Example call to python-wrapped from python

FLORISmodel.[c h]: c-implementation of FLORIS

"""
Note: (Pieter)
To set up FLORISSE standalone, set up python with numpy, cython and a compiler. On Windows 7, I used the Anaconda python package (32-bit version), which has the Cython and MinGW compiler included, and I used these instructions to set up the compler (http://docs.cython.org/src/tutorial/appendix.html, with the Anaconda installation path instead of standard Python path) and then the above installation line within the Anaconda Command Prompt. 64-bit has issues because MinGW does not support it yet, and Visual Studio / SDK does not work well.
