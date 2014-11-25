FLORISSE
========


Created by Pieter Gebraad and Paul Fleming.
Copyright (c) NREL. All rights reserved.

-This archive contains the c implementation of FLORIS, with python wrapper through CYTHON

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