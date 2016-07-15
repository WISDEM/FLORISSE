# EEover setup for windows
#
# FOR WINDOWS, BEFORE running this file with 'python setupWindows.py build_ext --inplace' .-
# - Download GSL Windows DLL and headers for either 32 and 64-bit from https://code.google.com/archive/p/oscats/downloads, For example, gsl-1.15-dev-win64.zip for the 64-bit. 
# - Place contents of zip-file in folder called GSL in the EEOver directory, such that EEOver/GSL/lib, EEOver/GSL/include, EEOver/GSL/bin are there.
#
# Then test by running: python callEE.py
# 
# Pieter Gebraad 2016

from distutils.core import setup, Extension

# define the extension module
ellipseEllipseOverlap = Extension('ellipseEllipseOverlap', sources=['pythonEEover.c', 'solvers.c',  'zsolve_quartic.c', 'Roots3And4.c'], libraries=['gsl', 'gslcblas'],library_dirs=['GSL/lib'], include_dirs=['GSL/include'])

# move GSL dll's to current directory
from shutil import copy
from glob import glob
from os import getcwd

print('Copying over DLLs')
dest_dir = getcwd()
for file in glob(r'GSL\bin\libgsl*.dll'):
    print file                                                                                                                                        
    copy(file, dest_dir)
print('Done')
	
# run the setup
setup(ext_modules=[ellipseEllipseOverlap])