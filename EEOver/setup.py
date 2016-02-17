from distutils.core import setup, Extension

# define the extension module
ellipseEllipseOverlap = Extension('ellipseEllipseOverlap', sources=['pythonEEover.c', 'solvers.c',  'zsolve_quartic.c', 'Roots3And4.c'], libraries=['gsl', 'gslcblas'])


# run the setup
setup(ext_modules=[ellipseEllipseOverlap])