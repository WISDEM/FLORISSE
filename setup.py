from distutils.core import setup
from distutils.extension import Extension
import numpy as np

try:
    USE_CYTHON = True
    from Cython.Build import cythonize
except Exception:
    USE_CYTHON = False


ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension('_floris', ['src/_floris'+ext, 'src/FLORISmodel.c'],
                        include_dirs=[np.get_include()])]

if USE_CYTHON:
    extensions = cythonize(extensions)

setup(
    name='Floris',
    description='',
    author='',
    author_email='',
    package_dir={'': 'src'},
    py_modules=['floris'],
    license='Apache License, Version 2.0',
    ext_modules=extensions
)