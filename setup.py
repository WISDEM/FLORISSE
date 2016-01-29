#!/usr/bin/env python
# encoding: utf-8

from setuptools import setup
from numpy.distutils.core import setup, Extension

module = Extension('_floris', sources=['src/floris.f90', 'src/adStack.c', 'src/adBuffer.f'],
                   extra_compile_args=['-O2', '-c'])
setup(
    name='FLORIS',
    version='',
    description='floris wake model with cosine factor',
    author='',
    author_email='',
    package_dir={'': 'src'},
    py_modules=['floris'],
    install_requires=['openmdao>=1.4'],
    test_suite='test/tests.py',
    license='Apache License, Version 2.0',
    zip_safe=False
)
