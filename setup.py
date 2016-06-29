#!/usr/bin/env python
# encoding: utf-8

from numpy.distutils.core import setup, Extension

module1 = Extension('_floris', sources=['src/florisse/floris.f90', 'src/florisse/adStack.c', 'src/florisse/adBuffer.f'],
                   extra_compile_args=['-O2', '-c'])

module2 = Extension('_florisDiscontinuous', sources=['src/florisse/florisDiscontinuous.f90', 'src/florisse/adStack.c', 'src/florisse/adBuffer.f'],
                   extra_compile_args=['-O2', '-c'])

setup(
    name='FLORISSE',
    version='0.0.0',
    description='differentiable floris wake model with cosine factor',
    install_requires=['openmdao>=1.5','akima>=1.0.0'],
    package_dir={'': 'src'},
    ext_modules=[module1, module2],
    dependency_links=['https://github.com/andrewning/akima/tarball/master#egg=akima'],
    packages=['florisse'],
    license='Apache License, Version 2.0',
)
