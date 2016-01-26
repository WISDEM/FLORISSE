README for FLORIS in OpenMDAOv1.x

Installation instructions MAC
-----------------------------
- system requirements
    gfortran
    gcc
    python 2.7.x
    numpy
    openmdao v1.x
- put all files in desired directory
- run the following commands:
    $ gfortran -c adBuffer.f
    $ gcc -c adStack.c
    $ f2py -c --opt=-O2 -m _floris floris.f90 adBuffer.o adStack.o
    
    
Installation instructions Windows
---------------------------------
- system requirements
    gfortran
    gcc
    mingw
    python 2.7.x
    numpy
    openmdao v01.x
- put all files in desired directory
- run the following commands:
    $ gfortran -c adBuffer.f
    $ gcc -c adStack.c
    $ python \your\path\to\f2py.py -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _floris floris.f90 adBuffer.o adStack.o
        (most likely your path is C:\python27\Scripts\f2py.py)
- if you get an error in the line "as=b['args']" try to update numpy 
    ($ pip install numpy --upgrade)
- run the example using
    $ python example_call.py
        

Installation instructions Marylou
---------------------------------
- module dependencies ($ module load <module name>)(purge all modules first)
    petsc/3.6.3
    python2/7
    compiler_gnu/4.9.2
    mpi/openmpi-1.8.4_gnu-4.9.2
- python dependencies ($ pip install --user <package name>)
    openmdao (use a clone of the repo and {$ pip install --user -e .} from top level of 
              acquired repo)
    mpi4py
    petsc4py    
- compiler FLORIS (clone with ssh on Marylou)
    $ gcc -fPIC -c adStack.c
    $ gfortran -fPIC -c adBuffer.f
    $ f2py -c --opt=-O2 -m _floris floris.f90 adBuffer.o adStack.o
- test installation (from within directory containing FLORIS)
    $ python tests.py
    $ python exampleCall.py
    $ python exampleOptimizationAEP.py 2
    $ mpirun -np 4 python exampleOptimizationAEP.py 2
