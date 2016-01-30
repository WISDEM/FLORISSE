README for FLORIS in OpenMDAOv1.5+

## Installation instructions MAC
- system requirements  
    gfortran  
    gcc  
    python 2.7.x  
    numpy  
    openmdao >= v1.5  
- from top repo directory run  
    $ python setup.py install --record installedfiles.txt  
- uninstall with  
    $ cat installedfiles.txt | xargs rm -r  

## Alternative Installation instructions MAC  
- system requirements  
    gfortran  
    gcc  
    python 2.7.x  
    numpy  
    openmdao >= v1.5  
- put all files in desired directory  
- run the following commands from src/florisse:  
    $ gfortran -c adBuffer.f  
    $ gcc -c adStack.c  
    $ f2py -c --opt=-O2 -m _floris floris.f90 adBuffer.o adStack.o  
    
    
## Installation instructions Windows  
- system requirements  
    gfortran  
    gcc  
    mingw  
    python 2.7.x  
    numpy  
    openmdao >= v01.5  
- put all files in desired directory  
- run the following commands from src\florisse:  
    $ gfortran -c adBuffer.f  
    $ gcc -c adStack.c  
    $ python \your\path\to\f2py.py -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _floris floris.f90 adBuffer.o adStack.o  
        (most likely your path is C:\python27\Scripts\f2py.py)  
- if you get an error in the line "as=b['args']" try to update numpy 
    ($ pip install numpy --upgrade)  
- run the example using  
    $ python test\example_call.py  
        

## Installation instructions Marylou  
- module dependencies ($ module load <module name>)(purge all modules first)  
    petsc/3.6.3  
    python2/7  
    compiler_gnu/4.9.2  
    mpi/openmpi-1.8.4_gnu-4.9.2  
- python dependencies ($ pip install --user <package name>)  
    openmdao >= v1.5 (use a clone of the repo and {$ pip install --user -e .} from top level of 
              acquired repo)  
    mpi4py  
    petsc4py      
- compiler FLORIS (clone with ssh on Marylou)  
    $ cd src/florisse  
    $ gcc -fPIC -c adStack.c  
    $ gfortran -fPIC -c adBuffer.f  
    $ f2py -c --opt=-O2 -m _floris floris.f90 adBuffer.o adStack.o  
- test installation (from within src/florisse)  
    $ python test/tests.py  
    $ python test/exampleCall.py  
    $ python test/exampleOptimizationAEP.py 2  
    $ mpirun -np 4 python test/exampleOptimizationAEP.py 2  
