README for FLORIS in OpenMDAOv1.x

Compile locally:
-enter containing directory
-$cc -c adStack.c
-$gfortran -c adBuffer.f

-$f2py -c --opt=-O2 -m _floris floris.f90 adBuffer.o adStack.o
or
-$f2py -c --opt=-O2 -m _floris floriscosine.f90 adBuffer.o adStack.o
- comment out necessary lines in floris_openmdao1.py depending on which file was compiled

Compile on Marylou
-install openmdao
--download the openmdao repo to your computer
--upload the openmdao repo to marylou in the desired directory
--$pip install --user -e .
--check with $nostests .
-compiler FLORIS
--$module unload compiler_gnu*
--$module load compiler_gnu/5.2.0
--$cc -fPIC -c adStack.c
--$gfortran -fPIC -c adBuffer.f

--$f2py -c --opt=-O2 -m _floris floris.f90 adBuffer.o adStack.o
or
--$f2py -c --opt=-O2 -m _floris floriscosine.f90 adBuffer.o adStack.o
-- comment out necessary lines in floris_openmdao1.py depending on which file was compiled

