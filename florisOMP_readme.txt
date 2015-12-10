NOTES ON COMPILING AND OTHER THINGS FOR RUNNING PARALLEL FLORIS ON MARYLOU



Before compiling
$module purge
$module load python/2.7.10 petsc/3.6.3
$module load petsc/3.6.3



Compiling
$ gcc -c -fPIC adStack.c
$ gfortran -c -fPIC adBuffer.f 
$ f2py -c --f90flags='-fopenmp' -lgomp --opt=-O2 -m _floris florisOMP.f90 adBuffer.o adStack.o 



Running on Marylou
submission_script.sh:

#!/bin/bash

#SBATCH --time=01:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=8   # number of nodes
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH --qos=standby
#SBATCH -Cavx # option required for

# Compatibility variables for PBS.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# for some reason export OMP_NUM_THREADS is not updating so this number is irrelevent in this case
# instead use -genv OMP_NUM_THREADS $SLURM_CPUS_PER_TASK in the mpi call

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# normally I think this should be -x OMP_NUM_THREADS, but this setup requires -genv (maybe due to hydra?)
mpirun -np 2 -genv OMP_NUM_THREADS $SLURM_CPUS_PER_TASK python exampleOptimizationAEP.py



Running on a workstation
$ export OMP_NUM_THREADS = 8
$ mpirun -np 2 python exampleOptimizationAEP.py
