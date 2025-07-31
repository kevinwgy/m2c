#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH -t 00:20:00
#SBATCH -p normal_q 
#SBATCH -A m2clab

# Add any modules you might require. 
module reset
module load PETSc/3.20.3-foss-2023a
module load gomkl/2023b

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

# Run the MPI program mpiProg. srun is Slurm's program for running on assigned resources.
### Solvers and parameters
M2C_SIZE=64
M2C_EXE=/home/kevinw3/tinkercliffs/zoey7/m2c
M2C_INPUT=input.st
AEROS_SIZE=16
AEROS_EXE=/home/kevinw3/tinkercliffs/FEMWorking/bin/aeros 
AEROS_INPUT=fem.in

### Run
mpiexec -n $M2C_SIZE $M2C_EXE $M2C_INPUT : -n $AEROS_SIZE $AEROS_EXE $AEROS_INPUT 2>&1 | tee log.out

exit;
