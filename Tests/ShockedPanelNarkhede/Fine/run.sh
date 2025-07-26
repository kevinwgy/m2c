#!/bin/sh

#!/bin/bash
#SBATCH --job-name=giordano           # Job name
#SBATCH --output=srun_log.out         # Output file
#SBATCH --error=srun_log.out          # Error file
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=128         # Number of tasks per node
#SBATCH --time=12:00:00               # Time limit hrs:min:sec
#SBATCH --partition=normal_q          # Partition or queue name
#SBATCH --account=m2clab              # Cluster account

export UCX_LOG_LEVEL=error

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

# Run the job
M2C_SIZE=120
M2C_EXE=~/tinkercliffs/m2c/m2c
M2C_INPUT=input.st

AEROS_SIZE=1
AEROS_EXE=~/tinkercliffs/FEMWorkingFoam/bin/aeros
AEROS_INPUT=fem.in

mpiexec -n $M2C_SIZE $M2C_EXE $M2C_INPUT :\
 -n $AEROS_SIZE $AEROS_EXE $AEROS_INPUT 2>&1 | tee log.out

