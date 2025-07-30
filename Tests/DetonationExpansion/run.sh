#!/bin/sh

#!/bin/bash
#SBATCH --job-name=test_run           # Job name
#SBATCH --output=srun_log.out         # Output file
#SBATCH --error=srun_log.out          # Error file
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=128         # Number of tasks per node
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --partition=normal_q          # Partition or queue name
#SBATCH --account=m2clab              # Cluster account

# Load necessary modules (e.g., MPI if required)

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

# Run the job
M2C_SIZE=16
M2C_EXE=~/tinkercliffs/zoey7/m2c
M2C_INPUT=input.st

mpiexec -n $M2C_SIZE $M2C_EXE $M2C_INPUT
