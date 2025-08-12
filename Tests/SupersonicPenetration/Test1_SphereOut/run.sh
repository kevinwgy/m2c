#!/bin/bash

#SBATCH --job-name=test_run           # Job name
#SBATCH --output=srun_log.out         # Output file
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -t 1:00:00
#SBATCH -p normal_q 
#SBATCH -A m2clab 

# Add any modules you might require. 
module reset
module load PETSc/3.20.3-foss-2023a
module load gomkl/2023b

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

# Run the MPI program mpiProg. srun is Slurm's program for running on assigned resources.
mpirun -n 32 /home/kevinw3/tinkercliffs/zoey7/m2c input.st 2>&1 | tee log.out

exit;
