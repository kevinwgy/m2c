#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH -t 24:00:00
#SBATCH -p normal_q 
#SBATCH -A m2clab 

# Add any modules you might require. 
module load foss/2023b

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

# Run the MPI program
node_list=(`scontrol show hostnames`)
mpiexec --bind-to none --host $node_list:96 -n 96 /home/kevinw3/owl/zoey7/m2c input.st 2>&1 | tee log.out

exit;
