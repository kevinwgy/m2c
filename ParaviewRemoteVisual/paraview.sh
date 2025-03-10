#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
###SBATCH --ntasks-per-node=32
#SBATCH -t 4:00:00
#SBATCH -p normal_q
###SBATCH -p dev_q
#SBATCH -A m2clab 
###SBATCH -A wang_aoe_lab 

# Add any modules you might require. 
module load ParaView/5.11.2-foss-2023a 
cd $SLURM_SUBMIT_DIR
mpirun -n 8 pvserver
#mpirun -n 32 pvserver
exit;
