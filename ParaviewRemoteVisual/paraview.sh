#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -t 3:00:00
###SBATCH -p normal_q
#SBATCH -p dev_q
#SBATCH -A wang_aoe_lab 

# Add any modules you might require. 
module load ParaView/5.9.1-foss-2021a-mpi
cd $SLURM_SUBMIT_DIR
mpirun -n 32 pvserver
exit;
