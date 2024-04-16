#!/bin/bash

#### Resource Request: ####
# Tinkercliffs: There are 308 "base compute nodes", each with 128 cores (AMD).
#
# Examples:
#   Request 4 cores (on any number of nodes)
#   #SBATCH --ntasks=4
#   Request exclusive access to all resources on 2 nodes 
#   #SBATCH --nodes=2 
#   #SBATCH --exclusive
#   Request 4 cores (on any number of nodes)
#   #SBATCH --ntasks=4
#   Request 2 nodes with 12 tasks running on each
#   #SBATCH --nodes=2
#   #SBATCH --ntasks-per-node=12
#   Request 12 tasks with 20GB memory per core
#   #SBATCH --ntasks=12 
#   #SBATCH --mem-per-cpu=20G
#   Request 5 nodes and spread 50 tasks evenly across them
#   #SBATCH --nodes=5
#   #SBATCH --ntasks=50
#   #SBATCH --spread-job
#
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128

#### Walltime ####
# Set the walltime, which is the maximum time your job can run in HH:MM:SS
# Note that if your job exceeds the walltime estimated during submission, the scheduler
# will kill it. So it is important to be conservative (i.e., to err on the high side)
# with the walltime that you include in your submission script. 
#SBATCH -t 2:00:00

#### Queue ####
# Queue name. Cascades has five queues:
#   normal_q        for production jobs on all Broadwell nodes
#   largemem_q      for jobs on the two 3TB, 60-core Ivy Bridge servers
#   dev_q           for development/debugging jobs. These jobs must be short but can be large.
#   v100_normal_q   for production jobs on Skylake/V100 nodes
#   v100_dev_q      for development/debugging jobs on Skylake/V100 nodes
#SBATCH -p normal_q 

#### Account ####
# This determines which allocation this job's CPU hours are billed to.
# Replace "youraccount" below with the name of your allocation account.
# If you are a student, you will need to get this from your advisor.
# For more on allocations, go here: http://www.arc.vt.edu/allocations
#SBATCH -A m2clab 

# Add any modules you might require. 
module load iimpi/2019b

# Change to the directory from which the job was submitted
cd $SLURM_SUBMIT_DIR

# Below here enter the commands to start your job. A few examples are provided below.
# Some useful variables set by the job:
#  $SLURM_SUBMIT_DIR   Directory from which the job was submitted
#  $PBS_NODEFILE       File containing list of cores available to the job
#  $PBS_GPUFILE        File containing list of GPUs available to the job
#  $SLURM_JOBID        Job ID (e.g., 107619.master.cluster)
#  $SLURM_NTASKS       Number of cores allocated to the job
# You can run the following (inside a job) to see what environment variables are available:
#  env | grep SLURM
#
# Some useful storage locations (see ARC's Storage documentation for details): 
#  $HOME     Home directory. Use for permanent files.
#  $WORK     Work directory. Use for fast I/O.
#  $TMPFS    File system set up in memory for this job. Use for very fast, small I/O
#  $TMPDIR   Local disk (hard drive) space set up for this job

# Run the MPI program mpiProg. srun is Slurm's program for running on assigned resources.
mpirun -n 256 ../../../m2c input.st 2>&1 | tee log.out
# The -n flag tells srun how many processes to use. $SLURM_NTASKS is an environment variable
# that holds the number of tasks you requested. So if you selected --nodes=2 and --ntasks-per-node=12
# above, $SLURM_NTASKS will hold 24.
# srun -n $SLURM_NTASKS ./mpiProg

exit;
