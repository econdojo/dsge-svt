#!/bin/bash
#SBATCH --job-name=sw07_bmm   # create a short name for your job
#SBATCH --partition=defq         # choose partition
#SBATCH --nodes=5               # node count
#SBATCH --ntasks=5              # total number of tasks across all nodes
#SBATCH --cpus-per-task=12       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --time=50:00:00       # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=fei.tan@slu.edu

module purge
module load slurm
module load MATLAB

# Execute jobs in parallel

srun -N 1 -n 1 matlab -nodisplay -nosplash -r "test_bmm(1)" &
srun -N 1 -n 1 matlab -nodisplay -nosplash -r "test_bmm(2)" &
wait