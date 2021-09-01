#!/bin/bash
#SBATCH --job-name=tarb          # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=15       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=40:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=tanf@slu.edu

module purge
module load matlab/R2019a

matlab -nodisplay -nosplash -r demo
