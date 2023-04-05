#!/bin/sh
#
#SBATCH --job-name=job-7RE8-1OF2
#SBATCH --output=jobfiles/thread-7RE8-1OF2.txt
#SBATCH --error=jobfiles/error-7RE8-1OF2.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=4
module load openmpi/gcc/
bash run_process.sh 7RE8 1OF2 run-dir cartesian
exit
