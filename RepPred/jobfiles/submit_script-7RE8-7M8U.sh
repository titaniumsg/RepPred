#!/bin/sh
#
#SBATCH --job-name=job-7RE8-7M8U
#SBATCH --output=jobfiles/thread-7RE8-7M8U.txt
#SBATCH --error=jobfiles/error-7RE8-7M8U.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=4
module load openmpi/gcc/
bash run_process.sh 7RE8 7M8U run-dir cartesian
exit
