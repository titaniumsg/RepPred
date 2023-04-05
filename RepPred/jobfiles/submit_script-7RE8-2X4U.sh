#!/bin/sh
#
#SBATCH --job-name=job-7RE8-2X4U
#SBATCH --output=jobfiles/thread-7RE8-2X4U.txt
#SBATCH --error=jobfiles/error-7RE8-2X4U.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=4
module load openmpi/gcc/
bash run_process.sh 7RE8 2X4U run-dir cartesian
exit
