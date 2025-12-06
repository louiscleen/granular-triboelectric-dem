#!/bin/bash
# Submission script for Lemaitre4
#SBATCH --job-name=
#SBATCH --no-requeue
#SBATCH --time=24:00:00 # hh:mm:ss
#
#SBATCH --ntasks=200
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=1000 # megabytes
#SBATCH --partition=batch
#
#SBATCH --mail-user=
#SBATCH --mail-type=ALL

module purge
module load OpenMPI

srun ./dem_LM4
