#!bin/bash
#SBATCH -J myMPI
#SBATCH -o myMPI.o%j
#SBATCH -n 2
#SBATCH -p normal
#SBATCH -t 01:00:00
#SBATCH -A hugo3272
set -x
ibrun ./mult_mpi 0 0 10 
