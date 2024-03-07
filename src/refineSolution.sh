#!/bin/bash
#
#SBATCH -p smp  # partition (queue)
#SBATCH -N 1  # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --cpus-per-task=12
#SBATCH --constraint=FC430
#SBATCH --mem 128G # memory pool for all cores
#SBATCH -t 2-23:00 # time (D-HH:MM)

module load matlab-2022b;
echo "Successfully loaded MATLAB!";


LOADDIR=$1
LOADFILE=$2;
N=$3;

matlab -nodesktop  -nosplash -r "refineSolution $LOADDIR $LOADFILE $N" 

