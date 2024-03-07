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

N=$1;
D=$2;
A=$3;
W=$4;
DA=$5;
AM=$6
matlab -nodesktop  -nosplash -r "doglegSolve $N $D $A $W $DA $AM" 

