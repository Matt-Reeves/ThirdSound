#!/bin/bash
#
#SBATCH -p smp  # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --constraint=FC430
#SBATCH --mem 5000 # memory pool for all cores
#SBATCH -t 0-23:00 # time (D-HH:MM)
#SBATCH -o gpuPVM.out # STDOUT
#SBATCH -e gpuPVM.err # STDERR
##SBATCH -e slurm.%N.%j.err # STDERR

module load matlab;
echo "Successfully loaded MATLAB!";


H=$1;
srun matlab -nodesktop -nojvm  -nosplash -r "analyse" 