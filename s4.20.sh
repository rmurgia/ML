#! /bin/bash

#SBATCH --partition=regular1
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread                                                                                   

#SBATCH --job-name=TAU
#SBATCH --time=12:00:00
#SBATCH --mail-user=rmurgia@sissa.it
#SBATCH --mail-type=ALL

python /home/rmurgia/ML/read_los_tau_ALL_420.py > TAU_LOS_1to3
