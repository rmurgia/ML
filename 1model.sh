#! /bin/bash

#SBATCH --partition=regular1
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread                                                                                   

#SBATCH --job-name=TAU_lcdm
#SBATCH --time=12:00:00
#SBATCH --mail-user=rmurgia@sissa.it
#SBATCH --mail-type=ALL
#SBATCH --mem=0

python /home/rmurgia/ML/read_los_tau_ALL_2021.py > TAU_lcdm
