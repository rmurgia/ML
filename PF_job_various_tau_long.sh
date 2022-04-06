#! /bin/bash

#SBATCH --partition=long2
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread                                                                                   

#SBATCH --job-name=PF_taus
#SBATCH --time=48:00:00
#SBATCH --mail-user=rmurgia@sissa.it
#SBATCH --mail-type=ALL

python /home/rmurgia/ML/compute_PF_various_tau_2021.py >> PF_output_taus_NEW
