#! /bin/bash

#SBATCH --partition=regular1
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --hint=nomultithread                                                                                   

#SBATCH --job-name=PF
#SBATCH --time=06:00:00
#SBATCH --mail-user=rmurgia@sissa.it
#SBATCH --mail-type=ALL

python /home/rmurgia/ML/compute_PF.py >> PF_output
