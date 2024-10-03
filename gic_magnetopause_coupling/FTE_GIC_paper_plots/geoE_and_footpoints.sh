#!/bin/bash -l
#SBATCH -t 01:00:00
#SBATCH -J geoE_and_footpoints
#SBATCH -p short
#SBATCH -n 1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=20G
#SBATCH -M carrington

#export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
#module load Python/3.7.2-GCCcore-8.2.0

module purge                 ## Purge modules for a clean start
module load Python/3.7.2-GCCcore-8.2.0
module load matplotlib
export PTNOLATEX=1 PTBACKEND=Agg PTNONINTERACTIVE=1 PYTHONPATH=$PYTHONPATH:/proj/ykempf/analysator/:/proj/ykempf/.local PYTHONUSERBASE=/proj/ykempf/.local

python geoE_and_footpoints.py

echo Job complete.
