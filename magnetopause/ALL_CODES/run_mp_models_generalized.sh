#!/bin/bash -l
#SBATCH -t 08:00:00
#SBATCH -J test_figure
#SBATCH -p short
#SBATCH -n 1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=20G
#SBATCH -M carrington

#module purge
export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
#module load matplotlib/3.2.1-foss-2020a-Python-3.8.2

python mp_models_generalized.py
echo Job complete.
