#!/bin/bash -l
###carrington:

#SBATCH -J ocb
#SBATCH -t 01:00:00
#SBATCH -M vorna
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=10G # memory per node 20G per task?
##SBATCH --no-requeue





export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
module load Python/3.7.2-GCCcore-8.2.0

time python plot_B_r.py

echo Job complete.


