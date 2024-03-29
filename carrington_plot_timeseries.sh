#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH -J test_figure
#SBATCH -p short
#SBATCH -n 1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=50G
#SBATCH -M carrington



export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
#module load Python/3.7.2-GCCcore-8.2.0




python carrington_plot_timeseries.py -run EGL 



echo Job complete.
