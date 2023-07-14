#!/bin/bash -l
###carrington:

#SBATCH -J validation_paper_plots
#SBATCH -t 6:00:00
#SBATCH -M carrington
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=160G # memory per node 20G per task?
##SBATCH --no-requeue


umask 007
ulimit -c unlimited
cd $SLURM_SUBMIT_DIR

t=8

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES

## Carrington: has 2x16 cores (i think)
#cores_per_node=32

## Vorna: has 2 x 8 cores
cores_per_node=16

# Hyperthreading
ht=2
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t

#--------------------------------------------------------------------



export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
#module load Python/3.7.2-GCCcore-8.2.0

#time python validation_paper_plots.py -nproc 64
time python validation_paper_plots.py -nproc 1

echo Job complete.
