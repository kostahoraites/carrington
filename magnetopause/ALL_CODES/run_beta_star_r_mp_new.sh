#!/bin/bash -l
###carrington:

#SBATCH -J beta_star_r_mp_new
#SBATCH -t 06:00:00
#SBATCH -M carrington
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32       # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=160G # memory per node 20G per task?
##SBATCH --no-requeue




##!/bin/bash

#######
##SLURM
#######

##SBATCH --time=1-00:00:00
##SBATCH --exclusive

##carrington:

##SBATCH -M carrington
##SBATCH --partition=long
##SBATCH --ntasks=1
##SBATCH --nodes=1
##SBATCH --cpus-per-task=32        # cpu-cores per task (>1 if multi-threaded tasks)
##SBATCH --mem=160G                # memory per node 160G/32 per task?


export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
#module load Python/3.7.2-GCCcore-8.2.0


srun python beta_star_r_mp_new.py -nproc 32

echo running...
echo Job complete.
