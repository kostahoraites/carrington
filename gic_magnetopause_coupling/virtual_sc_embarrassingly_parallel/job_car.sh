#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --job-name=virtual_sc
#SBATCH --partition=long
#SBATCH --exclusive
#SBATCH -M carrington
#SBATCH --nodes=1
#SBATCH -c 64                  # CPU cores per task
#SBATCH -n 1                # number of tasks (4xnodes)
##sbatCH --mem=55G # memory per node

umask 007
ulimit -c unlimited
cd $SLURM_SUBMIT_DIR

t=8

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES
#Vorna has 2 x 8 cores
cores_per_node=32
# Hyperthreading
ht=2
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t


#srun python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.{STARTFILE..ENDFILE}.vlsv -var vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_poynting.x vg_poynting.y vg_poynting.z -c coords_galaxy_15_btrace.txt -re -n 64 > OUTFILE.txt
#srun python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.{STARTFILE..ENDFILE}.vlsv -var vg_e_vol.x vg_e_vol.y vg_e_vol.z vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z proton/vg_precipitationdifferentialflux.6 -c /wrk-vakka/users/horakons/carrington/virtual_sc/txt_files/coords_grid.txt -re -n 64 > OUTFILE.txt
srun python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/bulk1.{STARTFILE..ENDFILE}.vlsv -var vg_j.x vg_j.y vg_j.z vg_j_parallel -c /wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/txt_files/fg_innerboundary_btrace_FHAFGB.txt -re -n 64 > OUTFILE.txt



