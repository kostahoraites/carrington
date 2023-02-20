#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --job-name=vlsvintpol_job_car
#SBATCH --partition=long
#SBATCH --exclusive
#SBATCH -M carrington
#SBATCH -c 1                  # CPU cores per task
### need to set the proper number of nodes and tasks.
##SBATCH --nodes=2
##SBATCH -n 128             # number of tasks (64xnodes)
##SBATCH --nodes=5
##SBATCH -n 320             # number of tasks (64xnodes)
#SBATCH --nodes=8
#SBATCH -n 512             # number of tasks (64xnodes)
##SBATCH --nodes=9
##SBATCH -n 576             # number of tasks (64xnodes)
##SBATCH --nodes=10
##SBATCH -n 640             # number of tasks (64xnodes)
##sbatCH --mem=55G # memory per node

umask 007
ulimit -c unlimited
cd $SLURM_SUBMIT_DIR
export PTNONINTERACTIVE=1 

module purge
module load mpi4py/3.0.0-intel-2018b-Python-3.7.0
module load matplotlib/3.2.1-foss-2020a-Python-3.8.2
module load OpenMPI/4.0.3-GCC-9.3.0


t=1

#--------------------------------------------------------------------
#---------------------DO NOT TOUCH-----------------------------------
nodes=$SLURM_NNODES
#Vorna has 2 x 8 cores, what about Carrington?
cores_per_node=32
# Hyperthreading
ht=2
#Change PBS parameters above + the ones here
total_units=$(echo $nodes $cores_per_node $ht | gawk '{print $1*$2*$3}')
units_per_node=$(echo $cores_per_node $ht | gawk '{print $1*$2}')
tasks=$(echo $total_units $t  | gawk '{print $1/$2}')
tasks_per_node=$(echo $units_per_node $t  | gawk '{print $1/$2}')
export OMP_NUM_THREADS=$t


# straight line x=-15 to x=15 (as per andreeova 2011)

#python coords_line.py

#egl
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_line_egl.txt -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -var vg_b_vol.magnitude vg_rho -c txt_files/coords_line.txt -re 
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_line_egl_x-15_15_y0_z0.txt -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -var vg_b_vol.magnitude vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_rho -c txt_files/coords_line.txt -re 

#egi
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_line_egi_x-15_15_y0_z0.txt -i /wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/bulk1.*.vlsv -var vg_b_vol.magnitude vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_rho -c txt_files/coords_line.txt -re 



#mpirun  python vlsvintpol_mpi.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.{0000651..0001750}.vlsv -var vg_b_vol.magnitude vg_rho  -c coords_galaxy_15_btrace.txt -re 

#mpirun  python vlsvintpol_mpi.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -var vg_b_vol.magnitude vg_rho  -c coords_galaxy_15_btrace.txt -re 


#mpirun  python vlsvintpol_mpi.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -var vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z vg_b_vol.x vg_b_vol.y vg_b_vol.z   -c coords_galaxy_15_btrace.txt -re 



#### interpolate along a traced field line:

#python coords_galaxy_15_btrace.py


#example variables:            vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z vg_b_vol.x vg_b_vol.y vg_b_vol.z -c coords_galaxy_15_btrace.txt -re


#egl
mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_dayside_magnetopause_1pt_egl.txt -i /wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -var vg_e_vol.x vg_e_vol.y vg_e_vol.z vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude vg_rho -c txt_files/coords_dayside_magnetopause_1pt.txt -re
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_btrace_egl.txt -i /wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -var vg_e_vol.x vg_e_vol.y vg_e_vol.z vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude vg_rho -c txt_files/coords_galaxy_15_btrace_egl.txt -re
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_dayside_magnetopause_btrace_egl.txt -i /wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -var vg_e_vol.x vg_e_vol.y vg_e_vol.z vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude vg_rho -c txt_files/coords_dayside_magnetopause_btrace_egl.txt -re

#egi
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_btrace_egi.txt -i /wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/bulk1.*.vlsv -var vg_e_vol.x vg_e_vol.y vg_e_vol.z vg_b_vol.x vg_b_vol.y vg_b_vol.z -c txt_files/coords_galaxy_15_btrace_egi.txt -re


### interpolate data on square grid:


#python coords_grid.py -nx 7 -ny 10 -xlim -15 20 -ylim -25 25 -zlim 0 0

#egl
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_squaregrid_egl_dayside.txt -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -c txt_files/coords_grid.txt -re -var vg_e_vol.x vg_e_vol.y vg_e_vol.z vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_squaregrid_egl_dayside_xy.txt -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.*.vlsv -c txt_files/coords_grid.txt -re -var vg_e_vol.x vg_e_vol.y vg_e_vol.z vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z proton/vg_precipitationdifferentialflux.6

##egi
#mpirun python vlsvintpol_mpi.py -outfile txt_files/virtual_sc_data_squaregrid_egi_dayside.txt -i /wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/bulk1.*.vlsv -c txt_files/coords_grid.txt -re -var vg_e_vol.x vg_e_vol.y vg_e_vol.z vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z
