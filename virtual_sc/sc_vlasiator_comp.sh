#!/bin/bash -l
#SBATCH -t 02:00:00
#SBATCH -J test_figure
#SBATCH -p short
#SBATCH -n 1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=20G
#SBATCH -M ukko

export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
#module load Python/3.7.2-GCCcore-8.2.0

##EGL (modify virtual_sc.py accordingly):
run="egl"
#mpirun python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.000{0621..1760}.vlsv -var vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude proton/vg_v.magnitude proton/vg_v.x proton/vg_v.y proton/vg_v.z proton/vg_rho proton/vg_temperature proton/vg_pdyn vg_e_vol.magnitude proton/vg_beta -c txt_files/coords_sc_vlasiator_comp.txt -re -n 6 | grep -v "Using" > txt_files/virtual_sc_data_vlasiator_comp_$run.txt
#python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.000{0621..1760}.vlsv -var vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude proton/vg_v.magnitude proton/vg_v.x proton/vg_v.y proton/vg_v.z proton/vg_rho proton/vg_temperature proton/vg_pdyn vg_e_vol.magnitude proton/vg_beta -c txt_files/coords_sc_vlasiator_comp.txt -re -n 6 | grep -v "Using" > txt_files/virtual_sc_data_vlasiator_comp_$run.txt
  # had to run a second time to get the last few files?
#python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.000{1743..1760}.vlsv -var vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude proton/vg_v.magnitude proton/vg_v.x proton/vg_v.y proton/vg_v.z proton/vg_rho proton/vg_temperature proton/vg_pdyn vg_e_vol.magnitude proton/vg_beta -c txt_files/coords_sc_vlasiator_comp.txt -re -n 6 | grep -v "Using" > txt_files/virtual_sc_data_vlasiator_comp_$run_part2.txt
#EGP (modify virtual_sc.py accordingly):
#run="egp"
#python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk/group/spacephysics/vlasiator/3D/EGP/bulk1/bulk1.000{0269..0506}.vlsv -var vg_b_vol.x vg_b_vol.y vg_b_vol.z vg_b_vol.magnitude proton/vg_v.magnitude proton/vg_v.x proton/vg_v.y proton/vg_v.z proton/vg_rho proton/vg_temperature proton/vg_pdyn vg_e_vol.magnitude proton/vg_beta -c coords_galaxy_15_anomaly.txt -re -n 8 | grep -v "Using" > virtual_sc_data_$run.txt




# notes on the flags, from vlsvintpol.py documentation
#parser.add_argument('-var', nargs='*', help="a list of variable.operator's to output, e.g., v.magnitude rho B.x " )
#parser.add_argument('-i', nargs='*', help="a list of vlsv files")
#parser.add_argument('-c', help="A file with coordinates (can also be give from stdin)")
#parser.add_argument('-re', action='store_true', help="Coordinates in RE, in meters by default")
#parser.add_argument('-n', help="Number of processes to use, default 1")

python sc_vlasiator_comp.py

echo Job complete.
