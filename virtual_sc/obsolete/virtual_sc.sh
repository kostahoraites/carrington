#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH -J test_figure
#SBATCH -p long
#SBATCH -n 1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=20G
#SBATCH -M carrington

export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
module load Python/3.7.2-GCCcore-8.2.0


#the old way of generating data. see job_car.txt which runs the updated, mpi version, use that instead. 
##NOTE: vlsvintpol.py won't work unless the coordinates file (e.g. coords.txt) only specifies coordinates that are within the vlasiator domain

#python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.000{0621..1760}.vlsv -var vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z vg_b_vol.x vg_b_vol.y vg_b_vol.z -c coords_galaxy_15_anomaly.txt -re -n 16 | grep -v "Using" > virtual_sc_data_$run.txt


##Galaxy 15 Anomaly. Approximate spacecraft position at [-6.5, 0.8, 0.5] R_E
#run="egi"
run="egl"


# traced field lines

#srun python virtual_sc.py -keogram -run $run -fileprefix timeseries_btrace -filein txt_files/virtual_sc_data_btrace_${run}.txt -poynting -filter -var Spar Sperp Sx Sy Sz
#srun python virtual_sc.py -keogram -norm -run $run -fileprefix timeseries_btrace
#srun python virtual_sc.py -keogram -filter -norm -run $run -fileprefix timeseries_btrace

# heat map on a square grid (defined by running coords_squaregrid.py)
#srun python virtual_sc.py -heatmap -run $run -fileprefix squaregrid_xz_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}.txt -r1 X_RE -r2 Z_RE -poynting -filter
srun python virtual_sc.py -heatmap -run $run -fileprefix squaregrid_xz_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside.txt -r1 X_RE -r2 Z_RE -poynting -filter -norm

#

#srun python virtual_sc.py -keogram -run $run -fileprefix aurora_phi_150 -filein csv_files/aurora_keogram_data_egl_phi_150.csv -xvar theta_GSE_deg -yvar t_sec -xlim 0 60
#srun python virtual_sc.py -keogram -run $run -fileprefix aurora_phi_145_2_ -filein csv_files/aurora_keogram_data_egl_phi_145.csv -xvar theta_GSE_deg -yvar t_sec -xlim 0 60
#srun python virtual_sc.py -keogram -run $run -fileprefix aurora_phi_145 -filein csv_files/aurora_keogram_data_egl_phi_145.csv -r1 theta_GSE_deg -t t_sec



echo Job complete.
