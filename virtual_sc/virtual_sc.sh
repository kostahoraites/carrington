#!/bin/bash

#######
##SLURM
#######

#SBATCH --time=1-00:00:00
#SBATCH --exclusive

##SBATCH -o proc_c64_final.out
##SBATCH --job-name=proc_c64_final

###carrington:

#SBATCH -M carrington
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=160G # memory per node 20G per task?



export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
module load Python/3.7.2-GCCcore-8.2.0


#the old way of generating data. see job_car.txt which runs the updated, mpi version, use that instead. 
##NOTE: vlsvintpol.py won't work unless the coordinates file (e.g. coords.txt) only specifies coordinates that are within the vlasiator domain

#python /proj/horakons/analysator/tools/vlsvintpol.py -i /wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.000{0621..1760}.vlsv -var vg_b_vol.magnitude vg_rho proton/vg_pdyn vg_v.x vg_v.y vg_v.z vg_b_vol.x vg_b_vol.y vg_b_vol.z -c coords_galaxy_15_anomaly.txt -re -n 16 | grep -v "Using" > virtual_sc_data_$run.txt


##Galaxy 15 Anomaly. Approximate spacecraft position at [-6.5, 0.8, 0.5] R_E
#run="egi"
run="egl"


# straight line x=-15 to x=15
#srun python virtual_sc.py -keogram -run $run -fileprefix x_line_x-15_15_y0_z0 -filein txt_files/virtual_sc_data_line_${run}_x-15_15_y0_z0.txt -norm -var vg_b_vol.magnitude vg_rho
#srun python virtual_sc.py -timeseries -run $run -fileprefix x_line_x-15_15_y6_z2 -filein txt_files/virtual_sc_data_line_${run}_x-15_15_y6_z2_${run}.txt -norm -var vg_b_vol.magnitude vg_rho
#srun python virtual_sc.py -timeseries -run $run -fileprefix x_line_x-15_15_y0_z0 -filein txt_files/virtual_sc_data_line_${run}_x-15_15_y0_z0.txt -var vg_b_vol.magnitude B.x B.y B.z


# single point
#srun python virtual_sc.py -timeseries -wavelet -run $run -fileprefix timeseries_1pt -filein txt_files/virtual_sc_data_dayside_magnetopause_1pt_${run}.txt -poynting -ta -var Spar vg_rho vg_b_vol.magnitude

# traced field lines

#srun python virtual_sc.py -timeseries -keogram -run $run -fileprefix timeseries_btrace -filein txt_files/virtual_sc_data_btrace_${run}.txt -poynting -filter -var Spar vg_rho vg_b_vol.magnitude
srun python virtual_sc.py -timeseries -keogram -ta -run $run -fileprefix timeseries_dayside_magnetopause_btrace -filein txt_files/virtual_sc_data_dayside_magnetopause_btrace_${run}.txt -poynting -filter -var Spar B.z vg_rho -ylim 800 1200




#aurora (carrington.py data)

#srun python virtual_sc.py -keogram -run $run -fileprefix aurora_phi_150 -filein csv_files/aurora_keogram_data_egl_phi_150.csv -xvar theta_GSE_deg -yvar t_sec -xlim 0 60
#srun python virtual_sc.py -keogram -run $run -fileprefix aurora_phi_145_2_ -filein csv_files/aurora_keogram_data_egl_phi_145.csv -xvar theta_GSE_deg -yvar t_sec -xlim 0 60
#srun python virtual_sc.py -keogram -run $run -fileprefix aurora_phi_145 -filein csv_files/aurora_keogram_data_egl_phi_145.csv -r1 theta_GSE_deg -t t_sec


# heat map on a square grid (defined by running coords_squaregrid.py)
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xz_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside.txt -r1 X_RE -r2 Z_RE -dt -abs -log -cmap nipy_spectral -var vg_b_vol.magnitude -vmin 1e-10 -vmax 5e-9
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xy_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xy.txt -r1 X_RE -r2 Y_RE -dt -abs -log -cmap nipy_spectral -var vg_b_vol.magnitude -vmin 1e-10 -vmax 5e-9
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xy_z0_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xy_z0_hires.txt -r1 X_RE -r2 Y_RE -dt -abs -cmap nipy_spectral -var vg_b_vol.magnitude -vmin 0 -vmax 2.5e-9
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xz_y0_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xz_y0_hires.txt -r1 X_RE -r2 Z_RE -dt -abs -cmap nipy_spectral -var vg_b_vol.magnitude -vmin 0 -vmax 2.5e-9
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xy_z0_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xy_z0_hires.txt -r1 X_RE -r2 Y_RE -dt -abs -cmap nipy_spectral -var vg_rho -vmin 0 -vmax 4e5
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xz_y0_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xz_y0_hires.txt -r1 X_RE -r2 Z_RE -dt -abs -cmap nipy_spectral -var vg_rho -vmin 0 -vmax 4e5
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xz_y0_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xz_y0_hires.txt -r1 X_RE -r2 Y_RE -dt -abs -log -cmap nipy_spectral -var proton/vg_precipitationdifferentialflux.6 -vmin 1e-6 vmax 1e2
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xy_z5_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xy_z5.txt -r1 X_RE -r2 Y_RE -dt -abs -log -cmap nipy_spectral -var vg_b_vol.magnitude -vmin 1e-10 -vmax 5e-9
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xy_z0_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xy_z0_hires.txt -r1 X_RE -r2 Y_RE -dt -abs -log -cmap plasma -var proton/vg_precipitationdifferentialflux.6 -vmin 1e-6 -vmax 1e2
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xz_y0_dayside -filein txt_files/virtual_sc_data_squaregrid_${run}_dayside_xz_y0_hires.txt -r1 X_RE -r2 Z_RE -dt -abs -log -cmap plasma -var proton/vg_precipitationdifferentialflux.6 -vmin 1e-6 -vmax 1e2
#srun python virtual_sc.py -nproc 16 -heatmap -run $run -fileprefix squaregrid_xz_nightside -filein txt_files/virtual_sc_data_squaregrid_${run}_nightside.txt -r1 X_RE -r2 Z_RE -poynting -filter


echo Job complete.
