#!/bin/bash -l
#SBATCH -t 02:00:00
#SBATCH -J test_figure
#SBATCH -p short
#SBATCH -n 1
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=250G
#SBATCH -M carrington



export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
module load Python/3.7.2-GCCcore-8.2.0


# key for -var flag
# 1. open vs. closed field boundary
# 2. proton differential energy flux
# 3. Field aligned currents (FACS)
# 4. Magnetopause position
# 5. dB/dt. Note: this needs to be evaluated in the ionosphere (once fully implemented by Urs)
# 6. other plots


## arguments
##    '-run', default='EGL', help="the Vlasiator run, in all caps " )
##    '-var', default = ['1','2','3','4','5','6'], nargs='*', help="a list of plot identifiers (numbers), set to which plots you want to make" )
##    '-nphi', default=360, help="the number of longitudes in the lon-lat grid" )
##    '-phimin', default = -180, help="minimum longitude, (range -180 to 180)" ) # note phi = (pi/180) * (phi_input-180)
##    '-phimax', default = 180, help="maximum longitude, (range -180 to 180)" )
##    '-nlat', default=180, help="the number of latitudes in the lon-lat grid" )
##    '-latmin', default = -90, help="minimum longitude, (range -90 to 90)" ) # note lat = -lat_input * (pi/180) 
##    '-latmax', default = 90, help="maximum longitude, (range -90 to 90)" )
##    '-startstop', nargs='*', help="2-element list, start and stop index (divided by deltanframes)" )
##    '-save', action='store_true', help="set this flag to save the data into a .csv format")  #default: save=False when -save flag not set
##    '-savefile', default='/wrk/users/horakons/carrington/data/test_data.csv', help="filename for the .csv data")
##   '-append', action='store_true', help="set this flag to save the data into a .csv format")  #default: save=False when -save flag not set





#python carrington.py 
#python carrington.py -run EGI -var 6
#python carrington.py -startstop 640 641 -run EGL -var 3 -nphi 1 -phimin 150 -phimax 150 -append -save -savefile $HOME_DIR/carrington/data/aurora_keogram_data_egl_phi_150.csv

run=EGL
start=81
stop=82
deltanframes=20  # loop from start*deltanframes to stop*deltanframes, step size deltanframes
python carrington.py -startstop $start $stop -deltanframes 20 -run EGL -var 1 2 3 4 -save -savefile "$HOME_DIR/carrington/data/$run/aurora_plot_data/aurora_plot_data_${run}_${start}_${end}.csv"

#python carrington.py -run EGI -var 1 2 3 4 6



echo Job complete.
