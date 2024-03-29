#!/bin/bash

#sbatch carrington_multithread.sh nproc start stop deltanframes

#example: sbatch carrington_multithread.sh 16 901 916


#######
##SLURM
#######

#SBATCH --time=06:00:00
#SBATCH --exclusive

##SBATCH -o proc_c64_final.out
##SBATCH --job-name=proc_c64_final

###carrington:

##SBATCH -M ukko
#SBATCH -M carrington
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=160G # memory per node 20G per task?

###vorna (not working yet):

##SBATCH -M vorna
##SBATCH --partition=short
##SBATCH --ntasks=1
##SBATCH --nodes=1
##SBATCH --cpus-per-task=4        # cpu-cores per task (>1 if multi-threaded tasks)
##SBATCH --mem=50G # memory per node 20G per task?




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


# key for -var flag
# 1. open vs. closed field boundary
# 2. proton differential energy flux
# 3. Field aligned currents (FACS)
# 4. Magnetopause position
# 5. dB/dt. Note: this needs to be evaluated in the ionosphere (once fully implemented by Urs)
# 6. other plots


## arguments

##   parser.add_argument('-nproc', default=1, help="number of processors to use " )
##   parser.add_argument('-startstop', nargs='*', help="2-element list, start and stop index (divided by deltanframes)" )
##   parser.add_argument('-deltanframes',  help="only analyze one in every delta_nframes" )
##   parser.add_argument('-run', default='EGL', help="the Vlasiator run, in all caps " )
##   parser.add_argument('-var', default = ['1','2','3','4','5','6'], nargs='*', help="a list of plot identifiers (numbers), set to which plots you want to make" )
##   parser.add_argument('-nphi', default=360, help="the number of longitudes in the lon-lat grid" )
##   parser.add_argument('-phimin', default = -180, help="minimum longitude, (range -180 to 180)" ) # note phi = (pi/180) * (phi_input-180)
##   parser.add_argument('-phimax', default = 180, help="maximum longitude, (range -180 to 180)" )
##   parser.add_argument('-nlat', default=180, help="the number of latitudes in the lon-lat grid" )
##   parser.add_argument('-latmin', default = -90, help="minimum longitude, (range -90 to 90)" ) # note lat = -lat_input * (pi/180) 
##   parser.add_argument('-latmax', default = 90, help="maximum longitude, (range -90 to 90)" )
##   parser.add_argument('-plot', action='store_true', help="set this flag to plot the data")  #default: plot=False when -plot flag not set
##   parser.add_argument('-save', action='store_true', help="set this flag to save the data into a .csv format")  #default: save=False when -save flag not set
##   parser.add_argument('-savefile', default='/wrk-vakka/users/horakons/carrington/data/test_data.csv', help="filename for the .csv data")
##   parser.add_argument('-append', action='store_true', help="set this flag to append to .csv file")  #default: append=False when -save flag not set





#run=EGI
#run=EGL
#run=EGP
run=FHAFGB  # note on FHAFBG (which contain fg_b data): simulation time is frame index *10
#run=EGILIKE2

#fileprefix=aurora_plot_data_south
fileprefix=aurora_plot_data
savdirtmp=$HOME_DIR/carrington/data/$run
savdir=$savdirtmp/$fileprefix

mkdir $savdir
echo $savdir

#python carrington.py 
#python carrington.py -run EGI -var 6
#python carrington.py -run EGI -var 1 2 3 4 6
#python carrington.py -nproc $1 -startstop $2 $3 -deltanframes 1 -run $run -var 1 -save -savefile "$savdir/$fileprefix.csv"
python carrington.py -nproc $1 -startstop $2 $3 -deltanframes 1 -run $run -var 1 2 3 -latmin 60 -latmax 90 -nlat 301 -nphi 361 -plot -save -savefile "$savdir/$fileprefix.csv"
#python carrington.py -nproc $1 -startstop $2 $3 -deltanframes 1 -run $run -var 1 2 3 -plot -latmin -90 -latmax -60 -nlat 301 -nphi 361 -save -suffix "_south" -savefile "$savdir/$fileprefix.csv"
#python carrington.py -nproc $1 -startstop $2 $3 -deltanframes 1 -run $run -var 3 -plot -latmin 60 -latmax 90 -nlat 31 -nphi 181


file=$savdir/$(printf "%s" $fileprefix)_$(printf "%04d" $2)_$(printf "%04d" $3).csv
echo $file
firstfile=$(printf "%s" $savdir)/$(printf "%s" $fileprefix)_$(printf "%04d" $2).csv
cat $firstfile > $file
rm $firstfile
i=$2
while [ $i -lt $3 ]
do
        echo $i
        ((i++))
        thisfile=$(printf "%s" $savdir)/$(printf "%s" $fileprefix)_$(printf "%04d" "$i").csv
        echo $thisfile
        tail -n +2 $thisfile >> $file    #copy everything except header line
        rm $thisfile
done






echo Job complete.
