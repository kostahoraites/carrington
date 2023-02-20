#!/bin/bash

#example call (carrington): ./carrington_multithread_batch.sh 16 641 768
#example call (vorna): ./carrington_multithread_batch.sh 4 1041 1044

# WORKFLOW TO MAKE PAPER PLOTS
#cdwrk
#cd carrington
#./carrington_multithread_batch.sh 16 621 1760    # check if files have same length:   wc data/EGL/aurora_plot_data/*.csv
# /wrk-vakka/users/horakons/carrington/data/cat_data.sh /wrk-vakka/users/horakons/carrington/data/EGL/aurora_plot_data/
# sbatch carrington_plot_timeseries.sh




nthreads=$1
start=$2        
stop=$3

nframes=$(($start-$stop+1))
njobs=$(($nframes/$nthreads+1))

thisstart=$start

while [ $thisstart -le $stop ]
do 
  thisend=$(($thisstart+$nthreads-1))
  echo $thisstart $thisend
  #sbatch carrington_multithread_aurora.sh $nthreads $thisstart $thisend
  sbatch carrington_multithread.sh $nthreads $thisstart $thisend
  thisstart=$(($thisend+1))
done




