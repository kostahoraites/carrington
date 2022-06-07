#!/bin/bash

#example call (carrington): ./carrington_multithread_batch.sh 16 641 656
#example call (vorna): ./carrington_multithread_batch.sh 4 1041 1044

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




