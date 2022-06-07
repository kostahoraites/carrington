#!/bin/bash -l
#SBATCH -t 02:00:00
#SBATCH -J array_movie
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=0-200
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=250G
#SBATCH -M carrington


# This script can be used on CSC/taito to generate an array job, 
# which renders multiple frames in order to e.g. make a movie.

run='EGL'

case $run in

  "EGI")
    frameStart=34  # set to initial frame   #EGI   662-1506 (multiply x20 because delta_nframes=20, see carrington.py)
    frameEnd=75 # set to the final frame    #EGI
    ;;

  "EGL")
    frameStart=32  # set to initial frame   #EGL   640-1760 (multiply x20 because delta_nframes =20, see carrington.py)
    frameEnd=88 # set to the final frame    #EGL
    #frameStart=36  # set to initial frame   #EGL   640-1760 (multiply x20 because delta_nframes =20, see carrington.py)
    #frameEnd=36 # set to the final frame    #EGL
    ;;

  "EGP")
    frameStart=269  # set to initial frame   #EGP  269-351 (delta_nframes = 1, no precipitation data)
    frameEnd=319 # set to the final frame    #EGP   
    #frameStart=352  # set to initial frame   #EGP  352-506 (delta_nframes = 1, includes precipitation data)
    #frameEnd=506 # set to the final frame    #EGP   
    ;;

  *) echo "no run specified!"
    ;;

esac




# How many jobs? SLURM_ARRAY_TASK_COUNT does not work on all systems
# so calculate job count (or set it manually to match the array
# argument given to sbatch).
jobcount=$(( $SLURM_ARRAY_TASK_MAX - $SLURM_ARRAY_TASK_MIN + 1 )) 

# find job array index
index=$(( $SLURM_ARRAY_TASK_ID - $SLURM_ARRAY_TASK_MIN ))

frameEndC=$(( $frameEnd + 1 )) # Need to iterate to 1 past final frame
totalFrames=$(( $frameEndC - $frameStart )) # Total frames to calculate
increment=$(( $totalFrames / $jobcount )) # amount of frames per job (rounded down)

# Calculate remainder
remainder=$(( $totalFrames - $jobcount * $increment ))

start=$(( $frameStart + $index * $increment ))
end=$(( $start + $increment ))

# Remainder frames are divvied out evenly among tasks
if [ $index -lt $remainder ];
then 
    start=$(( $start + $index ))
    end=$(( $end + $index + 1 ))
else
    start=$(( $start + $remainder ))
    end=$(( $end + $remainder ))
fi;

# Ensure final job gets correct last frame
if [ $SLURM_ARRAY_TASK_ID -eq $SLURM_ARRAY_TASK_MAX ];
then 
    echo Verifying final frame: $end $frameEndC
    end=$frameEndC
fi;

#echo Calculating a total of $totalFrames frames divided amongst $jobcount jobs.
#echo Using a remainder of $remainder.
#echo Current job id is $SLURM_ARRAY_TASK_ID and calculated frames are
#echo from $start to $end 


#module load mayavi2
module purge
module load Python/3.7.2-GCCcore-8.2.0
export PATH=/proj/jesuni/projappl/tex-basic/texlive/2020/bin/x86_64-linux:$PATH
#export PYTHONPATH=$PYTHONPATH:/appl/opt/Python/3.7.2-GCCcore-8.2.0/easybuild/python:/proj/horakons/analysator:/proj/horakons/.local/:/proj/horakons/.local/bin:/wrk/users/horakons/utils/:/wrk/users/horakons/carrington/
#module load matplotlib

export PTNONINTERACTIVE=1
#export PTNOLATEX=1
#export PTOUTPUTDIR=/wrk/users/markusb/Plots/

# key for -var flag
# 1. open vs. closed field boundary
# 2. proton differential energy flux
# 3. Field aligned currents (FACS)
# 4. Magnetopause position
# 5. dB/dt. Note: this needs to be evaluated in the ionosphere (once fully implemented by Urs)
# 6. other plots


#python carrington.py $start $end
#python carrington.py -run $run -var 1 2 3 4 6
#python carrington.py $start $end -run $run -var 6
echo $start
echo $end
python carrington.py -startstop $start $end -deltanframes 20 -run $run -var 1 2 3 4 -save -savefile $HOME_DIR/carrington/data/$run/aurora_plot_data/aurora_plot_data_$run_$start_$end.csv



#srun python carrington.py -nproc 8 -startstop $start $start+8 -deltanframes 1 -run EGL -var 2 -nphi 1 -phimin 150 -phimax 150 -append -save -savefile $HOME_DIR/carrington/data/aurora_keogram_data_egl_phi_150_8_julia.csv

#srun -n 6 carrington_multithread.sh $start $start+6
#carrington_multithread.sh $start $start+6


echo Job $SLURM_ARRAY_TASK_ID complete.
