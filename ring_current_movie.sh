#!/bin/bash -l
#SBATCH -t 00:20:00
#SBATCH -J array_movie
#SBATCH -p short
#SBATCH -n 1
#SBATCH --array=0-200
#SBATCH --no-requeue
#SBATCH --mem-per-cpu=16000

# This script can be used on CSC/taito to generate an array job, 
# which renders multiple frames in order to e.g. make a movie.

#frameStart=32  # set to initial frame
#frameStart=32  # set to initial frame
frameStart=32  # set to initial frame   #EGL   640-1760 (delta_nframes =20)
frameEnd=88 # set to the final frame    #EGL
#frameStart=269  # set to initial frame   #EGP  269-319 (delta_nframes = 1)
#frameEnd=319 # set to the final frame    #EGP   

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
module load matplotlib

export PTNONINTERACTIVE=1
#export PTNOLATEX=1
export PTNOLATEX=
export PTOUTPUTDIR=/wrk/users/markusb/Plots/

python ring_current.py $start $end
echo Job $SLURM_ARRAY_TASK_ID complete.
