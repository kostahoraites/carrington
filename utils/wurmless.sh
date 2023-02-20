# perform less on an input filename that comes from a list command (or list of filenames?)
# ex. while in this directory, from the command line, call: ls wu* | ./wurmless.sh
# this should perform less on this file, wurmless.sh
# the whole point is to be able to read slurm files easier, see wurm alias in .bashrc file
while read first $1;
   do less $first;
done
