startfile=621
endfile=1760
#startfile=901
#endfile=1000
nodes=10       #nodes must divide evenly into endfile-startfile+1
numfiles=$((($endfile-$startfile+1)/$nodes))

batchfile=job_car.sh

for (( i = 0; i < $nodes; ++i )); do
  d=$(( numfiles*i ))
  thisstartfile=$(( startfile+d ))
  thisendfile=$(( thisstartfile+numfiles-1 ))
  #thisendfile=$(( thisstartfile+numfiles ))
  thisstartfile=$(printf "%07d\n" $thisstartfile)
  thisendfile=$(printf "%07d\n" $thisendfile)
  echo $thisstartfile $thisendfile
  thisbatchfile=run_${i}.sh
  cp $batchfile $thisbatchfile
  sed -i s/STARTFILE/$thisstartfile/g $thisbatchfile
  sed -i s/ENDFILE/$thisendfile/g $thisbatchfile
  thisoutfile=run_${i}_outfile
  sed -i s/OUTFILE/$thisoutfile/g $thisbatchfile
  sbatch $thisbatchfile
done

rm run_*.sh

