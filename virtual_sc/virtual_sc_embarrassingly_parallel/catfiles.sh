cat run_0_outfile.txt | sed -n 3p > runs.txt
for i in run_*_outfile.txt
do
tail -n +4 $i  >> runs.txt
done


