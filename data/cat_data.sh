#EXAMPLE:
#./cat_data.sh directory_with_csv_files
#assumes all files have the same header

dir=$1
echo "concatenating files in input directory [ $dir ]"

rm $dir/data.csv
ls $dir/*.csv > $dir/tmp
head -n 1 $dir/tmp > $dir/tmp_firstfile
cat $dir/tmp_firstfile | while read f; do head -n 1 ${f} > $dir/data.csv; done;
cat $dir/tmp | while read f; do grep -v proton ${f} >> $dir/data.csv; done;
rm $dir/tmp
rm $dir/tmp_firstfile



