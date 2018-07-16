prefix=$1

/usr/bin/time -v -o $prefix.kprocessor.time ./Kprocessor count -i $prefix.fastq -o $prefix.mqf  -k31 -t 4 -n ${prefix}_k31.hist
./Kprocessor dump -k 31 -i $prefix.mqf -o $prefix.kprocessor.text
sort -i $prefix.kprocessor.text -o  $prefix.kprocessor.text --parallel=4 -S 3000M -T /media/mostafa/Media/Datasets/Kprocessor/DataSets/

md5sum $prefix.kprocessor.text $prefix.dsk.text
cat $prefix.kprocessor.time
