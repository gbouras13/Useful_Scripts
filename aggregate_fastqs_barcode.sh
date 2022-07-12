#!/bin/zsh

mkdir -p ../aggregated_fastqs

for dir in *barcode*

do

cat $dir/*.fastq.gz > ../aggregated_fastqs/${dir}.fastq.gz

done

### if error run again

# mkdir -p ../aggregated_fastqs
#
# for dir in *barcode*
# do
# echo $dir
# cd $dir
# cp  ${dir}.fastq.gz ../../aggregated_fastqs/${dir}.fastq.gz
# cd ..
# done
