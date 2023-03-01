#!/bin/zsh

mkdir -p ../200822_aggregated_fastqs






for dir in *barcode*

do

cat $dir/*.fastq.gz > ../../aggregated_fastqs/batch3/${dir}.fastq.gz

done







# for the complicated batch situations (Kenny16S)

mkdir -p 231122_aggregated_fastqs

for dir in *bat*/fastq_pass/*barcode*

do

barcode=$(basename $dir)

cat $dir/*.fastq.gz >> 231122_aggregated_fastqs/${barcode}.fastq.gz


done


mkdir -p 241122_aggregated_fastqs





for dir in *bat*/fastq_pass/*barcode*

do

barcode=$(basename $dir)

cat $dir/*.fastq.gz >> 241122_aggregated_fastqs/${barcode}.fastq.gz

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
