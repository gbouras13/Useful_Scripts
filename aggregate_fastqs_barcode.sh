#!/bin/zsh

mkdir -p ../aggregated_fastqs

for dir in *barcode*

do

cd $dir
cat *.fastq.gz > ${dir}.fastq.gz
cp  ${dir}.fastq.gz ../../aggregated_fastqs
cd ..
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
