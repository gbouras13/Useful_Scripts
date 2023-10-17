#!/bin/zsh

mkdir -p ../../160823_KAPA_aggregated_fastqs



mkdir -p ../170823_Mock_James_aggregated_fastqs


for dir in *barcode*

do

cat $dir/*.fastq.gz > ../170823_Mock_James_aggregated_fastqs/${dir}.fastq.gz

done


for dir in *barcode*

do

cat $dir/*.fastq.gz > ../../08_06_23_Plate2_aggregated_fastqs/${dir}.fastq.gz

done


mkdir -p ../140823_aggregated_fastqs_issy


for dir in *barcode*

do

cat $dir/*.fastq.gz > ../140823_aggregated_fastqs_issy/${dir}.fastq.gz

done

tcga_cptac_taxonomic_profiler binning --threads 1 --input tcga_trimnami_output/fastp   --output TCGA_output_all --rerun-triggers mtime --dry-run

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

ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Sinsheimervirus/annotation_releases



https://osf.io/download/tswbp/


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
