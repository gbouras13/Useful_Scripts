import argparse
from argparse import RawTextHelpFormatter
import os 
import shutil
import gzip
import csv

'''
requires 2 column input csv

`Sample,barcode`

If you have samples with multiple runs, the code will detect and append to the same {sample}.fastq.gz

'''


def get_input():
    usage = 'python3 aggregate_nanopore_fastqs.py ...'
    parser = argparse.ArgumentParser(description='script to aggreagted nanopore FASTQs and label with sample name. Will append already existing files, so be careful!', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', action="store", help='input csv file in csv format. 2 columns with sample names an barcode (integers). Needs 2 columns: Sample and Barcode. Example Sample,Barcode\nsample,01 ',  required=True)
    parser.add_argument('-d', '--directory', action="store", help='fastq_pass directory path',  required=True)
    parser.add_argument('-o', '--output', action="store", help='output directory where the aggregated FASTQs will be saved',  required=True)
    parser.add_argument('--dorado', action="store_true", help='use if dorado rebasecalled, will be only 1 fastq')
    args = parser.parse_args()

    return args

args = get_input()

# read in all the 
csv_file = args.infile
output_dir = args.output
fastq_path = args.directory
data_dict = {}

if os.path.isdir(output_dir) == False:
    os.mkdir(output_dir)

with open(csv_file, mode='r', encoding='utf-8-sig') as file:
    csv_reader = csv.DictReader(file)
    for row in csv_reader:
        sample_id = row['Sample']
        # to add barcode to the front and 2 digits for the number
        barcode = f"barcode{int(row['Barcode']):02}"
        data_dict[sample_id] = barcode


# loop over the dictionary
for sample_id, barcode in data_dict.items():

    print(f"aggregating {barcode} as sample {sample_id}")
    # Create the output file path
    output_file = os.path.join(output_dir, f'{sample_id}.fastq.gz')

    if args.dorado is True:
        # Find and concatenate the corresponding .fastq.gz files
        fastq_files = [f for f in os.listdir(fastq_path) if f.endswith('.fastq.gz') and barcode in f]
        if len(fastq_files) == 0:
            raise FileNotFoundError(f"No fastq file containing the barcode {barcode} found in {fastq_path}")
        else:
            barcode_path = fastq_path
            for f in fastq_files:
                print(f"file {f} found for {barcode}")            
    else: # normal off the sequencer
        # there will be lots of fastq files
        # Find and concatenate the corresponding .fastq.gz files
        barcode_path = os.path.join(fastq_path, barcode)
        fastq_files = [f for f in os.listdir(barcode_path) if f.endswith('.fastq.gz')]

    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"{output_file} exists. Appending the new FASTQs to this file.")
        # If it exists, open in "append" mode
        with gzip.open(output_file, 'ab') as output:
            for fastq_file in fastq_files:
                out_path = os.path.join(barcode_path, fastq_file)
                with gzip.open(out_path, 'rb') as input_file:
                    for line in input_file:
                        output.write(line)
    else:
        # If it doesn't exist, open in "write" mode
        with gzip.open(output_file, 'wb') as output:
            for fastq_file in fastq_files:
                out_path = os.path.join(barcode_path, fastq_file)
                with gzip.open(out_path, 'rb') as input_file:
                    for line in input_file:
                        output.write(line)
