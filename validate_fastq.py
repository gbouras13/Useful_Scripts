#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import os
import gzip
import sys

def get_input():
	usage = 'python validate_fastq.py ...'
	parser = argparse.ArgumentParser(description='Script to validate fastq files that have some corrupt records (length of sequence not matching quality).', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='Input file in fastq format - needs to be .fastq or .fastq.gz suffixed',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='Outfile file in gzipped fastq format',  required=True)
	args = parser.parse_args()
	return args

# https://groverj3.github.io/articles/2019-08-22_just-write-your-own-python-parsers-for-fastq-files.html

def validate_fastq(input_file):
    filename, file_extension = os.path.splitext(input_file)
    if file_extension == ".gz":
        with gzip.open(input_file, 'rt') as input_handle:
            record = []
            n = 0
            for line in input_handle:
                if n == 0 and line[:1] != "@": # if the first char isn't @, continue
                    record = []
                    continue
                n += 1
                record.append(line)
                if n == 4:
                    if len(record[1]) != len(record[3]): # if the lengths dont match
                        n = 0
                        record = []
                        continue
                    else:
                        yield record
                        n = 0
                    record = []
    elif file_extension == ".fastq":
        with open(input_file, 'r') as input_handle:
            record = []
            n = 0
            for line in input_handle:
                if n == 0 and line[:1] != "@": # if the first char isn't @, continue
                    n = 0
                    record = []
                    continue
                n += 1
                record.append(line)
                if n == 4:
                    if len(record[1]) != len(record[3]): # if the lengths dont match between sequence and quality
                        n = 0
                        continue
                    else:
                        yield record
                        n = 0
                    record = []
    else:
        sys.exit("neither .fastq or .fastq.gz input file detected")



if __name__ == "__main__":

    args = get_input()

    print("All output will be gzipped - setting an output file suffixed .fastq.gz is recommended.")

    # generator
    generator = validate_fastq(args.infile)

    ### writes out fastq gzipped by default
    with gzip.open(args.outfile, 'wt') as f:
        for x in generator:
            f.write(x[0])
            f.write(x[1])
            f.write(x[2])
            f.write(x[3])


