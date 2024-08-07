#!/usr/bin/env python3
from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter

def get_input():
	usage = 'python3 get_fasta_from_genbank.py ...'
	parser = argparse.ArgumentParser(description='script to quickly get FASTA from GBK', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in genbnak format',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='output file in FASTA format', required=True )
	args = parser.parse_args()

	return args

args = get_input()

def extract_fasta_from_gbk(gbk_file, fasta_file):
    """
    Extract FASTA sequences from a GenBank (GBK) file and save them to a FASTA file.
    :param gbk_file: Path to the input GBK file.
    :param fasta_file: Path to the output FASTA file.
    """
    with open(gbk_file, "r") as input_handle, open(fasta_file, "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        count = SeqIO.write(sequences, output_handle, "fasta")
        print(f"Converted {count} records")

# Usage example

extract_fasta_from_gbk(args.infile, args.outfile)
