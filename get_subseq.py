#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_input():
	usage = 'python3 get_subseq.py ...'
	parser = argparse.ArgumentParser(description='script to quickly get subsequence from fasta and coordinates', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-s', '--start', action="store", help='first coordinate' )
	parser.add_argument('-e', '--end', action='store', help='second coordinate' )
	args = parser.parse_args()

	return args

args = get_input()

record = SeqIO.read(args.infile, "fasta")

print(record.seq[(int(args.start) - 1):int(args.end)])