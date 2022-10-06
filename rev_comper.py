#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_input():
	usage = 'python3 rev_comper.py ...'
	parser = argparse.ArgumentParser(description='script to reverse compliment fasta (all contigs)', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

with open(args.outfile, 'w') as out_fa:
    for dna_record in SeqIO.parse(args.infile, "fasta"):
        dna_record.seq = dna_record.seq.reverse_complement()
        SeqIO.write(dna_record, out_fa, 'fasta')