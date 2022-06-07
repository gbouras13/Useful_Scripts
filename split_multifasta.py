#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

def get_input():
	usage = 'python3 split_multifasta.py ...'
	parser = argparse.ArgumentParser(description='script to split multifasta into individual fasta files based on accession', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='output directory',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)


for dna_record in SeqIO.parse(args.infile, "fasta"):

    id = dna_record.id
    seq = dna_record.seq
    filename = id + '.fasta'

    outfile = os.path.join(args.outdir, filename)

    single_record = SeqRecord(seq, id, description="")

    with open(outfile, 'w') as nt_fa:
        SeqIO.write(single_record, nt_fa, 'fasta')
