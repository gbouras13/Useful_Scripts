#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_input():
	usage = 'python3 chromosome_sta.py ...'
	parser = argparse.ArgumentParser(description='script to add sample name to contigs of fasta file', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

sample = args.infile.split('.')[0]

organism = "Staphylococcus aureus"

with open(args.outfile, 'w') as aa_fa:
    for dna_record in SeqIO.parse(args.infile, "fasta"):
        
        id_new = sample
        description_new = "[organism=" + organism +"] [location=chromosome] [topology=circular] [completeness=complete]"

        # write new record no desc
        dna_record_new = SeqRecord(dna_record.seq, id=id_new,  description=description_new, name="")
        SeqIO.write(dna_record_new, aa_fa, 'fasta')