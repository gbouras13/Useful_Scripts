#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_input():
	usage = 'python3 trim_uniprot_header.py ...'
	parser = argparse.ArgumentParser(description='script to get the accession only from uniprot multifasta', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

with open(args.outfile, 'w') as aa_fa:
    for dna_record in SeqIO.parse(args.infile, "fasta"):
        # split the list
        spl_list = dna_record.id.split("|")

        # get second item (accession)

        acc = spl_list[1]
                # write new record no description
        aa_record = SeqRecord(dna_record.seq, id=acc, description="")
        SeqIO.write(aa_record, aa_fa, 'fasta')