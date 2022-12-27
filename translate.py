#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_input():
	usage = 'python3 translate.py ...'
	parser = argparse.ArgumentParser(description='script to translate a multifasta nucleotide fasta to amino acid', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

with open(args.outfile, 'w') as aa_fa:
    for dna_record in SeqIO.parse(args.infile, "fasta"):
        #print(dna_record.id)
	# use both fwd and rev sequences
        dna_seqs = [dna_record.seq, dna_record.seq.reverse_complement()]

        # generate all translation frames
        aa_seqs = (s[i:].translate(to_stop=True) for i in range(3) for s in dna_seqs)

        # select the longest one
        max_aa = max(aa_seqs, key=len)

        # write new record
        aa_record = SeqRecord(max_aa, id=dna_record.id, description="translated sequence")
        SeqIO.write(aa_record, aa_fa, 'fasta')
