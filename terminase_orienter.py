#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_input():
	usage = 'python3 terminase_orienter.py ...'
	parser = argparse.ArgumentParser(description='script to orient phage genome to begin with large terminase subunit', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-s', '--strand', action="store", help='strandedness of terminase - pos or neg',  required=True)
	parser.add_argument('-t', '--terminase', action="store", help='terminase coordinate',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='output file',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

record = SeqIO.read(args.infile, "fasta")

# get length of the fasta
len = len(record.seq)

print(record.seq )

# positive it's easy

# reorient to start at the terminase  
if args.strand == "pos":
	start = record.seq[(int(args.terminase) - 1):len]
	end = record.seq[0:int(args.terminase)-1]
	total = start + end



# revese compliment if the strand is negative

if args.strand == "neg":
	record.seq = record.seq.reverse_complement()
	start = record.seq[(len - int(args.terminase)):len]
	end = record.seq[0:(len - int(args.terminase))]
	total = start + end




record.seq = total

print(record.seq )

SeqIO.write(record, args.outfile, "fasta")
	
