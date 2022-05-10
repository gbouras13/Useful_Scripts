import argparse
from argparse import RawTextHelpFormatter
from Bio.Seq import Seq
import Bio.motifs as motifs
from Bio import SeqIO



def get_input():
	usage = 'python3 motif_finder.py ...'
	parser = argparse.ArgumentParser(description='script to output coordinates of motif in a fasta file', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
	parser.add_argument('-m', '--motif', action="store", help='find motif',  required=True)
	#parser.add_argument('-o', '--outfile', action="store", help='output file',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

#record = SeqIO.read(args.infile, "fasta")
motif = args.motif


instances = [Seq(motif)]
m = motifs.create(instances)

reads = list(SeqIO.parse(args.infile, "fasta"))

for i in range(len(reads)):
    for pos, seq in m.instances.search(reads[i].seq):
        print("%i %i %s" % ((i+1), pos, seq))