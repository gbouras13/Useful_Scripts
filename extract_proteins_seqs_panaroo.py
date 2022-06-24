#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from Bio.Seq import Seq


def get_input():
	usage = 'python3 extract_protein_seqs_panaroo.py ...'
	parser = argparse.ArgumentParser(description='script to get AA sequences from panaroo output in fasta format', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file of desired rows from gene_data.csv in csv (with header row)',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
	args = parser.parse_args()

	return args

args = get_input()

# read the csv in 

colnames=['gff_file', 'scaffold_name', 'clustering_id', 'annotation_id', 'prot_sequence', 'dna_sequence', 'gene_name', 'description'] 
# if file is empty

try:
    gene_data_df = pd.read_csv(args.infile, delimiter= ',', index_col=False, header=None, names=colnames, skiprows = 1)
except pd.errors.EmptyDataError:
    print('csv is empty')


with open(args.outfile, 'w') as aa_fa:
    for index, row in gene_data_df.iterrows():
        sequence = Seq(row["prot_sequence"])
        aa_record = SeqRecord(seq=sequence, id=str(row["gff_file"])+","+str(row["annotation_id"]), description="")
        SeqIO.write(aa_record, aa_fa, 'fasta')
