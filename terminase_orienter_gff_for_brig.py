import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import pandas as pd
import numpy as np
import os

def get_input():
    usage = 'python3 terminase_orienter_gff_for_brig.py -i <input gff> -s <strand of terminase> -t <start of terminase coordinate> -p <prefix> -f <fasta> -o <output dir>...'
    parser = argparse.ArgumentParser(description='script to orient phage genome to begin with large terminase subunit', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', action="store", help='input file in gff format',  required=True)
    parser.add_argument('-s', '--strand', action="store", help='strandedness of terminase - pos or neg',  required=True)
    parser.add_argument('-t', '--terminase', action="store", help='terminase coordinate',  required=True)
    parser.add_argument('-p', '--prefix', action="store", help='output files prefix',  required=True) 
    parser.add_argument('-f', '--fasta', action="store", help='input fasta file',  required=True)
    parser.add_argument('-o', '--outdir', action="store", help='output directory',  required=True )
    args = parser.parse_args()

    return args

args = get_input()

# read in fasta
#record = SeqIO.read(args.fasta, "fasta")
record = SeqIO.read(args.fasta, "fasta")
# get length of the fasta
len = len(record.seq)

# read in cleaned gff
colnames=['contig', 'program', 'type', 'start', 'end', 'score', 'strand', 'n', 'description'] 
gff_df = pd.read_csv(args.infile, delimiter= '\t', index_col=False, header=None, names=colnames)

# remove hypothetical protein 
gff_df['description'] = gff_df['description'].str.replace('hypothetical protein', '', regex=True)

strand = args.strand
# add one
if strand == "neg":
    terminase = int(args.terminase) + 1
if strand == "pos":
    terminase = int(args.terminase) - 1


# cut seqeunce
gff_df['start2'] = np.where(gff_df['start'] > terminase, gff_df['start'] - terminase, gff_df['start']+len -terminase)
gff_df['end2'] = np.where(gff_df['end'] > terminase, gff_df['end'] - terminase, gff_df['end']+len -terminase)



if strand == "neg":
    # tmp variable for the exchange of start and end
    tmpend = len - gff_df['start2']
    tmpstart = len - gff_df['end2']
    gff_df['start2'] = tmpstart
    gff_df['end2'] = tmpend
    # change strandedness
    gff_df['strand'] = np.where(gff_df['strand'] == '+', '-', '+')

# sort
gff_df = gff_df.sort_values(by=['start2'])

# get trnas

outdir = args.outdir
prefix = args.prefix

if not os.path.exists(outdir):
  os.makedirs(outdir)

trna_df = gff_df[gff_df['type'] == 'trna']
trna_df = trna_df[["start2", "end2", "description"]]
trna_df.to_csv(os.path.join(outdir, prefix + '_trna.txt'), sep="\t", index=False, header=False)

# cds
cds_df = gff_df[gff_df['type'] != 'trna']
fwd_df = cds_df[cds_df['strand'] == '+']
rev_df = cds_df[cds_df['strand'] == '-']

fwd_df = fwd_df[["start2", "end2", "description"]]
rev_df = rev_df[["start2", "end2", "description"]]
fwd_df.to_csv(os.path.join(outdir, prefix + '_fwd.txt'), sep="\t", index=False, header=False)
rev_df.to_csv(os.path.join(outdir, prefix + '_rev.txt'), sep="\t", index=False, header=False)




	
