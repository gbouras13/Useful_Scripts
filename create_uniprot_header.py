#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv

def get_input():
    usage = 'python3 create_uniprot_header.py ...'
    parser = argparse.ArgumentParser(description='script to get the accession only from uniprot multifasta', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', action="store", help='input file in fasta format',  required=True)
    parser.add_argument('-o', '--outfile', action="store", help='outfile file in fasta format',  required=True)
    parser.add_argument('-c', '--csv', action="store", help='gene_presence.csv from roary',  required=True)
    parser.add_argument('-s', '--species', action="store", help='species of pangenome e.g. "Staphylococcus aureus" needs quotes for the space ',  required=True)
    parser.add_argument('-t', '--taxid', action="store", help='taxid of pangenome e.g. "1280" ',  required=True)
    parser.add_argument('-p', '--program', action="store", help='program used to create the pangenome: roary or panaroo',  required=True, default='roary')
    args = parser.parse_args()
    return args


# get input
args = get_input()

# read in gene_presence as csv
reader = csv.reader(open(args.csv))
# set as dictionary
result = {}
for row in reader:
    key = row[0]
    if key in result:
        # implement your duplicate row handling here
        pass
    result[key] = row[2]
#print(result)


# roary 

# open fasta

if args.program == 'roary':
    with open(args.outfile, 'w') as dna_fa:
        for dna_record in SeqIO.parse(args.infile, "fasta"):
            # split the list
            spl_list = dna_record.description.split(" ")
            #print(dna_record.description)
            # get first item (locustag)
            locus = spl_list[0]
            # get second item gene name
            identifier = spl_list[1]
            
            # get description from gene_presence.csv dictionary
            description = result[identifier]
            
            id_uniprot = "ro|"+identifier + "|" + locus + " " + description + " OS=" + args.species + " OX=" + args.taxid + " PE=4 SV=1"

            # translate 
            s = dna_record.seq
            aa_seq = s[0:].translate(to_stop=True) 

            # write new record no description
            dna_record_uniprot = SeqRecord(aa_seq, id=id_uniprot, description="")
            SeqIO.write(dna_record_uniprot, dna_fa, 'fasta')
if args.program == 'panaroo':
    with open(args.outfile, 'w') as dna_fa:
        for dna_record in SeqIO.parse(args.infile, "fasta"):
            # split the list
            spl_list = dna_record.description.split(" ")
            #print(dna_record.description)
            # get  gene name
            identifier = spl_list[0]

            # no locus tag for panaroo pangenome
            
            # get description from gene_presence.csv dictionary
            description = result[identifier]
            
            id_uniprot = "ro|"+identifier + "|"  + description + " OS=" + args.species + " OX=" + args.taxid + " PE=4 SV=1"

            # translate 
            s = dna_record.seq
            aa_seq = s[0:].translate(to_stop=True) 

            # write new record no description
            dna_record_uniprot = SeqRecord(aa_seq, id=id_uniprot, description="")
            SeqIO.write(dna_record_uniprot, dna_fa, 'fasta')

        
print("uniprot headers created")