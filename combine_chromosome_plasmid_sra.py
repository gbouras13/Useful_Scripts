#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import shutil
import sys



# https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/#batch_assignment
def get_input():
    usage = 'python3 plasmid_header_sra.py ...'
    parser = argparse.ArgumentParser(description='script to add descriptions to plasmids header of fasta file', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-c', '--chromdir', action="store", help='chromosome directory',  required=True)
    parser.add_argument('-o', '--outdir', action="store", help='output directory',  required=True)
    parser.add_argument('-p', '--plasdir', action="store", help='plasmid directory',  required=True)
    parser.add_argument('-f', '--force', help="Overwrites the output directory.", action="store_true" )
    args = parser.parse_args()
    return args

args = get_input()

samples = []




	# remove outdir on force
if args.force == True:
    if os.path.isdir(args.outdir) == True:
        shutil.rmtree(args.outdir)
    else:
        print("\n--force was specified even though the outdir does not already exist. Continuing \n")
else:
    if os.path.isdir(args.outdir) == True:
        sys.exit("\nOutput directory already exists and force was not specified. Please specify -f or --force to overwrite the output directory. \n")  

if os.path.isdir(args.outdir) == False:
		os.mkdir(args.outdir)


# gets  samples and prints chromosomes to new file
for chrom_file in os.listdir(args.chromdir):
    if chrom_file.endswith(".fasta"):
        sample = chrom_file.replace(".fasta", "")
        samples.append(sample)
        outfile = os.path.join(args.outdir, sample + ".fasta")
        with open(outfile, 'w') as dna_fa:
            for dna_record in SeqIO.parse(os.path.join(args.chromdir,chrom_file), "fasta"):
                SeqIO.write(dna_record, dna_fa, 'fasta')

for plas_file in os.listdir(args.plasdir):
    if plas_file.endswith(".fasta"):
        # get sample - will be the first split by underscore
        sample = plas_file.split("_")[0]
        # get file
        outfile = os.path.join(args.outdir, sample + ".fasta")
        # append plasmid
        with open(outfile, 'a') as dna_fa:
            for dna_record in SeqIO.parse(os.path.join(args.plasdir,plas_file), "fasta"):
                SeqIO.write(dna_record, dna_fa, 'fasta')









# def samplesFromDirectory(dir):
#     """Parse samples from a directory"""
#     samples, extensions = glob_wildcards(os.path.join(dir,'{sample}_R1{extensions}'))
#     if not extensions:
#         sys.stderr.write("\n"
#                          "    FATAL: We could not parse the sequence file names from the specified directory.\n"
#                          "    We are expecting {sample}_R1{extension}, and so your files should contain the \n"
#                          "    characters '_R1' in the fwd reads and '_R2' in the rev reads. \n"
#                          "    Alternatively you can specify a 3-column TSV file instead to declare the sample\n"
#                          "    names and corresponding R1/R2 files. e.g. \n"
#                          "    sample1\tpath/to/reads/sample1.1.fastq.gz\tpath/to/reads/sample1.2.fastq.gz\n"
#                          "    sample2\tpath/to/reads/sample2.1.fastq.gz\tpath/to/reads/sample2.2.fastq.gz\n"
#                          "    ..."
#                          "    See https://hecatomb.readthedocs.io/en/latest/usage/#read-directory for more info\n"
#                          "\n")
#         sys.exit(1)
#     else:
#         for sample in samples:
#             outDict[sample] = {}
#             R1 = os.path.join(dir,f'{sample}_R1{extensions[0]}')
#             R2 = os.path.join(dir,f'{sample}_R2{extensions[0]}')
#             if os.path.isfile(R1) and os.path.isfile(R2):
#                 outDict[sample]['R1'] = R1
#                 outDict[sample]['R2'] = R2
#             else:
#                 sys.stderr.write("\n"
#                                  "    FATAL: Error globbing files. One of:\n"
#                                  f"    {R1} or\n"
#                                  f"    {R2}\n"
#                                  "    does not exist. Ensure consistent _R1/_R2 formatting and file extensions."
#                                  "\n")
#                 sys.exit(1)
#     return outDict


# organism = "Staphylococcus aureus"

# with open(args.outfile, 'w') as aa_fa:
#     for dna_record in SeqIO.parse(args.infile, "fasta"):
        
#         id_new = sample
#         description_new = "[organism=" + organism +"] [plasmid-name=" + sample + "] [topology=circular] [completeness=complete]"

#         # write new record no desc
#         dna_record_new = SeqRecord(dna_record.seq, id=id_new,  description=description_new, name="")
#         SeqIO.write(dna_record_new, aa_fa, 'fasta')


# # >contig02 [organism=Clostridium difficile] [strain=ABDC] [plasmid-name=pABDC1] [topology=circular] [completeness=complete]