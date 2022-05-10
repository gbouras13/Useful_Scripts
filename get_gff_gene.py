import argparse
from argparse import RawTextHelpFormatter
import pandas as pd



def get_input():
    usage = 'python3 get_gff_gene.py ...'
    parser = argparse.ArgumentParser(description='script to output the gene of a coordinate in a gff', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--infile', action="store", help='input file in gff format',  required=True)
    parser.add_argument('-c', '--contig', action="store", help='contig',  required=True)
    parser.add_argument('-p', '--position', action="store", help='position',  required=True)
    args = parser.parse_args()

    return args

args = get_input()

#record = SeqIO.read(args.infile, "fasta")
contig = args.contig
position = int(args.position)
gff = args.infile


# gff="/Users/a1667917/Documents/Total_Staph/gffs/D32.gff"
# contig=1
# position=10005

colnames=['contig', 'evidence', 'type', 'start', 'end', 'score', 'strand', 'n', 'description'] 
# if file is empty


try:
    gff_df = pd.read_csv(gff, delimiter= '\t', index_col=False, header=None, names=colnames)
except pd.errors.EmptyDataError:
    print('gff is empty')

# get only the rows off gff 
    types = ['CDS', 'tRNA', 'rRNA', 'tmRNA']
    gff_df = gff_df[gff_df['type'].isin(types)]

# find the row that has the position and contig

# selecting rows based on condition
selected_row_df = gff_df.loc[(gff_df['start'] <= position) & (gff_df['end'] >= position) & (gff_df['contig'] == str(contig) ) ]

# define if empty
empty = False

# if empty - non coding region
if selected_row_df.empty:
    empty = True
    print('DataFrame is empty!')

# function to get the locus tag (between ID= and the first semi colon)
def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


if empty == False:
# get is_name and locus_tage
    selected_row_df['locus_tag'] = selected_row_df['description'].apply(lambda x: find_between(x,"ID=", ";"  ) )
    selected_row_df['Gene_Name'] = selected_row_df['description'].apply(lambda x: find_between(x,"Name", ";"  ) )
    selected_row_df[['description','product']] = selected_row_df['description'].str.split('product=',expand=True)
    print(selected_row_df)

if empty == True:
    print("This position is in a non-coding region")

# # write to csv
# abundance_df.to_csv(output, sep=",", index=False, header=False)
# else:
# colnames=['contig', 'evidence', 'type', 'start', 'end', 'score', 'strand', 'n', 'description', 'locus_tag', 'IS_name', 'product'] 
# abundance_df = pd.DataFrame(columns=colnames)
# abundance_df.to_csv(output, sep=",", index=False, header=False)


