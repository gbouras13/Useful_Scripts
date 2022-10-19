
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from typing import List
import argparse
import pandas as pd
from BCBio import GFF




def get_input():
	parser = argparse.ArgumentParser(description='calculate coding density from gbk')
	parser.add_argument('-i', '--infile', action="store", help='Input genome file in genbank format.',  required=True)
	parser.add_argument('-o', '--outfile', action="store", help='file to write the output to .', required=True )

	args = parser.parse_args()

	return args


def parse_gff(filename):
    """
    Parse a gff file and return a Bio::Seq object
    """
    handle = open(filename, 'r')
    return GFF.parse(handle) 
    

def get_coding_density(seqiorec: SeqRecord, ftypes: List[str] = ['CDS', 'tRNA', 'rRNA']) -> float:
    """
    Get coding density for a SeqRecord considering given features types
    :param seqiorec: SeqRecord object
    :param ftypes: a list of feature types
    :return:
    Taken from https://github.com/npbhavya/Phicore/tree/main/PhicoreModules 
    """

    cdcov = np.zeros(len(seqiorec.seq))
    for feature in seqiorec.features:
        if feature.type in ftypes:
            start, stop = map(int, sorted([feature.location.start, feature.location.end]))
            cdcov[start:stop] += 1
    return sum([1 if x > 0 else 0 for x in cdcov]) / len(seqiorec.seq)


# get the input
args = get_input()

cds =[]
contigs =[]

for record in parse_gff(args.infile):
    contig = record.id
    cd = get_coding_density(record, ['CDS']) * 100
    contigs.append(contig)
    cds.append(cd)


cd_df = pd.DataFrame(
{'contig': contigs,
    'coding_density': cds,
})

cd_df.to_csv(args.outfile, sep="\t", index=False)

