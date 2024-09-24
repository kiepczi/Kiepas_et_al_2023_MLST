"""This script was used to generate data which will provide
information about total number of all and unique copies of marker
alleles per each genome.
"""

#Set Up
import pandas as pd
from pathlib import Path

#Load MLST tool outputs and combining them
MLST_all = pd.read_csv(Path('../../supplementary_file_4/output/revised_scheme_run/revised_scheme_run.txt'), sep='\t', usecols=['FILE', 'SCHEME', 'ST', '16S', 'atpD', 'gyrB', 'recA', 'rpoB', 'trpB'])



del MLST_all['SCHEME']
del MLST_all['ST']


#Next getting accession for each row
def get_accession(row):
    """Return genome accession"""

    try:
        accession = '_'.join(row['FILE'].split('/')[-1].split('_')[:2])
    except IndexError:
        accession = 'NA'
    
    return accession
MLST_all['accession'] = MLST_all.apply(get_accession, axis=1)

del MLST_all['FILE']


#First getting total copies of each marker per genome
def count_all_copies(row, col):
    """Return count of total number of marker genes"""
    try:
        if '-' in str(row[col]):
            copies = 0
        elif '?' in str(row[col]):
            copies = 0
        else:
            copies = str(row[col]).count(',') + 1
    except IndexError:
        copies = 'NA'
    
    return copies



#Next getting total number of unique copies of marker genes per genome
def count_unique_copies(row, col):
    """Return count of unique copies of marker genes"""
    try:
        if '-' in str(row[col]):
            copies = 0
        elif '?' in str(row[col]):
            copies = 0
        elif '~' in str(row[col]):
            copies = 0
        else:
            copies = len(set([_ for _ in str(row[col]).split(',')]))
    except IndexError:
        copies = 'NA'
    
    return copies


markers = ['16S', 'atpD', 'gyrB', 'recA', 'rpoB', 'trpB']

for _ in markers:
    MLST_all[f'total_{_}_copies'] = MLST_all.apply(count_all_copies, col=_, axis=1)
    MLST_all[f'unique_{_}_copies'] = MLST_all.apply(count_unique_copies, col=_, axis=1)
    del MLST_all[_]


MLST_all.to_csv(Path("../output/Genomes_allele_count.csv").expanduser(), index=False)


