"""This script was used to verify whether the alleles sequences obtained by MLST tool
contain protein coding sequences. 

Here, novel alleles returned by MLST with no protein sequence (eg. pseudogenes) will be discarded.
"""

#Set Up
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import os
import re

###Generating dataframe which will provide as a guide for blastn analysis

#Loading and combining data from most recent MLST analysis
MLST_output = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/strep_mlst11.txt").expanduser(), sep='\t')
MLST_replaced_genomes = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/mlst_replaced_genomes.txt").expanduser(), sep='\t')

MLST_data = pd.concat([MLST_output, MLST_replaced_genomes])

#Removing columns that are not of interest
del MLST_data['16S']
del MLST_data['ST']
del MLST_data['SCHEME']

#melting data
MLST_data = pd.melt(MLST_data, id_vars=["FILE"], 
                  var_name="marker", value_name='variant')



#unique 16S copies
def unique_copies(row, col):
    """Return unique allele copies """
    try:
        rRNA = str(set(int(s) for s in re.findall(r'-?\d+\.?\d*', str(row[col]))) if ',' in str(row[col]) else str(row[col])).replace('{', '').replace('}', '')
    except IndexError:
        rRNA = 'NA'

    return rRNA


MLST_data['variant'] = MLST_data.apply(unique_copies, col='variant', axis=1)

#Removing rows with genomes that have partial, novel or missing allele
MLST_data = MLST_data[~MLST_data['variant'].isin(list(set([_ for _ in MLST_data['variant'] if '?' in str(_) or '~' in str(_) or '-' in str(_) or ',' in str(_)])))] 

#Running blastn
datadir = Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/").expanduser()
filenames = sorted(datadir.glob("strep*/GCF*/*.fna"))


def run_blastn(row):
    """Retrun alleles with partial alleles (42?) replaces as '-'
    for specified df col"""


    accession = '_'.join(row['FILE'].split('/')[-1].split('_')[:2])
    marker = row['marker']+'_'+str(row['variant'])
    file = [_ for _ in filenames if accession in str(_) and 'from' not in str(_)][0]
    records = [_ for _ in list(SeqIO.parse(Path(f"~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/schemes/strep_mlst_round_11/{row['marker']}.tfa").expanduser(), 'fasta')) if _.description == marker]
    if len(records) == 0:
        print(accession, marker)
    SeqIO.write(records, 'query.fasta', 'fasta')

    os.system(f'blastn -query query.fasta -subject {file} -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task blastn -perc_identity 100 >> ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/all_alleles_blastn/blastn_output.txt')
    return accession


MLST_data['FILE'] = MLST_data.apply(run_blastn, axis=1)

MLST_data.to_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/MLST_alleles.csv").expanduser(), index=False)

