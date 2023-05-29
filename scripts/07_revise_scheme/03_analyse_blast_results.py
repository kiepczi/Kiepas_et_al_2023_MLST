"""This scripts was used to revise the scheme. 

Here, we discared novel alleles (and STs) from the scheme if the novel sequence
is a pseudogene without a coding sequence.
"""

import pandas as pd
from pathlib import Path
from Bio import SeqIO


#Loading blastn output
blast_data = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/all_alleles_blastn/blastn_output.txt").expanduser(), sep='\t', header=None, names=['query', 'subject', 'identity', 'algn_length', 'mismatches', 'gap_opens', 'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit_score', 'q_cov'])
blast_data = blast_data.loc[(blast_data['identity'] == 100) & (blast_data['q_cov'] == 100) & (blast_data['mismatches']==0) & (blast_data['gap_opens'] ==0) & (blast_data['evalue'] == 0)]


print(blast_data)

#path to troublesome genomes; getting a dictionary with genomic accession, and aseembly accession
#Running blastn
datadir = Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/").expanduser()
filenames = sorted(datadir.glob("strep*/GCF_*/*.fna"))

genomics_vs_assembly_accession = {}
for f in filenames:
    if 'from' not in str(f):
        records = SeqIO.parse(f, 'fasta')
        for record in records:
            genomics_vs_assembly_accession[record.description.split(' ')[0]] = '_'.join(str(f).split('/')[-1].split('_')[0:2])

blast_data['assembly_accession'] = blast_data['subject'].map(genomics_vs_assembly_accession)

#path to all feature tables
datadir = Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes").expanduser()
filenames = sorted(datadir.glob("strep*/GCF*/GCF*feature_table.txt"))



#Removing rows that belong to supressed genomes
NCBI_status = [genome for genome, status in pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/MLST_genome_status.csv").expanduser()).set_index('accession').to_dict()['status'].items() if status == 'available']
NCBI_status.extend(['GCF_002879675.3', 'GCF_003268535.2', 'GCF_003947265.2', 'GCF_012184365.2'])

blast_data = blast_data[blast_data['assembly_accession'].isin(NCBI_status)]

def get_accession(row, feature_of_interest):
    """Get info for missing seq"""

    try:
        acc = row['subject']
        s_start = row['s_start']
        s_end = row['s_end']
        feature = ''
        for fname in filenames:
            if str(fname).split('/')[-2] == row['assembly_accession']:
                table_data = pd.read_csv(fname, sep='\t')
                feature = [_ for _ in table_data.loc[(table_data['genomic_accession'] == row['subject']) & (table_data['start'] <= s_start) & (table_data['end']>=s_end) & (table_data['# feature'] == 'CDS')][feature_of_interest]][0]

    except IndexError:
        feature = 'NA'


    return feature


blast_data['name'] = blast_data.apply(get_accession, feature_of_interest='name', axis=1)
blast_data['symbol'] = blast_data.apply(get_accession, feature_of_interest='symbol', axis=1)
blast_data['attributes'] = blast_data.apply(get_accession, feature_of_interest='attributes', axis=1)
blast_data['status'] = blast_data.apply(get_accession, feature_of_interest='class', axis=1)
blast_data['product_accession'] = blast_data.apply(get_accession, feature_of_interest='product_accession', axis=1)


blast_data.to_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/MLST_alleles_status.csv").expanduser(), index=False)




