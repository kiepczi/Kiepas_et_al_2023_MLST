#Set up
import pandas as pd
from pathlib import Path
import re
from Bio import SeqIO


#Load data
df = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/mlst_pubmlst_strep_mincov80.txt").expanduser(),sep='\t')



#Making a column with accession numbers
def get_accession(row):
    """Return accession number from 'FILE' column."""

    try:
        acc = row['FILE'].split('/')[1]
    except IndexError:
        acc = 'NA'

    return acc

df["accession"] = df.apply(get_accession, axis=1)



#Partial alleles
def partial_alleles(row, col):
    """Retrun alleles with partial alleles (42?) replaces with
    accession number."""

    try:
        partial = row[col].replace(row[col], row['accession']) if "?" in row[col] else row[col] 
    except IndexError:
        partial = 'NA'

    return partial



df['atpD_new'] = df.apply(partial_alleles, col='atpD', axis=1)
df['gyrB_new'] = df.apply(partial_alleles, col='gyrB', axis=1)
df['recA_new'] = df.apply(partial_alleles, col='recA', axis=1)
df['rpoB_new'] = df.apply(partial_alleles, col='rpoB', axis=1)
df['trpB_new'] = df.apply(partial_alleles, col='trpB', axis=1)


#Novel alleles
def novel_alleles(row, col):
    """Retrun alleles where novel alleles (~42) replaced with 
    accession number"""
    try:
        novel = row[col].replace(row[col], row['accession']) if "~" in row[col] else row[col]
    except IndexError:
        novel = 'NA'

    return novel

df['atpD_new'] = df.apply(novel_alleles, col='atpD_new', axis=1)
df['gyrB_new'] = df.apply(novel_alleles, col='gyrB_new', axis=1)
df['recA_new'] = df.apply(novel_alleles, col='recA_new', axis=1)
df['rpoB_new'] = df.apply(novel_alleles, col='rpoB_new', axis=1)
df['trpB_new'] = df.apply(novel_alleles, col='trpB_new', axis=1)

#Missing alleles
def missing_alleles(row,col):
    """Return missing alleles replaces 
    with accession number. """
    try:
        missing = row[col].replace(row[col], row['accession']) if "-" in row[col] else row[col]
    except IndexError:
        missing = 'NA'

    return missing

df['atpD_new'] = df.apply(missing_alleles, col='atpD_new', axis=1)
df['gyrB_new'] = df.apply(missing_alleles, col='gyrB_new', axis=1)
df['recA_new'] = df.apply(missing_alleles, col='recA_new', axis=1)
df['rpoB_new'] = df.apply(missing_alleles, col='rpoB_new', axis=1)
df['trpB_new'] = df.apply(missing_alleles, col='trpB_new', axis=1)


#unique 16S copies
def get_16S(row):
    """Return unique copies of 16S rRNA for column '16S'"""
    try:
        rRNA = str(set(int(s) for s in re.findall(r'-?\d+\.?\d*', row['16S'])) if ',' in row['16S'] else row['16S']).replace('{', '').replace('}', '')
    except IndexError:
        rRNA = 'NA'

    return rRNA

df['16S_new'] = df.apply(get_16S, axis=1)
df['16S_new'] = df.apply(partial_alleles, col='16S_new', axis=1)
df['16S_new'] = df.apply(novel_alleles, col='16S_new', axis=1)
df['16S_new'] = df.apply(missing_alleles, col='16S_new', axis=1)

###CREATING LOCAL DATABASE
def extract_seq(filename):
    """Return list of sequences

    :param filename: FASTA file
    """
    seqlist1 = list(SeqIO.parse(filename, "fasta"))

    
    records = []
    for seq in seqlist1:
        seq.description = f"GCF_{str(filename).split('/')[-1].split('_')[1]}/{seq.description}"
        seq.id = f"{seq.description}/{str(filename).split('/')[-1]}"
        records.append(seq)

    return records

def build_local_database(col):
    """Return FASTA file with sequences for blast serach. 

    :param col: column with missing/partial/novel alleles replaced with genome accessions
    
    """
    datadir = Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/streptomyces_genomes").expanduser()
    filenames = sorted(datadir.glob("GCF*/GCF*genomic.fna"))

    files_of_ineterst = [_ for _ in df[col] if 'GCF' in _] #extract genome accessions from column
    blast_records = [] #Hold a list for sequences of interest
    for _ in files_of_ineterst:
        for filename in filenames:
            if _ in str(filename) and 'from' not in str(filename):
                records = extract_seq(filename) #Extracting list of sequences
                blast_records.extend(records)

    return SeqIO.write(blast_records, Path(f'~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/genomes_for_BLAST/{col}_blast.fasta').expanduser(), "fasta")


build_local_database('trpB_new')
build_local_database('16S_new')
build_local_database('atpD_new')
build_local_database('recA_new')
build_local_database('rpoB_new')
build_local_database('gyrB_new')


