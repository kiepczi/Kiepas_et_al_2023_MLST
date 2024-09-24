#Set up
import pandas as pd
from pathlib import Path
import re
from Bio import SeqIO
from collections import defaultdict
import numpy as np
from shutil import copy
import os

#Load data MLST data
df = pd.read_csv(Path("../../supplementary_file_4/output/genome_filtration_run/mlst_pubmlst_strep_mincov80.txt"),sep='\t')

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
        partial = row[col].replace(row[col], row['accession']) if "~" in row[col] else row[col]
    except IndexError:
        partial = 'NA'

    return partial

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

#Copy df to include only acession column and all 'new' alleles
df2 = df[['accession', '16S_new', 'atpD_new', 'recA_new', 'rpoB_new', 'gyrB_new', 'trpB_new']].copy()

#When allele was found by mlst return + otherwise -
def is_allele_found(row, col):
    """Return '+' when allele already found otherwise '-'. """
    try:
        partial = row[col].replace(row[col], '+') if "GCF" not in row[col] else '-'
    except IndexError:
        partial = 'NA'

    return partial


df2['atpD_new'] = df2.apply(is_allele_found, col='atpD_new', axis=1)
df2['gyrB_new'] = df2.apply(is_allele_found, col='gyrB_new', axis=1)
df2['recA_new'] = df2.apply(is_allele_found, col='recA_new', axis=1)
df2['rpoB_new'] = df2.apply(is_allele_found, col='rpoB_new', axis=1)
df2['trpB_new'] = df2.apply(is_allele_found, col='trpB_new', axis=1)
df2['16S_new'] = df2.apply(is_allele_found, col='16S_new', axis=1)


#Load data from blastn and tblastn searches
header_list = ['query_accession', 'subject_accesion', 'identity_pc', 'alignment_length', 'mismatches', 'gap_opens', 'query_start', 'query_end', 'subject_start', 'subject_end', 'e_value', 'bit_score', 'qcovs']
df_16S = pd.read_csv(Path("../output/BLAST_output/blastn_16S.txt"),sep='\t', names=header_list)
df_atpD = pd.read_csv(Path("../output/BLAST_output/tblastn_atpD.txt"),sep='\t', names=header_list)
df_recA = pd.read_csv(Path("../output/BLAST_output/tblastn_recA.txt"),sep='\t', names=header_list)
df_rpoB = pd.read_csv(Path("../output/BLAST_output/tblastn_rpoB.txt"),sep='\t', names=header_list)
df_gyrB = pd.read_csv(Path("../output/BLAST_output/tblastn_gyrB.txt"),sep='\t', names=header_list)
df_trpB = pd.read_csv(Path("../output/BLAST_output/tblastn_trpB.txt"),sep='\t', names=header_list)


#check if novel allele should be found by mlst tool.
#This is when minimum identity is 80 and minimum coverage is 80
def minid80_mincov80(row):
    """Return '+' when minimum idetity and coverage was 80%. """

    try:
        allele = '+' if row['identity_pc'] >80 and row['qcovs'] >80 else '-'
    except IndexError:
        allele = 'NA'

    return allele

df_16S["novel_allele"] = df_16S.apply(minid80_mincov80, axis=1)
df_atpD["novel_allele"] = df_atpD.apply(minid80_mincov80, axis=1)
df_recA["novel_allele"] = df_recA.apply(minid80_mincov80, axis=1)
df_rpoB["novel_allele"] = df_rpoB.apply(minid80_mincov80, axis=1)
df_gyrB["novel_allele"] = df_gyrB.apply(minid80_mincov80, axis=1)
df_trpB["novel_allele"] = df_trpB.apply(minid80_mincov80, axis=1)



def novel_allele_to_dict(dataframe):
    """Return dictionary with '+' when novel allele should be found by mlst
    and '-', when not, keyed by genome accession number
    """
    lists = [_ for _ in dataframe['novel_allele']]
    keys = [_.split('/')[0] for _ in dataframe['subject_accesion']]


    dicts = defaultdict(list)
    for k, v in zip(keys, lists):
        dicts[k].append(v)

    for k, v in dicts.items():
        if '+' in v:
            dicts[k] = '+'
        else:
            dicts[k] = '-'

    return dicts

rRNA = novel_allele_to_dict(df_16S)
atpD = novel_allele_to_dict(df_atpD)
recA = novel_allele_to_dict(df_recA)
rpoB = novel_allele_to_dict(df_rpoB)
gyrB = novel_allele_to_dict(df_gyrB)
trpB = novel_allele_to_dict(df_trpB)


df2['16S_new'] = df2['accession'].map(dict(rRNA)).fillna(df2['16S_new'])
df2['atpD_new'] = df2['accession'].map(dict(atpD)).fillna(df2['atpD_new'])
df2['recA_new'] = df2['accession'].map(dict(recA)).fillna(df2['recA_new'])
df2['rpoB_new'] = df2['accession'].map(dict(rpoB)).fillna(df2['rpoB_new'])
df2['gyrB_new'] = df2['accession'].map(dict(gyrB)).fillna(df2['gyrB_new'])
df2['trpB_new'] = df2['accession'].map(dict(trpB)).fillna(df2['trpB_new'])

#Drop genomes when at least one allele cannot be found
df2 = df2.replace('-', np.nan)
df2 = df2.dropna(axis=0)

#Getting genome acession numbers 
genomes = [_ for _ in df2['accession']]



datadir = Path("../../supplementary_file_2/data/streptomyces_genomes")
filenames = sorted(datadir.glob("GCF*/GCF*genomic.gbff"))



if not os.path.exists(Path('../output/all_genomes_for_mlst_extension')):
    os.makedirs(Path('../output/all_genomes_for_mlst_extension').expanduser())

blastn_trpB = []
for _ in genomes:
    for filename in filenames:
        if _ in str(filename):
            copy(filename, Path('../output/all_genomes_for_mlst_extension').expanduser())



