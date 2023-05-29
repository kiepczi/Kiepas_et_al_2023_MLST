#Set up
import pandas as pd
import numpy as np
import re
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import sys

#Assign variable to user specified argument
most_recent_mlst_output = sys.argv[1]
most_recent_mlst_scheme = sys.argv[2]
save_new_scheme_to = sys.argv[3]
save_organism_info_to = sys.argv[4]

#Load data. Mlst (v2.19.0) output
df = pd.read_csv(f'{most_recent_mlst_output}', sep='\t')


#Making a column with accession numbers
def get_accession(row):
    """Return accession number from 'FILE' column."""

    try:
        acc = 'GCF_' + row['FILE'].split('/')[1].split('_')[1]
    except IndexError:
        acc = 'NA'

    return acc

df["accession"] = df.apply(get_accession, axis=1)


datadir = Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/streptomyces_genomes").expanduser()
filenames = sorted(datadir.glob("GCF*/*report.txt"))

#Extract organism names from each file and add this info to df
def extract_species_and_refseq(file):

    """Return dictionary with organism name as a value
    keyed by accession number from NCBI assembly report"""

    with open(file) as f:
        for line in f.readlines():
            if "Organism name" in line:
                organism_name = line.split(" ")[4:6]
                species = " ".join(organism_name)
            if "RefSeq assembly accession" in line:
                refseq_assembly = line.strip("\n").split(" ")[-1]
    return refseq_assembly, species
 
species_dict = {}
for filename in filenames:
    refseq_assembly, species = extract_species_and_refseq(filename)
    species_dict[refseq_assembly] = species

df["organism"] = df["accession"].map(species_dict)


#unique 16S copies
def unique_copies(row, col):
    """Return unique allele copies """
    try:
        rRNA = str(set(int(s) for s in re.findall(r'-?\d+\.?\d*', row[col])) if ',' in row[col] else row[col]).replace('{', '').replace('}', '')
    except IndexError:
        rRNA = 'NA'

    return rRNA

alleles = ['16S', 'atpD', 'recA', 'rpoB', 'trpB', 'gyrB']
for _ in alleles:
    df[f'{_}_new'] = df.apply(unique_copies, col=_, axis=1)


#Partial alleles
def partial_alleles(row, col):
    """Retrun alleles with partial alleles (42?) replaces as '-'
    for specified df col"""

    try:
        partial = row[col].replace(row[col], "-") if "?" in row[col] else row[col] 
    except IndexError:
        partial = 'NA'



    return partial

for _ in alleles:
    df[f'{_}_new'] = df.apply(partial_alleles, col=f'{_}_new', axis=1)




#Novel alleles
def novel_alleles(row, col):
    """Retrun alleles where novel alleles (~42) replaced with 
    accession number for specified df col"""
    try:
        partial = row[col].replace(row[col], row['accession']) if "~" in row[col] else row[col]
    except IndexError:
        partial = 'NA'

    return partial

for _ in alleles:
    df[f'{_}_new'] = df.apply(novel_alleles, col=f'{_}_new', axis=1)


#Number of 16S copies per genome
def count_16S_copies(row):
    """Return count of 16S rRNA copies in ['16S] column"""
    try:
        if '-' in row['16S']:
            copies = 0
        elif '?' in row['16S']:
            copies = 0
        else:
            copies = row['16S'].count(',') + 1
    except IndexError:
        copies = 'NA'
    
    return copies

df['16S_copies'] = df.apply(count_16S_copies, axis=1)


#Drop columns that have more than 2 diffrent copies of 16S
def drop_multiple_allelles(df, col):
    """Drop rows if multiple diffrent alleles present"""

    data_points_val = len(df)
    df = df[~df[col].str.contains(",")]
    print(f'Number of genomes with multiple non identical {col} copies: {data_points_val - len(df)}')

    return df

for _ in alleles:
    df = drop_multiple_allelles(df=df, col=f'{_}_new')


#Replacing missing data with np.nan
df = df.replace('-', np.nan)

#Dropping rows if at least one housekeeping gene is not present
data_points_val = len(df)
df = df.dropna(subset=['16S_new', 'atpD_new', 'gyrB_new', 'recA_new', 'rpoB_new', 'trpB_new'])
print(f'Number of genomes with at least 1 missing allele: {data_points_val - len(df)}')


#Drop rows if novel allele still found! THESE WILL BE ADDED TO MLST AND RUN AGAIN UNTIL NO NOVEL ALLELES FOUND!
#This is only to see how the minimum spanning tree will look like for the current data
def drop_novel_allelles(df, col):
    """Drop rows if multiple diffrent alleles present"""

    data_points_val = len(df)
    df = df[~df[col].str.contains("GCF")]
    print(f'Number of genomes with novel alleles {col} copies: {data_points_val - len(df)}')

    return df

for _ in alleles:
    df = drop_novel_allelles(df=df, col=f'{_}_new')


#Assign existing pubMLST ST
current_scheme = pd.read_csv(f'{most_recent_mlst_scheme}',sep='\t', dtype="object")
if 'clonal_complex' in current_scheme:
    current_scheme = current_scheme.drop(['clonal_complex'], axis=1)


current_scheme.columns = ['ST', '16S_new', 'atpD_new', 'gyrB_new', 'recA_new', 'rpoB_new', 'trpB_new']

df = pd.merge(df, current_scheme, how='left', on=['16S_new', 'atpD_new', 'gyrB_new', 'recA_new', 'rpoB_new', 'trpB_new'])
df = df.drop(['ST_x'], axis=1)
#rename columns
df.columns = ['ST' if x=='ST_y' else x for x in df.columns]



#create pubmlst profiles
#Extracting columns of interest
pubmlst_new_profiles = df[['ST', '16S_new', 'atpD_new', 'gyrB_new', 'recA_new', 'rpoB_new', 'trpB_new']]


# #Rename columns so they match pubMLST format
pubmlst_new_profiles.columns = ['ST', '16S', 'atpD', 'gyrB', 'recA', 'rpoB', 'trpB']

#Change NaN to 0 int
pubmlst_new_profiles = pubmlst_new_profiles.fillna(0)
pubmlst_new_profiles = pubmlst_new_profiles.astype({"ST": int})

#Drop duplicates
pubmlst_new_profiles = pubmlst_new_profiles.drop_duplicates()


#Extract rows where there is no ST value
pubmlst_new_profiles.drop(pubmlst_new_profiles[pubmlst_new_profiles['ST'] >= 1].index, inplace = True)


# #Load existing profiles
latest_scheme = pd.read_csv(f'{most_recent_mlst_scheme}',sep='\t', dtype="object")
if 'clonal_complex' in latest_scheme:
    latest_scheme = latest_scheme.drop(['clonal_complex'], axis=1)


total_new_profiles = len(pubmlst_new_profiles)
total_pubmlst_profiles = len(latest_scheme) + 1
total_profiles = total_new_profiles + total_pubmlst_profiles


#Assigning new st
new_st = []
for _ in range(total_pubmlst_profiles, total_profiles):
    new_st.append(_)


pubmlst_new_profiles['ST'] = new_st


#Merging pubMLST and new ST files together 
merged = pd.concat([latest_scheme, pubmlst_new_profiles])
merged = merged.reset_index(drop=True)
merged.to_csv(f'{save_new_scheme_to}', sep="\t", index=False)


#Load updated_pubMLST_profiles and change column names so new ST can be assigned
merged = pd.read_csv(f'{save_new_scheme_to}',sep='\t', dtype=object)
merged.columns = ['ST', '16S_new', 'atpD_new', 'gyrB_new', 'recA_new', 'rpoB_new', 'trpB_new']


# Creating anciallary data
# Match ST from a pubMLST csv
df = pd.merge(df, merged, how='inner', on=['16S_new', 'atpD_new', 'gyrB_new', 'recA_new', 'rpoB_new', 'trpB_new'])


df = df.drop(['ST_x'], axis=1)
#rename columns
df.columns = ['ST' if x=='ST_y' else x for x in df.columns]


info_data = df[['ST', 'accession', 'organism', '16S_copies']]
info_data = info_data.astype({"ST": int})
info_data = info_data.sort_values(by=['ST'])
info_data.to_csv(f"{save_organism_info_to}", index=False)
