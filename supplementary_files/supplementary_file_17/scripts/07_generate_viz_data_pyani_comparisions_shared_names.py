"""This script was used to generate dataframe for visualization of pyani analysis for genomes which share identical 16S rRNA sequences.
"""


#Set Up
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict


def comparision_type(row, col, thereshold):
    """Return inter if the ID% is more than 95%"""

    try:
        val = 'intra' if row[col] >=thereshold else 'inter'
    except IndexError:
        val = 'NA'

    return val


#Path to all files
datadir = Path("../output/pyani_matrices_shared_names/")
filenames = sorted(datadir.glob("*.tab"))

#Create empty df
data = pd.DataFrame(columns=['run_id','genome1', 'genome2', 'coverage', 'identity'])

# #Hold defaultdict, so cluster ids can be assigned per group
# cluster_ids = defaultdict(list)

#Matching coverage and identity values
for i in range(117,197):
    #Loading data
    identity = pd.read_csv(Path(f"../output/pyani_matrices_shared_names/matrix_identity_{i}.tab"), sep='\t').rename(columns={'Unnamed: 0': 'genome1'})
    coverage = pd.read_csv(Path(f"../output/pyani_matrices_shared_names/matrix_coverage_{i}.tab"), sep='\t').rename(columns={'Unnamed: 0': 'genome1'})

    #Melting data
    identity_melt = pd.melt(identity, id_vars=['genome1'], value_vars=[_ for _ in identity if _ != 'genome1'], var_name='genome2', value_name='identity')
    coverage_melt = pd.melt(coverage, id_vars=['genome1'], value_vars=[_ for _ in identity if _ != 'genome1'], var_name='genome2', value_name='coverage')
    #Combine data
    combined = pd.merge(identity_melt, coverage_melt,  how='left', left_on=['genome1', 'genome2'], right_on = ['genome1','genome2'])
    #append to the final dataframe
    data = data.append(combined)
    #adding run_id
    data['run_id'] = data['run_id'].fillna(int(i))


data = data[data['genome1'] != data['genome2']]
data['comparision_type_species'] = data.apply(comparision_type, col='identity', thereshold=0.95, axis=1)
data['comparision_type_genus'] = data.apply(comparision_type, col='coverage', thereshold=0.50, axis=1)

#Getting number of connected components in which each name is present
datadir = Path("../input/shared_NCBI_names/")
filenames = sorted(datadir.glob("*/custom_classes.txt"))

connected_component_count_per_species = {}
for _ in filenames:
    connected_component_count_per_species[str(_).split('/')[-2].replace('_', ' ')] = len(set([_ for _ in pd.read_csv(_, sep='\t', names=['hash', 'genome', 'connected_comp'])['connected_comp']]))


#Getting number of STs represented by each name
datadir = Path("../input/shared_NCBI_names/")
filenames = sorted(datadir.glob("*/custom_labels.txt"))

ST_count_per_species = {}
for _ in filenames:
    ST_count_per_species[str(_).split('/')[-2].replace('_', ' ')] = len(set([_.split('ST ')[-1] for _ in pd.read_csv(_, sep='\t', names=['hash', 'genome', 'label'])['label']]))




#We can also add column with ST for genome 1 and genome 2
def get_ST(row, col):
    """Return accession number from 'FILE' column."""

    try:
        acc = row[col].split('- ST ')[-1].split(':')[0]
    except IndexError:
        acc = 'NA'

    return acc

data["genome1_ST"] = data.apply(get_ST, col='genome1', axis=1)
data["genome2_ST"] = data.apply(get_ST, col='genome2', axis=1)


#We can also add column with ST for genome 1 and genome 2
def get_name(row, col):
    """Return accession number from 'FILE' column."""

    try:
        acc = ' '.join(row[col].split('|')[-1].split(' ')[1:3])
    except IndexError:
        acc = 'NA'

    return acc
data["organism"] = data.apply(get_name, col='genome1', axis=1)

data['connected_components_count'] = data['organism'].map(connected_component_count_per_species)
data['ST_count_per_species'] = data['organism'].map(ST_count_per_species)

#Assign cluster ids gruped by the unique number of taxa present
unique_species_cluster_id = defaultdict(list)
for cluster_id, unique_species in data.set_index('run_id').to_dict()['connected_components_count'].items():
    unique_species_cluster_id[unique_species].append(cluster_id)

data['cluster_id'] = data['run_id'].map({_[1]:_[0] for k, v in unique_species_cluster_id.items() for _ in enumerate(v)})

data.to_csv(Path("../../supplementary_file_20/pyani_coverage_identity_shared_NCBI_names_df.csv").expanduser(), index=False)
