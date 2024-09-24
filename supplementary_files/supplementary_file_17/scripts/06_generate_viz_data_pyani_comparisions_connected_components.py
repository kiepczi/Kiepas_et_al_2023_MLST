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
datadir = Path("../output/pyani_matrices_connected_components/")
filenames = sorted(datadir.glob("*.tab"))

#Create empty df
data = pd.DataFrame(columns=['run_id','genome1', 'genome2', 'unique_taxa_names_per_connected_component', 'coverage', 'identity'])

# #Hold defaultdict, so cluster ids can be assigned per group
# cluster_ids = defaultdict(list)


unique_tax_names = {}
#Matching coverage and identity values
for i in range(1,117):
    #Loading data
    identity = pd.read_csv(Path(f"../output/pyani_matrices_connected_components/matrix_identity_{i}.tab"), sep='\t').rename(columns={'Unnamed: 0': 'genome1'})
    coverage = pd.read_csv(Path(f"../output/pyani_matrices_connected_components/matrix_coverage_{i}.tab"), sep='\t').rename(columns={'Unnamed: 0': 'genome1'})

    #Getting number of unique specie name per connected component
    unique_tax_names[i] = len(set([' '.join(_.split('|')[-1].split(' ')[1:3]) for _ in identity if _ != 'genome1']))
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
data['unique_taxa_names_per_connected_component'] = data['run_id'].map(unique_tax_names)


#Assign cluster ids gruped by the unique number of taxa present
unique_species_cluster_id = defaultdict(list)
for cluster_id, unique_species in data.set_index('run_id').to_dict()['unique_taxa_names_per_connected_component'].items():
    unique_species_cluster_id[unique_species].append(cluster_id)

data['cluster_id'] = data['run_id'].map({_[1]:_[0] for k, v in unique_species_cluster_id.items() for _ in enumerate(v)})


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

data.to_csv(Path("../../supplementary_file_20/pyani_coverage_identity_df.csv").expanduser(), index=False)

##Generating data for comparisions between the same STs

#Selecting rows for intra STs comparions
ST_data = data.loc[data['genome1_ST'] == data['genome2_ST']]


#Getting number of unique_taxa_names_per_ST
MLST_data_info = pd.read_csv(Path("../../supplementary_file_8/output/Genome_ST_info.csv"))

ST_names = {}
organism_name = MLST_data_info.set_index('accession').to_dict()['organism']
organism_ST = MLST_data_info.set_index('accession').to_dict()['ST']
for accession, organism in organism_name.items():
    for accession_2, ST in organism_ST.items():
        if accession == accession_2:
            ST_names[ST] = organism


#Getting number of unique_taxa_names_per_ST
MLST_data_info = pd.read_csv(Path("../../supplementary_file_8/output/Genome_ST_info.csv"))

ST_names = defaultdict(list)
organism_name = MLST_data_info.set_index('accession').to_dict()['organism']
organism_ST = MLST_data_info.set_index('accession').to_dict()['ST']
for accession, organism in organism_name.items():
    for accession_2, ST in organism_ST.items():
        if accession == accession_2:
            ST_names[ST].append(organism)



ST_data['unique_taxa_names_per_ST'] = ST_data['genome2_ST'].map({str(k):len(set(v)) for k, v in ST_names.items()})


#Getting new cluster_id for better visualization
#Assign cluster ids gruped by the unique number of taxa representig a single ST
unique_ST_cluster_id = defaultdict(list)
for cluster_id, unique_species in ST_data.set_index('genome1_ST').to_dict()['unique_taxa_names_per_ST'].items():
    unique_ST_cluster_id[unique_species].append(cluster_id)

ST_data['cluster_id'] = ST_data['genome1_ST'].map({str(int(_[1])):_[0] for k, v in unique_ST_cluster_id.items() for _ in enumerate(v)})

ST_data.to_csv(Path("../../supplementary_file_20/pyani_coverage_identity_ST_comparisions_df.csv").expanduser(), index=False)

