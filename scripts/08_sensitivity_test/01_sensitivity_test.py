#Set UP
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns 
from pathlib import Path
import re
import sys


#Assign variable to user specified argument
allele_to_drop = sys.argv[1]
save_matrix_to = sys.argv[2]

### Loading all data
#Data 1: Current STs assigned to genomes
MLST_current_STs = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/Genomes_ST_info.csv").expanduser(), dtype=object)

#Data 2: Most recent scheme
scheme = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/revised_scheme/revised_scheme.txt").expanduser(), dtype=object, sep='\t')
accessions = [_ for _ in MLST_current_STs['accession']]
#Data 3: Most recent MLST analysis output
#Here, we will retain rows of genomes of interest only (eg. those that have an assigned ST)
MLST_current_output = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/strep_mlst_revised_scheme.txt").expanduser(), sep='\t')
MLST_current_output['accession'] = MLST_current_output['FILE'].map({_:'_'.join(_.split('/')[1].split('_')[:2]) for _ in MLST_current_output['FILE']})
MLST_current_output = MLST_current_output.loc[MLST_current_output['accession'].isin(accessions)]

def unique_copies(row, col):
    """Return unique allele copies """
    try:
        rRNA = str(set(int(s) for s in re.findall(r'-?\d+\.?\d*', row[col])) if ',' in row[col] else row[col]).replace('{', '').replace('}', '')
    except IndexError:
        rRNA = 'NA'

    return rRNA

alleles = ['16S', 'atpD', 'recA', 'rpoB', 'trpB', 'gyrB']
for _ in alleles:
    MLST_current_output[f'{_}'] = MLST_current_output.apply(unique_copies, col=_, axis=1)

### Generate dataframe with new theoretical STs after dropping one of the alleles
def generate_df_after_exclude_allele(excluded_allele):

    scheme_data = scheme.drop_duplicates([_ for _ in alleles if _ != excluded_allele], keep='first')

    data_no_allele = pd.merge(MLST_current_output, scheme_data, how='inner', on=[_ for _ in alleles if _ != excluded_allele])
    data_no_allele.columns = ['ST' if x=='ST_y' else x for x in data_no_allele.columns]
    data_no_allele = data_no_allele[['ST', 'accession']]



    return data_no_allele

no_allele = generate_df_after_exclude_allele(f'{allele_to_drop}')




#Get dictionary with accession as keys and a list with all genome accessions with the same ST as a value
def common_genomes(dataframe):
    """Return dictionary with accessions as keys, 
    and a list with all genome accessions with the same ST as a value.
    :param data_path: path to the dataset (str)
    
    Output example:
        'genome1': ['genome1', 'genome2', 'genome5'],
        'genome2': ['genome1', 'genome2', 'genome5'],
        'genome5': ['genome1', 'genome2', 'genome5'],
        'genome3': ['genome3', 'genome4', 'genome6'],
        'genome4': ['genome3', 'genome4', 'genome6'],
        'genome6': ['genome3', 'genome4', 'genome6'],
        'genome7': ['genome7']
    """
    ST_accession = dataframe.groupby('ST')['accession'].agg(list).to_dict()

    accession_all_genomes = {}

    for k, v in ST_accession.items():
        for _ in v:
            accession_all_genomes[_] = v

    
    return accession_all_genomes



no_allele = common_genomes(no_allele)
all_alleles = common_genomes(MLST_current_STs)



#Create a matrix

def create_matrix(data_1):
    """Return matrix with accessions from data_1 as columns
    and accessions from data_2 are idicies. Fill each row[col]
    to be the same as idices. 
    :param data_1: list of accessions for first data
    :para, data_2: list of accessions for second data 2
    """
    
    data = []

    matrix = pd.DataFrame(data, columns=['data_2'])
    matrix = matrix.append(pd.DataFrame(data_1, columns=['data_2']), ignore_index=True)
    matrix['data_2'] = data_1
    for _ in data_1:
        matrix[_] = data_1
    matrix = matrix.set_index('data_2')

    return matrix


all_vs_all = create_matrix(accessions)
all_vs_no_allele = create_matrix(accessions)


#Assign all genomes that share the same ST as the genome in idx
for _ in accessions:
    all_vs_all[_] = all_vs_all[_].map(all_alleles)


for _ in accessions:
    all_vs_no_allele[_] = all_vs_no_allele[_].map(no_allele)


def sensitivity_test(row, col):
    """Return 1 if the genome accession in colum name and indices, 
    are found in row[col], otherwise 0. 
    """

    try:
        if str(col) in row[col]:
            row[col] = 1
        else:
            row[col] = 0

    except IndexError:
        row[col] = 'NA'

    return row[col]

# # ALL vs ALL
for _ in accessions:
    all_vs_all[_] = all_vs_all.apply(sensitivity_test, col=_, axis=1)
for _ in accessions:
    all_vs_all[_] = all_vs_all[_].astype(int)

all_vs_all_matrix = all_vs_all[accessions].to_numpy()

for _ in accessions:
    all_vs_no_allele[_] = all_vs_no_allele.apply(sensitivity_test, col=_, axis=1)

for _ in accessions:
    all_vs_no_allele[_] = all_vs_no_allele[_].astype(int)
all_vs_no_allele_matrix = all_vs_no_allele[accessions].to_numpy()

#Get data. lower triangle for no allele, and upper triangle with the allele
all_vs_all_upper = np.triu(all_vs_all_matrix, k=0)
all_vs_no_allele_lower = np.tril(all_vs_no_allele_matrix, k=-1)

combination = all_vs_all_upper + all_vs_no_allele_lower



def sensitivity_visualisation(data, save_as):
    """Return a heatmap. 
    """

    plt.figure(figsize = (50,50))
    sns.set(font_scale=4)
    ax = sns.heatmap(data, cmap='mako', xticklabels=50, yticklabels=50, cbar=False, mask=(data==0))
    
    plt.xlabel('Genome ID', fontsize = 40) # x-axis label with fontsize 15
    plt.ylabel('Genome ID', fontsize = 40) # y-axis label with fontsize 15
    plt.savefig(save_as)

    return plt.savefig(save_as)


sensitivity_visualisation(combination, f'{save_matrix_to}')