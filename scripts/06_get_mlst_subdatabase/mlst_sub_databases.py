#SetUp
import pandas as pd
import re
from pathlib import Path
from shutil import copy
import sys 

#Passing given arguments to variables
current_MLST_output = sys.argv[1]
current_scheme = sys.argv[2]
save_subdatabase_to = sys.argv[3]


# Load MLST output
df = pd.read_csv(f"{current_MLST_output}",sep='\t')

# Making column with accession numbers
def get_accession(row):
    """Return accession number from 'FILE' column."""

    try:
        acc = '_'.join(row['FILE'].split('/')[-1].split('_')[0:2])
    except IndexError:
        acc = 'NA'

    return acc

df["accession"] = df.apply(get_accession, axis=1)

stats = pd.DataFrame()

alleles = ['16S', 'rpoB', 'atpD', 'gyrB', 'recA', 'trpB']
#How many partial sequences have we found for each allele?
for _ in alleles:
    stats[f'partial_{_}'] = df[f'{_}'].str.count('\?')



#How many partial sequences in total?
total_partial = stats['partial_atpD'].sum() + stats['partial_rpoB'].sum() + stats['partial_16S'].sum() + stats['partial_recA'].sum() + stats['partial_gyrB'].sum() + stats['partial_trpB'].sum()
print(f'Total number of partial alleles: {total_partial}')

#How many missing alleles in total?
for _ in alleles:
    stats[f'missing_{_}'] = df[f'{_}'].str.count('\-')


total_missing = stats['missing_atpD'].sum() + stats['missing_rpoB'].sum() + stats['missing_16S'].sum() + stats['missing_recA'].sum() + stats['missing_gyrB'].sum() + stats['missing_trpB'].sum()
print(f'Total number of missing alleles: {total_missing}')



# Partial alleles
def partial_alleles(row, col):
    """Retrun alleles with partial alleles (42?) replaces with
    accession number."""

    try:
        partial = row[col].replace(row[col], "-") if "?" in row[col] else row[col]
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
    accession number for specified df col"""
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


#Drop multiple alleles
def drop_multiple_allelles(df, col):
    """Drop rows if multiple diffrent alleles present"""

    data_points_val = len(df)
    df = df[~df[col].str.contains(",")]
    print(f'Number of genomes with multiple non identical {col} copies: {data_points_val - len(df)}')

    return df
df = drop_multiple_allelles(df=df, col='trpB_new')
df = drop_multiple_allelles(df=df, col='rpoB_new')
df = drop_multiple_allelles(df=df, col='recA_new')
df = drop_multiple_allelles(df=df, col='gyrB_new')
df = drop_multiple_allelles(df=df, col='atpD_new')
df = drop_multiple_allelles(df=df, col='16S_new')


#Load ST profiles (output from mlst_analysis.py)
pubMLST = pd.read_csv(f"{current_scheme}",sep='\t', dtype="object")
pubMLST.columns = ['ST', '16S_new', 'atpD_new', 'gyrB_new', 'recA_new', 'rpoB_new', 'trpB_new']

#Assign ST values to genomes 
df = pd.merge(df, pubMLST, how='left', on=['16S_new', 'atpD_new', 'gyrB_new', 'recA_new', 'rpoB_new', 'trpB_new'])
df = df.drop(['ST_x'], axis=1)
#rename columns
df.columns = ['ST' if x=='ST_y' else x for x in df.columns]



#Change NaN to 0 int
df = df.fillna(0)
df = df.astype({"ST": int})


# #Extract rows where there is no ST value
df.drop(df[df['ST'] > 1].index, inplace = True)

#Get a list of genomes where there was no ST value assigned
new_database =df['accession'].tolist()


datadir = Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/mlst_databases_for_scheme_extension/all_genomes_for_mlst_extension").expanduser()
filenames = sorted(datadir.glob("GCF*genomic.gbff"))


x = []
#Extract genbank files for these genomes and run MLST again
for _ in new_database:
    for filename in filenames:
        if _ in str(filename):
            x.append(_)
        #     copy(filename, Path(f"{save_subdatabase_to}"))


print(len(set(x)))

data = ['_'.join(_.split('/')[-1].split('_')[0:2]) for _ in pd.read_csv("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/strep_mlst2.txt", sep='\t')['FILE']]

for _ in x:
    if _ not in data:
        print(_)

