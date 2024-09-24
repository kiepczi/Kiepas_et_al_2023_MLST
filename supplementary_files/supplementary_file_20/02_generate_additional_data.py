"""This script was used to generate additional data necessary for generating supporting figures.
"""

#Set Up
import pandas as pd
from pathlib import Path
from collections import defaultdict
from collections import Counter

#Loading data with each ST information. 
MLST_data = pd.read_csv(Path("../supplementary_file_8/output/Genome_ST_info.csv"))

print(MLST_data)

#Here, we will generate data to shown number of total genomes, and number of unique NCBI species names per ST
genome_ST = MLST_data.set_index('accession').to_dict()['ST']
genome_organism = MLST_data.set_index('accession').to_dict()['organism']

organism_per_ST = defaultdict(list)

for genome, ST in genome_ST.items():
    organism_per_ST[ST].append(genome_organism[genome])


df = pd.DataFrame(columns=('ST', 'Genome_Count', 'Unique_Species_Names_Count'))
for ST, organisms in organism_per_ST.items():
    df.loc[len(df)] = [ST, len(organisms), len(set(organisms))]

df.to_csv("ST_ancillary_info.csv", index=False)


print(Counter(organism_per_ST[241]))

# # Here, we will generate data that will be used to generate a graph showing number of total genomes per ANI species vs total number of unique species names per ANI species


# def get_ANI_groups(df, column1, column2):
#     shared_values = {}

#     for index, row in df.iterrows():
#         value1 = row[column1]
#         value2 = row[column2]

#         if pd.notnull(value1) and pd.notnull(value2):
#             if value2 in shared_values:
#                 shared_values[value2].append(value1)
#             else:
#                 shared_values[value2] = [value1]

#     return shared_values


# ANIm_genome_name_count = get_ANI_groups(MLST_data, 'organism', 'pyani_species_ID')


# df2 = pd.DataFrame(columns=('pyANI_species_ID', 'Genome_Count', 'Unique_names_Count'))

# for k, v in ANIm_genome_name_count.items():
#     df2.loc[len(df2)] = [k, len(v), len(set(v))]

# ANIm_ST_count = get_ANI_groups(MLST_data, 'ST', 'pyani_species_ID')

# ANIm_ST_count = {k:len(set(v)) for k, v in ANIm_ST_count.items()}

# df2['unique_ST_count'] = df2['pyANI_species_ID'].map(ANIm_ST_count)

# df2.to_csv("pyANI_species_count_info.csv", index=False)