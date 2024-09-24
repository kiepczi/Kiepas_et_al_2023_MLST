import pandas as pd
from collections import defaultdict 


data = pd.read_csv("Genome_ST_info.csv")

organism_dict = defaultdict(list)

# Iterate over the DataFrame and populate the defaultdict
for index, row in data.iterrows():
    organism_dict[row['organism']].append(row['ST'])



distribution_all = {}

distribution_old = {}

distribution_novel = {}

for organism, STs in organism_dict.items():
    distribution_all[organism] = len(set(STs))
    distribution_old[organism] = len(set([_ for _ in STs if _ <=237]))
    distribution_novel[organism] = len(set([_ for _ in STs if _ >=238]))



# Combine dictionaries into a single DataFrame
df = pd.DataFrame({'all': distribution_all, 'existing': distribution_old, 'novel': distribution_novel})

# Reset the index to make 'organism' a column
df = df.reset_index().rename(columns={'index': 'organism'})

# Melt the DataFrame to get the desired format
final_df = df.melt(id_vars='organism', var_name='type', value_name='count')

# Drop rows with NaN values
final_df = final_df.dropna()

# Display the final DataFrame
final_df.to_csv("species_st_distibution.csv", index=False)