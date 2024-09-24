"""This scripts was used to generate dataframe for annotations of the MST.

Here, for each genome with an indetified known or novel STs we will get information like:

- Strain
- Taxonomic name status according to LSPN
- Taxonomic type strain according to LSPN
"""


#Set Up
import pandas as pd
from pathlib import Path
from collections import defaultdict


#Load data (here we load most current data with annotation info)
df = pd.read_csv(Path("../../supplementary_file_6/output/revised_scheme_genome_ST_info.csv"), dtype=object)

#Path to NCBI report files
datadir = Path("../../supplementary_file_2/data/")
filenames = sorted(datadir.glob("strep*/GCF*/*report.txt"))

#Step 1: Getting stain information for each genome
def extract_species_and_refseq(file):

    """Return dictionary with organism name as a value
    keyed by accession number from NCBI assembly report"""

    with open(file) as f:
        for line in f.readlines():
            if "Infraspecific name" in line:
                strain = line.split('strain=')[-1].strip()
            elif 'Streptomyces sp.' in line:
                strain = 'NA'
            elif '# Isolate: ' in line:
                strain = line.split(' ')[-1].strip()


    return strain


strain_dict = {}
for filename in filenames:
    accession = str(filename).split('/')[-2]
    strain = extract_species_and_refseq(filename)
    strain_dict[accession] = strain


df["strain"] = df["accession"].map(strain_dict)


#Step 2: Get assembly level for each genome
def get_assembly_level(file):

    """Return dictionary with organism name as a value
    keyed by accession number from NCBI assembly report"""

    with open(file) as f:
        for line in f.readlines():
            if "Assembly level" in line:
                assembly_level = line.split('level: ')[-1].strip()


    return assembly_level

assembly_status = {}
for filename in filenames:
    accession = str(filename).split('/')[-2]
    status = get_assembly_level(filename)
    assembly_status[accession] = status


df["assembly_status"] = df["accession"].map(assembly_status)


#Step 3: Get information about the assigned NCBI names
#Adding LSPN tax status and type strain info 
import csv
class Taxon:
    
    genus = ""
    species = ""
    subspecies = ""
    status_kind = ""
    status_proposal = ""
    status_nom = ""
    status_tax = ""
    type_strains = []
    
    def has_type_strain(self, strain):
        """Returns true if passed strain is a type strain"""
        if strain.lower() in [_.lower() for _ in self.type_strains]:
            return True
        
        return False
    
    def __init__(self):
        pass
    
    def __repr__(self):
        """This is called when we represent the object as a string"""
        return f"<Taxon: {self.genus}, {self.species}, {self.subspecies}>"
    
    def __str__(self):
        """Return contents of Taxon as a string"""
        return f"{self.genus}, {self.species}, {self.subspecies}"

taxa = []

with open(Path("../input/LSPN_taxonomy.csv")) as ifh:
    reader = csv.reader(ifh, delimiter=",")
    next(reader, None)  # skip header
    for row in reader:
        # Here is where we do the hard work
        taxon = Taxon()  # create new Taxon instance
        
        # Get nomenclature
        taxon.genus = row[0]
        taxon.species = row[1]
        taxon.subspecies = row[2]
        
        # Get status
        # All status is in column 4, but stored as several items
        status = [_.strip() for _ in row[4].split(";")]
        taxon.status_kind = status[0]
        taxon.status_proposal = status[1]
        taxon.status_nom = status[2]
        taxon.status_tax = status[3]
        
        # Get type strains
        taxon.type_strains = row[8].split("; ")
        
        # Add taxon to list
        taxa.append(taxon)

#Extracting items that only belong to genus Strep
strep_lpsn = [_ for _ in taxa if _.genus == 'Streptomyces']
#Creating dictionary with taxon name/species as key and status tax as values
strep_status_tax = {}
for _ in strep_lpsn:
    strep_status_tax[f'{_.genus} {_.species}'] = _.status_tax
{v for k, v in strep_status_tax.items()}

#Map dictionary to our dataframe and create a column `status_tax`
df["status_tax"] = df["organism"].map(strep_status_tax)
df = df.fillna('-')

#Change status_tax to `no species info` if species are *Streptimyces sp.*, and `invalid name` if species name not found in LPSN data
def is_valid_info(row, col, col2):
    try:
        validity = 'not validly described' if 'Streptomyces sp.' in row[col] else 'no record in LPSN entry' if row[col2] == '-' else row[col2]
    except IndexError:
        validity = 'NA'
        
        
    return validity
df['status_tax'] = df.apply(is_valid_info, col='organism', col2='status_tax', axis=1)

#Step 3: Adding strain type info
all_type_strains = []
for _ in strep_lpsn:
    all_type_strains.extend(_.type_strains)

#Get dictionaries with strain name as value, keyed by ST for each row
def ST_strain_to_dict(row):
    """Return a dictionary with strain epithet as a value,
    keyed by ST for each row.
    """

    try:
        row_dict = {}
        key = row['ST']
        value = row['strain']
        row_dict[key] = value
    except IndexError:
        row_dict = 'NA'

    return row_dict
    

st_strains = df.apply(ST_strain_to_dict, axis=1)


#Create a dictionary with ALL organism names as a value, keyed by ST
all_strains_per_st = defaultdict(list)

for k, v in st_strains.items():
    for k1, v1 in v.items():
        all_strains_per_st[int(k1)].append(v1)
LPS_type_strain = defaultdict(list)
for k, v in all_strains_per_st.items():
    for _ in v:
        if _ in all_type_strains:
            LPS_type_strain[k].append('Yes')
        else:
            LPS_type_strain[k].append('No')

def add_type_strain_info(row, col):
    try:
        if row[col] in all_type_strains:
            type_strain = 'Yes'
        else:
            type_strain = 'No'
    except IndexError:
        type_strain = 'NA'
        
        
    return type_strain
df['Type_Strain'] = df.apply(add_type_strain_info, col='strain', axis=1)

df.to_csv(Path("../../supplementary_file_8/output/Genome_ST_info.csv").expanduser(), index=False)