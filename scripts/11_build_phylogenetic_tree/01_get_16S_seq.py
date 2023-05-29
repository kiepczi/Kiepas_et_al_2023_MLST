#SetUp
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path


#Load data with found STs and their accession
df = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/Genomes_ST_info.csv").expanduser(), dtype=object)

#Create dictionary with ST values keyed by accession
accession_ST = df.set_index('accession').to_dict()['ST']



#Load profiles
df2 = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/revised_scheme/revised_scheme.txt").expanduser(), sep='\t' ,dtype=object)

#Create dictionary with 16S allele number keyed by ST
ST_16S = df2.set_index('ST').to_dict()['16S']


#Create dictionary with 16S allele keyed by accession
accession_16S = {}
for k, v in accession_ST.items():
    for k2, v2 in ST_16S.items():
        if v == k2:
            accession_16S[k] = '16S_' + v2

#Load seq
records = list(SeqIO.parse(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/revised_scheme/16S.tfa").expanduser(), "fasta"))

new_recordss = {}
for k, v in accession_16S.items():
    for record in records:
        if record.description == v:
            new_recordss[k] = record.seq

new = []
for k, v in new_recordss.items():
    record = SeqRecord(
        v, 
        id=k,
        name=k,
        description=k)
    new.append(record)



SeqIO.write(new, Path('~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/nucleotide/16S_nucletide_found_STs.fasta').expanduser(), "fasta")