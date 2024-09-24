#SetUp
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path


#Load data with found STs and their accession
df = pd.read_csv(Path("../../supplementary_file_8/output/Genome_ST_info.csv"), dtype=object)

#Create dictionary with ST values keyed by accession
accession_ST = df.set_index('accession').to_dict()['ST']



#Load profiles
df2 = pd.read_csv(Path("../../supplementary_file_8/output/revised_scheme/revised_scheme.txt"), sep='\t' ,dtype=object)

#Create dictionary with 16S allele number keyed by ST
ST_16S = df2.set_index('ST').to_dict()['16S']


#Create dictionary with 16S allele keyed by accession
accession_16S = {}
for genome, ST in accession_ST.items():
            accession_16S[genome] = '16S_' + ST_16S[ST]



#Load seq
records = list(SeqIO.parse(Path("../../supplementary_file_8/output/revised_scheme/16S.tfa"), "fasta"))

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



SeqIO.write(new, Path('../output/nucleotide/16S_nucletide_found_STs.fasta'), "fasta")