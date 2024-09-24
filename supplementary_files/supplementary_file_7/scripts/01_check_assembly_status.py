"""It was brought to our attention that some of the genomes were suppressed from the NCBI database.
This script was used to identify which genomes were suppresed and deal.
"""

#Set Up
import pandas as pd
from pathlib import Path

#Load data 
NCBI_status = pd.read_csv(Path("../input/NCBI_status/assembly_summary_refseq_historical.txt"), sep='\t', skiprows=1)

#Create a dictionary with RefSeq accession ids as keys, and their status as values
genome_status = NCBI_status.set_index('#assembly_accession').to_dict()['version_status']

#Create a dictionary with RefSeq accession ids as keys, and their date of supression
genome_date = NCBI_status.set_index('#assembly_accession').to_dict()['asm_not_live_date']

#Create a dictionary with RefSeq accession ids as keys, and reason for supression
suppression_reason = NCBI_status.set_index('#assembly_accession').to_dict()['excluded_from_refseq']


#List of all genomes used in the MLST analysis (1938)
datadir = Path("../../supplementary_file_5/output/all_genomes_for_mlst_extension")
filenames = sorted(datadir.glob("GCF*"))
all_genomes = ['_'.join(str(_).split('/')[-1].split('_')[:2]) for _ in filenames]

#We can now identify all gneomes that have been suppressed or replaced
troublesome_genomes = {genome:status for genome, status in genome_status.items() if genome in all_genomes}

#Create dataframe which will contain status information for all genome used in the analysis and their status
genome_status = pd.DataFrame(all_genomes, columns =['accession'])
genome_status['status'] = genome_status['accession'].map(troublesome_genomes)
genome_status = genome_status.fillna('available')
genome_status['date'] = genome_status['accession'].map(genome_date)
genome_status['reason'] = genome_status['accession'].map(suppression_reason)


#saving data
genome_status.to_csv(Path('../output/MLST_genome_status.csv'), index=False)
