"""It was brought to our attention that some of the genomes were suppressed from the NCBI database,
or some of the novel alleles identified by MLST tool do not have a coding sequences,
This script was used to ensure that no novell allele markers or ST are being incorporated 
into the extended streptomyces scheme if there is no other genome that would represent them.
"""
#Set Up
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import re
import shutil

#Loading data with allele status (THIS DATA already excluded suppressed genomes)
data = pd.read_csv("MLST_alleles_status.csv")
data = data[~data.duplicated(['query', 'assembly_accession'])]

#Generate lists of alleles that have and do not have a correspoding protein sequence
alleles_with_protein_seq = [_ for _ in data[data.status == 'with_protein']['query']]
false_novel = [_ for _ in data[data.status != 'with_protein']['query'] if _ not in alleles_with_protein_seq]


#Revising 16S alleles (here we remove novel alleles that have been reported from suppressed genomes)
#Loading the most recent MLST output, and adding file column 
MLST_output = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/strep_mlst11.txt").expanduser(), sep='\t')
MLST_replaced_genomes = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/mlst_replaced_genomes.txt").expanduser(), sep='\t')

MLST_data = pd.concat([MLST_output, MLST_replaced_genomes])
MLST_data['accession'] = MLST_data['FILE'].map({_:'_'.join(_.split('/')[1].split('_')[:2]) for _ in MLST_data['FILE']})
MLST_data['status'] = MLST_data['accession'].map(pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/MLST_genome_status.csv").expanduser()).set_index('accession').to_dict()['status']).fillna('available')

genomes_of_interest = [_ for _ in MLST_data[MLST_data.status == 'available']['accession']]


true_novel_16S = list(set([]))
for _ in MLST_data[MLST_data.status == 'available']['16S']:
    if '?' in _ or '-' in _ or '~' in _:
        pass
    else:
        true_novel_16S.extend(['16S_'+_ for _ in _.split(',')])



for _ in MLST_data[MLST_data.status != 'available']['16S']:
    if '?' in _ or '-' in _ or '~' in _:
        pass
    else:
        for _ in ['16S_'+_ for _ in _.split(',')]:
            if _ not in true_novel_16S:
                false_novel.append(_)



def revise_alleles(allele):
    """Check the alleles, and return new set
    by discrading false positive novel alleles."""

    #Current alleles
    current_alleles = list(SeqIO.parse(Path(f"~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/schemes/strep_mlst_round_11//{allele}.tfa").expanduser(), "fasta"))
    
    #pubMLST_alleles
    pubMLST_alleles = list(SeqIO.parse(Path(f"~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/pubMLST_streptomyces_scheme//{allele}.tfa").expanduser(), "fasta"))
    
    #Hold a list for revised_alleles
    true_novel_allele_seq = []
    old_vs_new_allele = {}
        
    current = 0
    for record in current_alleles:
        if int(record.description.split('_')[-1]) <= len(pubMLST_alleles):
            true_novel_allele_seq.append(record)
            current +=1
            old_vs_new_allele[int(current)] = int(current)

        elif int(record.description.split('_')[-1]) > len(pubMLST_alleles) and record.description not in false_novel:
            current +=1
            old_vs_new_allele[int(record.description.split('_')[-1])] = int(current)
            record.id = allele +'_'+str(current)
            record.name = allele +'_'+str(current)
            record.description = allele +'_'+str(current)
            true_novel_allele_seq.append(record)


    SeqIO.write(true_novel_allele_seq, Path(f'~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/revised_scheme/{allele}.tfa').expanduser(), 'fasta')
    return old_vs_new_allele


#Revisng scheme (Here we remove novel STs that have been reported from supprssed genomes, of genomes that contain at least one allele marker as a pseudo gene)
#Load current scheme, and update the the alleles
scheme = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/schemes/strep_mlst_round_11/strep_mlst_round_11.txt").expanduser(), sep='\t')
alleles = ['16S', 'atpD', 'recA', 'rpoB', 'trpB', 'gyrB']
for _ in alleles:
    x = revise_alleles(_)
    scheme[_] = scheme[_].map(x).fillna(0).astype(int)


#First we can drop the ST, that miss at least one allele
scheme = scheme[(scheme != 0).all(1)]


#Now we can get a list of ST to keep
MLST_data = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/organism_info/strep_mlst_data_16S_round_11.csv").expanduser())
MLST_data['status'] = MLST_data['accession'].map(pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/MLST_genome_status.csv").expanduser()).set_index('accession').to_dict()['status']).fillna('available')
genomes_without_all_coding_seq = [_ for _ in data[data.status != 'with_protein']['assembly_accession']]
genomes_with_ST = [_ for _ in MLST_data[MLST_data.status == 'available']['accession'] if _ not in genomes_without_all_coding_seq]
genomes_with_ST.append('GCF_003947265.1')

found_MLST_data = MLST_data[MLST_data['accession'].isin(genomes_with_ST)]
found_st = sorted(list(set([_ for _ in found_MLST_data['ST']])))

pubMLST_not_found = [_ for _ in range(1, 237) if _ not in found_st]
all_true_st = found_st + pubMLST_not_found

old_vs_new_ST = {}

current = 0
for _ in sorted(all_true_st):
    if _ <= 236:
        current +=1
        old_vs_new_ST[_] = _
    else:
        current +=1 
        old_vs_new_ST[_] = current


#drop where the ST is no longer found
scheme['ST'] = scheme['ST'].map(old_vs_new_ST).fillna(0).astype(int)
scheme = scheme[(scheme != 0).all(1)]
scheme.to_csv(Path('~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/revised_scheme/revised_scheme.txt').expanduser(), sep='\t', index=False)


#Getting local database to re run-mlst with new scheme
datadir = Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/").expanduser()
filenames = sorted(datadir.glob("strep*/GCF*/*.gbff"))
for _ in filenames:
    if str(_).split('/')[-2] in genomes_of_interest and 'from' not in str(_):
        shutil.copy(_, Path(f"~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/database_for_revised_scheme").expanduser())
