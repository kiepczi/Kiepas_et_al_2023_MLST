"""It was brought to our attention that some of the genomes were suppressed from the NCBI database,
or some of the novel alleles identified by MLST tool do not have a coding sequences,
This script was used to ensure that no novell allele markers or ST are being incorporated 
into the extended streptomyces scheme if there is no other genome that would represent them.

In this scipt we also check if the allele sequences in the updated scheme were reported in the extened scheme,
and make adequate changes. 
"""


#Set Up
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import re
import shutil




#Revising 16S alleles (here we remove novel alleles that have been reported from suppressed genomes)
#Loading the most recent MLST output, and adding file column 
MLST_output = pd.read_csv(("../../supplementary_file_4/output/extension_scheme_round_7/extension_scheme_round_7.txt"), sep='\t', usecols=['FILE', 'SCHEME', 'ST', '16S', 'atpD', 'gyrB', 'recA', 'rpoB', 'trpB'])
MLST_output['accession'] = MLST_output['FILE'].map({_:'_'.join(_.split('/')[-1].split('_')[:2]) for _ in MLST_output['FILE']})
MLST_replaced_genomes = pd.read_csv(Path("../../supplementary_file_4/output/replaced_genomes/streptomyces_replaced_genomes.txt"), sep='\t', usecols=['FILE', 'SCHEME', 'ST', '16S', 'atpD', 'gyrB', 'recA', 'rpoB', 'trpB'])


MLST_data = pd.concat([MLST_output, MLST_replaced_genomes])
MLST_data['accession'] = MLST_data['FILE'].map({_:'_'.join(_.split('/')[-1].split('_')[:2]) for _ in MLST_data['FILE']})
MLST_data['status'] = MLST_data['accession'].map(pd.read_csv(Path("../../supplementary_file_7/output/MLST_genome_status.csv")).set_index('accession').to_dict()['status']).fillna('available')


genomes_of_interest = [_ for _ in MLST_data[MLST_data.status == 'available']['accession']]



#list of markers:
markers = ['16S', 'atpD', 'recA', 'rpoB', 'trpB', 'gyrB']


def identify_true_novel(gene):
    """retun list of alleles reported from genomes
    which were not supressed, or replaced. 
    """

    true_novel = list(set([]))

    for _ in MLST_data[MLST_data.status == 'available'][gene]:
        if '?' in str(_) or '-' in str(_) or '~' in str(_):
            pass
        else:
            true_novel.extend([f'{gene}_'+str(_) for _ in str(_).split(',')])

    return true_novel


all_true_novel_sequences = []

for _ in markers:
    all_true_novel_sequences.extend(identify_true_novel(_))


def identify_false_novel(gene):

    false_novel = []
    for _ in MLST_data[MLST_data.status != 'available'][gene]:
        if '?' in _ or '-' in _ or '~' in _:
            pass
        else:
            for _ in [f'{gene}_'+_ for _ in _.split(',')]:
                if _ not in all_true_novel_sequences:
                    false_novel.append(_)

    return false_novel

all_false_novel_sequences = []

for _ in markers:
    all_false_novel_sequences.extend(identify_false_novel(_))



def revise_alleles(allele):
    """Check the alleles, and return new set
    by discrading false positive novel alleles."""

    #Current alleles
    current_alleles = list(SeqIO.parse(Path(f"../../supplementary_file_6/output/schemes/scheme_round_6/{allele}.tfa"), "fasta"))
    
    #pubMLST_alleles
    pubMLST_alleles = list(SeqIO.parse(Path(f"../../supplementary_file_3/{allele}.tfa").expanduser(), "fasta"))
    
    #Hold a list for revised_alleles
    true_novel_allele_seq = []
    old_vs_new_allele = {}
        
    current = 0
    for record in current_alleles:
        if int(record.description.split('_')[-1]) <= len(pubMLST_alleles):
            true_novel_allele_seq.append(record)
            current +=1
            old_vs_new_allele[int(current)] = int(current)

        elif int(record.description.split('_')[-1]) > len(pubMLST_alleles) and record.description not in all_false_novel_sequences:
            current +=1
            old_vs_new_allele[int(record.description.split('_')[-1])] = int(current)
            record.id = allele +'_'+str(current)
            record.name = allele +'_'+str(current)
            record.description = allele +'_'+str(current)
            true_novel_allele_seq.append(record)


    SeqIO.write(true_novel_allele_seq, Path(f'../output/revised_scheme/{allele}.tfa'), 'fasta')
    return old_vs_new_allele



#Load and revise the scheme
scheme = pd.read_csv(Path("../../supplementary_file_6/output/schemes/scheme_round_6/scheme_round_6.txt"), sep='\t')

for _ in markers:
    marker_new = revise_alleles(_)
    scheme[_] = scheme[_].map(marker_new).fillna(0).astype(int)

#We can now get the list of current STs to keep.
#To do this we will load the genome_ST_info and retain rows/STs that were assigned for genomes that are still avaliable
current_genome_ST_info = pd.read_csv("../../supplementary_file_6/output/extension_round_7_genome_ST_info.csv")
# Filter the DataFrame based on the 'accession' column. Here we also add 3 genomes, which were replaced, but their updated version still had the same ST
current_genome_ST_info = current_genome_ST_info[current_genome_ST_info['accession'].isin(genomes_of_interest + ['GCF_000612545.1', 'GCF_002224125.1', 'GCF_003947265.1'])]


# Get the unique 'ST' values from the filtered DataFrame
st_to_keep = current_genome_ST_info['ST'].unique().tolist()

#Adding pubMLST STs that are nore represented in NCBI
for _ in range(1,238):
    if _ not in st_to_keep:
        st_to_keep.append(_)


scheme = scheme[scheme['ST'].isin(st_to_keep)]


old_vs_new_ST = {}

current = 0
for _ in sorted([_ for _ in scheme['ST']]):
    if _ <= 237:
        current +=1
        old_vs_new_ST[_] = _
    else:
        current +=1 
        old_vs_new_ST[_] = current


scheme['ST'] = scheme['ST'].map(old_vs_new_ST).fillna(0).astype(int)

print(scheme)



scheme.to_csv(Path('../output/revised_scheme/revised_scheme.txt'), sep='\t', index=False)


# Getting local database to re run-mlst with new scheme
datadir = Path("../../supplementary_file_2/data/")
filenames = sorted(datadir.glob("strep*/GCF*/*.gbff"))
for _ in filenames:
    if str(_).split('/')[-2] in genomes_of_interest and 'from' not in str(_):
        shutil.copy(_, Path(f"../output/database_for_revised_scheme"))

