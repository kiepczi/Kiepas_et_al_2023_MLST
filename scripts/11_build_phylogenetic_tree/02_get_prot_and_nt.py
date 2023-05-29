#SetUp
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from random import sample
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
import sys



#Loading data containing status information about each marker allele (eg. location, preudo genes etc. generated during the 07_revision_scheme step)
df = pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/MLST_alleles_status.csv").expanduser(), dtype=object)

#Getting a list of genomes with either known or novel STs, as only these will be included in the phylogenetic analysis
genomes_of_interest = [_ for _ in pd.read_csv(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_revision/Genomes_ST_info.csv").expanduser(), dtype=object)['accession']]

#Selecting rows from df1 that contain information about the genomes of interest
df = df.loc[df['assembly_accession'].isin(genomes_of_interest)]

def get_gene(row):
    """Return gene name from the query column."""

    try:
        acc = row['query'].split('_')[0]
    except IndexError:
        acc = 'NA'

    return acc

df["gene"] = df.apply(get_gene, axis=1)





#Retain columns of intrest (eg. assembly_accession, gene, and product_accession)
for _ in df:
    if _ not in ['assembly_accession', 'gene', 'product_accession']:
        del df[_]

#Reshape dataframe so each column has product_accession for one of the 5 markers
df = df.pivot_table(index='assembly_accession', columns='gene', values='product_accession', aggfunc='first')

#Get function to trim stop codon sequences
def trim_stop_codons(sequence):
    """Trim sequences to remove stop codons.
    """
    stop_codons = ["TAG", "TAA", "TGA"]


    if sequence[-3:] in stop_codons:
         sequence = sequence[:-3]
    else:
        sequence = sequence
    
    return sequence


def trim_seq(prot, nt):

    """Trim protein and nucleotide sequences,
    so that the length of nucleotide sequence
    is exactly the 3x the length of the protein
    """
    if len(nt) - len(prot*3) < 0:
        prot_idx = len(nt) - len(prot*3)
        prot = prot[:prot_idx]
        nt_idx = len(prot*3) - len(nt)
        nt = nt[:nt_idx]
    elif len(nt) - len(prot*3) >= 1:
        nt_idx = len(prot*3)- len(nt)
        nt = nt[:nt_idx]
    else:
        nt = nt
        prot = prot

    return prot, nt




def extract_sequences(gb_file, protein_id):
    """Extract protein and nucleotide sequences
    if the protein_id found in gene.qualifiers['prot_id]
    """
    gb = SeqIO.parse(gb_file, 'genbank')
    accession = str(gb_file).split('/')[-2]





    for genome in gb:
            for gene in genome.features:
                if gene.type=="CDS":
                    if 'protein_id' in gene.qualifiers:
                        if protein_id in gene.qualifiers['protein_id'][0]:
                            #Getting portein sequence
                            prot_seqs = Seq(gene.qualifiers['translation'][0])
                            #Getting protein sequences

                            #Getting nucleotide_sequence
                            start = gene.location.nofuzzy_start
                            end = gene.location.nofuzzy_end
                            if gene.strand == -1:   #if the gene is on reverse strand, get a reverse_complement()
                                seqs = genome.seq[start:end].reverse_complement()
                                nt_seqs = trim_stop_codons(seqs)

                            else:
                                seqs = genome.seq[start:end] #if the gene is on forward strand, reverse_complement() is not needed
                                nt_seqs = trim_stop_codons(seqs)

    


    prot_seq, nt_seq = trim_seq(prot_seqs, nt_seqs)

    prot_record = SeqRecord(
                    prot_seq,
                    id=accession,
                    name=accession,
                    description=accession,
                            
                )
    record = SeqRecord(
                    nt_seq,
                    id=accession,
                    name=accession,
                    description=accession,
                            
                )


    return prot_record, record


#Path to genbank files
datadir = Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/").expanduser()
filenames = sorted(datadir.glob("strep*/GCF*/GCF*genomic.gbff"))

genomes = ["GCF_002154385.1"]

#Get dictionaries for each marker gene with protein_id keyed by genome accession
atpD = df.to_dict()['atpD']
gyrB = df.to_dict()['gyrB']
recA = df.to_dict()['recA']
rpoB = df.to_dict()['rpoB']
trpB = df.to_dict()['trpB']

#Hold lists to which new records will be appened
atpD_prot = []
atpD_nt = []

gyrB_prot = []
gyrB_nt = []

recA_prot = []
recA_nt = []

rpoB_prot = []
rpoB_nt = []

trpB_prot = []
trpB_nt = []

for fname in filenames:
    genome = str(fname).split('/')[-2]
    if genome in genomes_of_interest:
        print(genome)
        #atpD
        atpD_prot_seq, atpD_nt_seq = extract_sequences(fname, atpD[genome])
        atpD_prot.append(atpD_prot_seq)
        atpD_nt.append(atpD_nt_seq)
        #gyrB
        gyrB_prot_seq, gyrB_nt_seq = extract_sequences(fname, gyrB[genome])
        gyrB_prot.append(gyrB_prot_seq)
        gyrB_nt.append(gyrB_nt_seq)
        #recA
        recA_prot_seq, recA_nt_seq = extract_sequences(fname, recA[genome])
        recA_prot.append(recA_prot_seq)
        recA_nt.append(recA_nt_seq)
        #rpoB
        rpoB_prot_seq, rpoB_nt_seq = extract_sequences(fname, rpoB[genome])
        rpoB_prot.append(rpoB_prot_seq)
        rpoB_nt.append(rpoB_nt_seq)
        #trpB
        trpB_prot_seq, trpB_nt_seq = extract_sequences(fname, trpB[genome])
        trpB_prot.append(trpB_prot_seq)
        trpB_nt.append(trpB_nt_seq)



SeqIO.write(atpD_prot, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/protein/atpD_prot.fasta").expanduser(), "fasta")
SeqIO.write(atpD_nt, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/nucleotide/atpD_nt.fasta").expanduser(), "fasta")

SeqIO.write(gyrB_prot, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/protein/gyrB_prot.fasta").expanduser(), "fasta")
SeqIO.write(gyrB_nt, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/nucleotide/gyrB_nt.fasta").expanduser(), "fasta")


SeqIO.write(recA_prot, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/protein/recA_prot.fasta").expanduser(), "fasta")
SeqIO.write(recA_nt, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/nucleotide/recA_nt.fasta").expanduser(), "fasta")

SeqIO.write(rpoB_prot, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/protein/rpoB_prot.fasta").expanduser(), "fasta")
SeqIO.write(rpoB_nt, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/nucleotide/rpoB_nt.fasta").expanduser(), "fasta")

SeqIO.write(trpB_prot, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/protein/trpB_prot.fasta").expanduser(), "fasta")
SeqIO.write(trpB_nt, Path("~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/nucleotide/trpB_nt.fasta").expanduser(), "fasta")