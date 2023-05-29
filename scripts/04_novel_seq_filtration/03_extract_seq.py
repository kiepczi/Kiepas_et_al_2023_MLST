
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pathlib import Path


#Extracting sequences of interest and writing them into fasta file
def get_sequences(gb_file, gene_of_interst, max_start, max_end):
    """Return nt sequences for gene_of_intrest. 
     
    :param gb_file: NCBI GenBank File 
    :param gene_of_inetrest: str of gene for which the protein sequence 
    will  be extracted
    """
    gbank = SeqIO.parse(gb_file, "genbank")

    prot_seq = []
    for genome in gbank:
        for gene in genome.features:
             if 'gene' in gene.qualifiers.keys():
                if gene_of_interst in gene.qualifiers['gene'][0] and 'product' in gene.qualifiers:
                    start = gene.location.nofuzzy_start
                    end = gene.location.nofuzzy_end
                    if start <= max_start and end >= max_end:
                    
                        if gene.strand == -1:   #if the gene is on reverse strand, get a reverse_complement()
                            seqs = genome.seq[start:end].reverse_complement()
                        else:
                            seqs = genome.seq[start:end]#if the gene is on forward strand, reverse_complement() is not needed
                        record = SeqRecord(
                            seqs,
                            id=genome.id,
                            name=genome.name,
                            description= f'{gene_of_interst}_{genome.description}',
                            
                        )
    
                    prot_seq.append(record)

        SeqIO.write(prot_seq, Path(f"~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/whole_{gene_of_interst}.fasta").expanduser(), "fasta")
    return prot_seq
    
atpD = get_sequences(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/streptomyces_genomes/GCF_009862685.1/GCF_009862685.1_ASM986268v1_genomic.gbff").expanduser(), "atpD", 22971, 23564)
trpB = get_sequences(Path("~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/streptomyces_genomes/GCF_016860545.1/GCF_016860545.1_ASM1686054v1_genomic.gbff").expanduser(), "trpB", 31746, 32312)