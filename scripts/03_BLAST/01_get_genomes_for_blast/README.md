# README.md - get_genomes_for_blast

# **Overview**
This repository contains a script used to extract genomic sequences for which blastn and tblastn search needs to be carried out. 

###  **How genomes for blast search were identified:**
For alleles with missing alleles (indicated by `-`), partial alleles (indicated by `?`) or novel alleles (indicated by `~`), `.genomic.fna` sequences were used to build local databses. 

## **Files and subfolders**

The curent set of input files include:
- `~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/mlst_pubmlst_strep_mincov80.txt"`: [MLST tool](https://github.com/tseemann/mlst) tab-separated output file provided by `--legacy` parameter and `--scheme` parameter. Headings contain:

    - the filename
    - scheme specified with `--scheme` parameter
    - the ST (sequence type)
    - the allele IDs


```bash 
FILE	SCHEME	ST	16S	atpD	gyrB	recA	rpoB	trpB
bacteria/GCF_000092385.1/GCF_000092385.1_ASM9238v1_genomic.gbff.gz	streptomyces	-	~39	50	~49	~48	~50	~55
bacteria/GCF_000154905.1/GCF_000154905.1_ASM15490v1_genomic.gbff.gz	streptomyces	-	96?	151?	94?	-	31?	-
bacteria/GCF_000154925.1/GCF_000154925.1_ASM15492v1_genomic.gbff.gz	streptomyces	-	98?	126	123	130	125	-
bacteria/GCF_000242715.1/GCF_000242715.1_ASM24271v2_genomic.gbff.gz	streptomyces	-	138?	~157	141?	~102	~97	-
```
- `../../../data/raw_data/NCBI_genomes/streptomyces_genomes`: File with subfolders containing genomic sequences for genomes used in MLST analysis


The current set of output subfolders include:
- `../../data/blastn_tblastn_databases/'`: FASTA files for genomic sequences  of genus *Streptomyces*, that were choosen for blast searches. NOTE: each marker may contain diffrent genomic sequences. 


