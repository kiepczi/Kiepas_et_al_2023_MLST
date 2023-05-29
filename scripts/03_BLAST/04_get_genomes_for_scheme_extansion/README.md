# **README.md - build_local_database_for_mlst**

## **Overview**

This repository contains a script used to build local database for mlst analysis. 

### **How genomes for MLST analysis were chosen**
Genomic sequenced were added to a local mlst database when:
- all 6 markers with exact matches were found by MLST tool
- blast serach found all 6 markers with minimum coverage of 80% and identity of  80% per genome. Otherwise genomes were discarded. 


## **Files and subfolders**
The current set of input files and subfolders include:
- `~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/mlst_pubmlst_strep_mincov80.txt`: [MLST tool](https://github.com/tseemann/mlst) tab-separated output file provided by `--legacy` parameter and `--scheme` parameter. Headings contain:

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
- `~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/streptomyces_genomes`: File with subfolders containing genbank files for genomes used in first MLST analysis

- `~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/*.txt`: File containing tblastn and blastn outputs provided by `outfmt "6 std qcovs"` parameter

The current set of output subfolder include:

- `~/Desktop/Kiepas_et_al_2023_MLST/data/mlst_databases_for_scheme_extension/all_genomes_for_mlst_extension`: File containing genbank files with genomes eligible for MLST analysis. See `How genomes for MLST analysis were chosen` section for more information. 
