# README.md - mlst_sub_databases

## **Overview**

This repository contains a script used to built a mlst databases for genomes that were partically classified by MLST tool. 

*Partically classified genemes include:*
- Genomes with at least 1 partial allele
- Geomes with at least 1 missing allele
- Genomes with at least 1 novel allele

## **Arguments**
To run the following script, pass the following arguments:

The current set of input folders include:
- `arg1` - Most recent [MLST tool](https://github.com/tseemann/mlst) tab-separated output file provided by `--legacy` parameter and `--scheme` parameter. Headings contain:

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

- `arg2` - Most recent MLST scheme; tab-separated file containing one ST definitaion. 

```bash 
ST	16S	atpD	gyrB	recA	rpoB	trpB	
1	1	1	1	1	1	1	
2	2	2	2	2	2	2	
3	3	3	3	3	3	3	
4	4	4	4	4	4	4	
5	5	5	5	5	5	5	
6	6	6	6	6	6	6	
7	7	7	7	7	7	7	
8	8	6	8	8	8	8	
9	9	8	9	9	9	9	
10	3	9	10	10	10	10
```

- `arg3` -  Path to the outut file

## Example Usage

```bash
python 
