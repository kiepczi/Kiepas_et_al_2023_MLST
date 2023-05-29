# **README.md - mlst_ST_assignment**

## **Overview**

This repository cotains a script used to identify and assign new ST to an already existing pubMLST databse. This repository is also used to produce `csv` file, which contains genome accession number, organism name and 16S copies associated with each ST. 

## **Arguments**
To run the script, pass the following arguments:


- `arg1`: Most recent MLST tool](https://github.com/tseemann/mlst) tab-separated output file provided by `--legacy` parameter and `--scheme` parameter. Headings contain:

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


- `argv2` - Most recent MLST scheme, tab-separated file containing one ST definitaion per row. 

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

- `arg3` -  Path to a file where a new extended scheme will be saved to 

- `arg4` - Path to a file to generate csv output with accession number, organism name and 16S copies associated with each ST

To run this script, path to genbank reports is required, this is to exract the organism names. 

- `bacteria`: File containing subfolders with assembly reports. 

Example useage for the most recent scheme:
```bash
python mlst_strep_analysis.py ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/mlst_pubmlst_strep_mincov80.txt ~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/pubMLST_streptomyces_scheme/streptomyces.txt  ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/schemes/strep_mlst_round_1.txt ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/organism_info/strep_mlst_data_round_1.csv
``