# **README.md** - **Kiepas_et_al_2023_MLST**

This repository contains all supplementary information for analyses describing inconsistencies between taxonomies inferred using MLST and whole-genome identities in *Streptomyces*.

## **MLST and *Streptomyces***
This repository is provided to enable both reproduction and independed exploration of the analysis reported in this manuscript.


# Table of contents

1. [Reporting Problems](#reporting-problems)
2. [Contributors](#contributors)
3. [Contact Us](#contact-us)
4. [Downloading Repository](#downloading-repository)
5. [Set Up](#set-up)
6. [Repository Files](#repository-files)
7. [Reproducing analyses](#reproducing-analyses)



## **Reporting Problems**
Please report any issues and problem with this repository [here](https://github.com/kiepczi/Kiepas_et_al_2023_MLST/issues).

## **Contributors**
This manuscript has the following contributors:
- [Angelika B Kiepas](https://github.com/kiepczi) - PhD Candidate, Univeristy of Strathclyde
- [Dr Leighton Pritchard](https://github.com/widdowquinn) - Strathclyde Chancellor's Fellow, Univeristy of Strathclyde
- [Prof Paul A Hoskisson](https://github.com/PaulHoskisson) - Professor, Univeirsty of Strathclyde

## **Contact Us**
How to reach us:
- Angelika B Kiepas:
    - Email: angelika.kiepas@strath.ac.uk
    - Twitter: [@kiepczi](https://twitter.com/kiepczi?lang=en)
- Leighton Pritchard:
    - Email: leighton.pritchard@strath.ac.uk
    - Twitter: [@widdowquinn](https://twitter.com/widdowquinn)
- Paul A Hoskisson:
    - Email: paul.hoskisson@strath.ac.uk
    - Twitter: [@PaulHoskisson](https://twitter.com/PaulHoskisson?ref_src=twsrc%5Egoogle%7Ctwcamp%5Eserp%7Ctwgr%5Eauthor)


## **Downloading Repository**
If you wish to indepedently explore, repoduce and/or validate the analyses reported in the manuscipt, you can use `git` to clone this reporitory to your machine's Desktop directory. 

```bash
git clone https://github.com/kiepczi/Kiepas_et_al_2023_MLST.git
```
Alternatively, click [here](https://github.com/kiepczi/Kiepas_et_al_2023_MLST/archive/refs/heads/main.zip) to download this repository as a `.zip` file.

To use this repository please enusre that it is located on your local machine's Desktop directory, as we used absolute paths. 


## **Set Up**
We strongly recommend to create a `conda` enviroment specific for this activity, for example using the commands:

```
bash
conda create --name mlst_strep python=3.8 -y
conda activate mlst_strep
conda install --file requirements.txt -y
```

Please also download the following software, and follow the installation instructions are appropriare for each program:

- [pyANI v0.3](https://github.com/widdowquinn/pyani)
- [Cytoscape v3.9.0](https://cytoscape.org)
- [MLST v2.22.1](https://github.com/tseemann/mlst)



## **Repository Files**
Here you can find a list of all supplementary files provided in this repository, and current set of subfolders include:

**Supplementary File 1**: Schematic representation of the pipeline used to extend the pubMLST *Streptomyces* scheme.  (PDF 85KB)

**Supplementary File 2**: Download and store genomes. This ZIP file contains bash script used to download all *Streptomyces*, *Kitasatospora*, *Streptacidiphilus* and *Streptoalloteichus* genomes used in this manuscript. The ZIP file also contains five separate txt files; streptomyces_genomes.txt listing all genome accessions for initial *Streptomyces* candidates, streptomyces_replaced_geonomes.txt file with accessions for replaced *Streptomyces* genomes, kitasatospora_genomes.txt for *Kitasatospora* genomes, streptacidiphilus_genomes.txt for *Streptacidiphilus* genomes and streptoalloteichus_genomes.txt for *Streptoalloteichus* genomes. (ZIP 16KB)

**Supplementary File 3**: The canonical *Streptomyces* pubMLST scheme. This ZIP file contains the canonical *Streptomyces* pubMLST scheme manually downloaded from https://pubmlst.org. (ZIP 63KB)

**Supplementary File 4**: Run MLST tool. This ZIP file contains bash script used to run MLST tool, and all generated output for scheme extension, scheme revision, sister genera analysis and *Streptomyces* replaced genomes. (ZIP 911KB)

**Supplementary File 5**: Run BLAST. This ZIP file contains all Python and bash scripts used to filtrate *Streptomyces* genomes using BLAST analysis, and all generated output. (ZIP 457KB)

**Supplementary File 6**: Extend pubMLST scheme. This ZIP file contains all Python scripts used to extend the pubMLST *Streptomyces* scheme, and all generated outputs. (ZIP 1.6MB) 

**Supplementary File 7**. Check assembly status. This ZIP file contains a Python script, NCBI assembly report, and all generated outputs for checking assembly status. (ZIP 5.1MB)

**Supplementary File 8**. Revise the extended scheme. This ZIP file contains all Python scripts used to revise the extended scheme. Additionally, the ZIP file contains the revised scheme itself, along with a txt file containing genome accessions for the genomes used to extended the scheme with exclusion of the suppressed genomes. (ZIP 291KB)

**Supplementary File 9**. Sensitivity Test. This ZIP file contains Python script used for sensitivity test. (ZIP 3KB)

**Supplementary File 10**. Sensitivity test heatmap for 16S marker gene. Black cells in the heatmap correspond to genomes sharing the same STs. The upper triangular region represents genomes sharing the same STs before the marker was excluded, and the lower triangular region represents genomes which are represented by the same STs after the exclusion of the 16S gene from the scheme. (PDF 9.5MB)


**Supplementary File 11**. Sensitivity test heatmap for atpD marker gene. Black cells in the heatmap correspond to genomes sharing the same STs. The upper triangular region represents genomes sharing the same STs before the marker was excluded, and the lower triangular region represents genomes which are represented by the same STs after the exclusion of the atpD gene from the scheme. (PDF 9.5MB)

**Supplementary File 12**. Sensitivity test heatmap for *recA* marker gene. Black cells in the heatmap correspond to genomes sharing the same STs. The upper triangular region represents genomes sharing the same STs before the marker was excluded, and the lower triangular region represents genomes which are represented by the same STs after the exclusion of the *recA* gene from the scheme. (PDF 9.5MB)

**Supplementary File 13**. Sensitivity test heatmap for *gyrB* marker gene. Black cells in the heatmap correspond to genomes sharing the same STs. The upper triangular region represents genomes sharing the same STs before the marker was excluded, and the lower triangular region represents genomes which are represented by the same STs after the exclusion of the *gyrB* gene from the scheme. (PDF 9.5MB)

**Supplementary File 14**. Sensitivity test heatmap for *trpB* marker gene. Black cells in the heatmap correspond to genomes sharing the same STs. The upper triangular region represents genomes sharing the same STs before the marker was excluded, and the lower triangular region represents genomes which are represented by the same STs after the exclusion of the *trpB* gene from the scheme. (PDF 9.5MB)

**Supplementary File 15**. Sensitivity test heatmap for *rpoB* marker gene. Black cells in the heatmap correspond to genomes sharing the same STs. The upper triangular region represents genomes sharing the same STs before the marker was excluded, and the lower triangular region represents genomes which are represented by the same STs after the exclusion of the *rpoB* gene from the scheme. (PDF 9.5MB)

**Supplementary File 16**. Minimum Spanning Tree Analysis. This ZIP file contains all Python scripts and jupyter notebooks used to represent the extended Streptomyces scheme by calculating Minimum Spanning Tree, conducting initial investigations, and storing the generated outputs. Additionally, this ZIP file contains the XML file with Minimum Spanning Tree with fixed node positions, and LSPN taxonomy status. (ZIP 3MB)

**Supplementary File 17**. ANI analysis. This ZIP file contains all Python and bash scripts used to determine taxonomic boundaries among **Streptomyces** genomes sharing at least one identical MLST allele markers, and among Streptomyces genomes sharing the same names in NCBI. This ZIP file also contains all generated outputs and pyANI log files, and jupyter notebook used to resolve genus and species boundaries. (ZIP 31.1MB)

**Supplementary File 18**. Check pseudo genes. This ZIP file contains all Python scripts to check which marker sequences were reported from pseudo genes, and all generated outputs. (ZIP 4.8MB)

**Supplementary File 19**. MLSA analysis. This ZIP file contains all Python and bash scripts used for phylogenetic reconstruction of MLSA, and all generated outputs and log files. The ZIP file also contains a jupyter notebook used to investigate congruence between the Minimum Spanning Tree and MLSA phylogenetic tree. (ZIP 8.3MB)

**Supplementary File 20**: Generate and store figures using Python and R. This ZIP file containing data, additional Python and R scripts used to generate figures used in this manuscript. Additionally, the ZIP file contains PDF files of all figures used in the main text of this manuscript. (ZIP 10.1MB) 

**Supplementary File 21**. Relationship between the number of randomly sampled genomes 10-90% and number of disjoint graphs from the artificial scheme. (PDF 46KB)

**Supplementary File 22**. Distribution of relative connected component sizes from randomly samples 10-90% genomes from the artificial scheme. (PDF 26KB)

**Supplementary File 23**. Interactive Minimum Spanning Tree. HTML file containing interactive Minimum Spanning Tree of the extended pubMLST *Streptomyces* scheme with each node colour corresponding to the number of node/STs connections/degrees. (HTML 4.1MB)

**Supplementary File 24**. Interactive Minimum Spanning Tree. HTML file containing interactive Minimum Spanning Tree of the extended pubMLST *Streptomyces* scheme showing GenBank represented STs (green) and non-GenBank represented STs (blue). (HTML 4.1MB)

**Supplementary File 25**. Interactive Minimum Spanning Tree. HTML file containing interactive Minimum Spanning Tree of the extended pubMLST *Streptomyces* scheme showing existing pubMLST STs (green) and novel STs (blue). (HTML 4.1MB)

**Supplementary File 26**. Interactive Minimum Spanning Tree. HTML file containing interactive Minimum Spanning Tree of the extended pubMLST *Streptomyces* scheme showing number of unique genus per connected component. Each candidate genus is represented as a single node colour within a connected component. STs lacking a representative genome in NCBI are shown as grey nodes. (HTML 4.1MB)

**Supplementary File 27**. Interactive Minimum Spanning Tree. HTML file containing interactive Minimum Spanning Tree of the extended pubMLST *Streptomyces* scheme showing number of unique species per connected component. Each candidate species is represented as a single node colour within a connected component. STs lacking a representative genome in NCBI are shown as grey nodes. (HTML 4.1MB)

**Supplementary File 28**. Genomes sharing identical STs are assigned different names in NCBI. (PDF 7KB)

**Supplementary File 29**. ANIm genome coverage analysis of genomes found in the same group of connected STs. PDF file consists of a scatter plots showing genome coverage for all connected components comprising of at least two sequenced genomes. The number of unique species names found in a connected components is shown at the top of each plot. Within species comparisons (>=95% genome identity) are shown in blue, whereas between species comparisons   (>95% genome identity) are shown in red. The red horizontal line indicates the whole-genome genus threshold (50%). (PDF 1.2MB)

**Supplementary File 30**. ANIm genome identity analysis of genomes found in the same group of connected STs. Genome identity scatter plots for all connected components comprising of at least two sequenced genomes with the numbers at the top of each plots correspond to the number of unique species names found in a connected component. Within genus comparisons (>=50% genome identity) are shown in orange, whereas between genus comparisons (>50% genome identity) are shown in purple. The red horizontal line indicates the whole-genome species threshold (95%). (PDF 1.1MB)

**Supplementary File 31**.  Multiple STs are used to describe single Streptomyces species (>=50% genome coverage; >=95% genome identity). (PDF 9KB)

**Supplementary File 32**.  Multiple NCBI names are used to describe single Streptomyces species (>=50% genome coverage; >=95% genome identity). (PDF 8KB)

**Supplementary File 33**.  MLSA phylogenetic tree with original branch lengths. (PDF 96KB)




