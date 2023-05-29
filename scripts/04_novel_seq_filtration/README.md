# README.md - novel_alleles_pubmlst

## **Overview**
This repository contains a script used to extend an already existing pubMLST allele sequences FASTA file with [MLST tool](https://github.com/tseemann/mlst) acceptable format. 




## **Arguments**
To run the following script pass the following arguments:
novel_sequences = sys.argv[1]
current_marker_alleles = sys.argv[2]
save_new_marker_alleles = sys.argv[3]



- `pubMLST_streptomyces_alleles`: FASTA files (`.fas`) with existing pubMLST allele sequences for genus *Streptomyces*. These files have been downloaded from [pubMLST](https://pubmlst.org/bigsdb?db=pubmlst_streptomyces_seqdef&page=downloadAlleles).

- `arg1`- FASTA file with novel alleles identified by [MLST tool](https://github.com/tseemann/mlst). 

- `arg2` - path to a folder with current set of marker allele sequences

- `arg3` - path to a file, where extended version of FASTA with all novel and current marker allele sequences will be saved to. 



## **Example:**

There are 154 16S rRNA sequences for genus *Streptomyces* in [pubMLST](https://pubmlst.org/bigsdb?db=pubmlst_streptomyces_seqdef&page=downloadAlleles) database. There were 31 novel 16S rRNA sequences identified by [MLST tool](https://github.com/tseemann/mlst). After running `novel_alleles.py`, it will return a FASTA file with already existing allele sequences for 16S rRNA locus (named as `>16S_1` to `>16S_154`) and novel allele sequences (names as `>16S_155` to `>16S_184`). 


<br>
<br>

**NOTE**:

The [MLST tool](https://github.com/tseemann/mlst) allows to insert a new scheme into a local `mlst` database. These files must have `.tfa` extension and the FASTA sequence IDs must be names as `>allele_number`. 

Novel alleles that contain ambiguity symbols are discarded!


# Example usage

```bash
python 01_novel_seq_filtration.py ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_tool_output/strep_mlst1_novel_genes.fasta ~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/pubMLST_streptomyces_scheme/ ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/schemes/strep_mlst_round_1 ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences
```