#!/bin/bash
# trim_gaps.sh
# This scrip was used to trim the alignmets gaps with trimal with -colnumbering parameter and to generate file with the matchin columns in the new and old alignment


# Triming the alignment
# arg1 - input alignment
# arg2 - output alignment
# arg3 - output colmatch

trimal -in ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/alignments/concatenated/concatenated_alignment.fasta -out ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/alignments/concatenated/trimmed_concantenated_alignment.fasta -automated1 -colnumbering > ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/alignments/concatenated/trimmed_columns_info.txt