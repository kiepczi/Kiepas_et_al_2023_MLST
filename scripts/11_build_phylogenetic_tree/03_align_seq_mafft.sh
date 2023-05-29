#!/bin/bash
# align_sco_mafft.sh
# Align marker sequences using mafft


# Align each set of SCOs
for fname in ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/protein/*; do

    mafft $fname > ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/alignments/protein/` basename ${fname%%_prot.fasta}`_aln.fasta
done

mafft ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/nucleotide/16S_nucletide_found_STs.fasta > ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/alignments/nucleotide/16S_aln.fasta