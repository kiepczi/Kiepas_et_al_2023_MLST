#!/bin/bash
# backtranslate_alignments.sh
# Backthread nucleotide sequences onto protein alignments

# Create output directory

# Align each set of SCOs
for fname in ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/marker_sequences/nucleotide/*; do
    if [[ $fname != *'16S'* ]]; then
        t_coffee -other_pg seq_reformat -in ${fname} -in2 ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/alignments/protein/` basename ${fname%%_nt.fasta}`_aln.fasta -output fasta -action +thread_dna_on_prot_aln -output fasta > ~/Desktop/Kiepas_et_al_2023_MLST/data/phylogeny/alignments/nucleotide/` basename ${fname%%_nt.fasta}`_nt_aln.fasta

    fi

done
