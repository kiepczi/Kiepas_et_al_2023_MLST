#!/bin/bash
# align_subsequences_mafft.sh
# Srcipt used to align subsequences and extracted whole sequences

#atpD sequences
cat ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/*atpD*fasta > ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/atpD_all.fasta
mafft --auto ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/atpD_all.fasta > ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/atpD_aligned.fasta
rm -r ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/atpD_all.fasta

#trpB sequences
cat ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/*trpB*fasta > ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/trpB_all.fasta
mafft --auto ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/trpB_all.fasta > ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/trpB_aligned.fasta
rm -r ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/trpB_all.fasta