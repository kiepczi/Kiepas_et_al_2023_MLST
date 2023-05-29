#!/bin/bash
# blastn_and_tblastn.sh
# Srcipt used to run both blastn and tblastn


#Arguments for tblastn 
blastn -query ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/atpD_subsequences.fasta -subject ~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/streptomyces_genomes/GCF_009862685.1/GCF_009862685.1_ASM986268v1_genomic.fna -outfmt "6 std qcovs" -max_target_seqs 1000000000 > ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/subsequences/blastn_atpD_subsequence.txt
blastn -query ~/Desktop/Kiepas_et_al_2023_MLST/data/MLST_scheme_extension/allele_marker_subsequences/trpB_subsequences.fasta -subject ~/Desktop/Kiepas_et_al_2023_MLST/data/raw_data/NCBI_genomes/streptomyces_genomes/GCF_016860545.1/GCF_016860545.1_ASM1686054v1_genomic.fna -outfmt "6 std qcovs" -max_target_seqs 1000000000 > ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/subsequences/blastn_trpB_subsequence.txt