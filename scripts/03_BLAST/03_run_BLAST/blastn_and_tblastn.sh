#!/bin/bash
# blastn_and_tblastn.sh
# Srcipt used to run both blastn and tblastn
 

blastn -query ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blastn_tblastn_marker_genes/16S.fasta -subject ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/databases_for_BLAST/16S_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task blastn > ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/genome_filtration/blastn_16S.txt

# tblastn -query ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blastn_tblastn_marker_genes/atpD.fasta -subject ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/databases_for_BLAST/atpD_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/tblastn_atpD.txt

tblastn -query ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blastn_tblastn_marker_genes/gyrB.fasta -subject ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/databases_for_BLAST/gyrB_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/genome_filtration/tblastn_gyrB.txt
tblastn -query ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blastn_tblastn_marker_genes/recA.fasta -subject ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/databases_for_BLAST/recA_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/genome_filtration/tblastn_recA.txt
tblastn -query ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blastn_tblastn_marker_genes/rpoB.fasta -subject ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/databases_for_BLAST/rpoB_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/genome_filtration/tblastn_rpoB.txt
tblastn -query ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blastn_tblastn_marker_genes/trpB.fasta -subject ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/databases_for_BLAST/trpB_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ~/Desktop/Kiepas_et_al_2023_MLST/data/BLAST_analysis/blast_output/genome_filtration/tblastn_trpB.txt