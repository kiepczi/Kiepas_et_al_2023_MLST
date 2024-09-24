#!/bin/bash
# blastn_and_tblastn.sh
# Srcipt used to run both blastn and tblastn
 

blastn -query ../output/blastn_tblastn_marker_genes/16S.fasta -subject ../output/genomes_for_BLAST/16S_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task blastn > ../output/BLAST_output/blastn_16S.txt

tblastn -query ../output/blastn_tblastn_marker_genes/atpD.fasta -subject ../output/genomes_for_BLAST/atpD_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ../output/BLAST_output/tblastn_atpD.txt
tblastn -query ../output/blastn_tblastn_marker_genes/gyrB.fasta -subject ../output/genomes_for_BLAST/gyrB_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ../output/BLAST_output/tblastn_gyrB.txt
tblastn -query ../output/blastn_tblastn_marker_genes/recA.fasta -subject ../output/genomes_for_BLAST/recA_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ../output/BLAST_output/tblastn_recA.txt
tblastn -query ../output/blastn_tblastn_marker_genes/rpoB.fasta -subject ../output/genomes_for_BLAST/rpoB_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ../output/BLAST_output/tblastn_rpoB.txt
tblastn -query ../output/blastn_tblastn_marker_genes/trpB.fasta -subject ../output/genomes_for_BLAST/trpB_new_blast.fasta -outfmt "6 std qcovs" -max_target_seqs 1000000000 -task tblastn > ../output/BLAST_output/tblastn_trpB.txt