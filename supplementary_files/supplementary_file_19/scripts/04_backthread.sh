#!/bin/bash
# backtranslate_alignments.sh
# Backthread nucleotide sequences onto protein alignments

# Note: atpD was the only one marker with no pseudo sequences so will be saved into a diffrent directory
t_coffee -other_pg seq_reformat -in ../output/nucleotide/atpD_nt.fasta -in2 ../output/alignments/all_sequences/protein/atpD_pro_aln.fasta -output fasta -action +thread_dna_on_prot_aln -output fasta > ../output/alignments/all_sequences/nucleotide/atpD_nt_aln.fasta 
#Rest of the markers will be save to a diffrent folder so that we can add rest of the sequences
t_coffee -other_pg seq_reformat -in ../output/nucleotide/gyrB_nt.fasta -in2 ../output/alignments/no_pseudo_genes/protein/gyrB_pro_aln.fasta -output fasta -action +thread_dna_on_prot_aln -output fasta > ../output/alignments/no_pseudo_genes/nucleotide/gyrB_nt_no_pseudo.fasta
t_coffee -other_pg seq_reformat -in ../output/nucleotide/recA_nt.fasta -in2 ../output/alignments/no_pseudo_genes/protein/recA_pro_aln.fasta -output fasta -action +thread_dna_on_prot_aln -output fasta > ../output/alignments/no_pseudo_genes/nucleotide/recA_nt_no_pseudo.fasta
t_coffee -other_pg seq_reformat -in ../output/nucleotide/rpoB_nt.fasta -in2 ../output/alignments/no_pseudo_genes/protein/rpoB_pro_aln.fasta -output fasta -action +thread_dna_on_prot_aln -output fasta > ../output/alignments/no_pseudo_genes/nucleotide/rpoB_nt_no_pseudo.fasta
t_coffee -other_pg seq_reformat -in ../output/nucleotide/trpB_nt.fasta -in2 ../output/alignments/no_pseudo_genes/protein/trpB_pro_aln.fasta -output fasta -action +thread_dna_on_prot_aln -output fasta > ../output/alignments/no_pseudo_genes/nucleotide/trpB_nt_no_pseudo.fasta