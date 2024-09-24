#!/bin/bash
# Adding new sequences to existing alignment

mafft --add ../output/nucleotide_for_pseudo_genes/recA_nt_pseudo.fasta ../output/alignments/no_pseudo_genes/nucleotide/recA_nt_no_pseudo.fasta > ../output/alignments/all_sequences/nucleotide/recA_nt_aln.fasta
mafft --add ../output/nucleotide_for_pseudo_genes/gyrB_nt_pseudo.fasta ../output/alignments/no_pseudo_genes/nucleotide/gyrB_nt_no_pseudo.fasta > ../output/alignments/all_sequences/nucleotide/gyrB_nt_aln.fasta
mafft --add ../output/nucleotide_for_pseudo_genes/rpoB_nt_pseudo.fasta ../output/alignments/no_pseudo_genes/nucleotide/rpoB_nt_no_pseudo.fasta > ../output/alignments/all_sequences/nucleotide/rpoB_nt_aln.fasta
mafft --add ../output/nucleotide_for_pseudo_genes/trpB_nt_pseudo.fasta ../output/alignments/no_pseudo_genes/nucleotide/trpB_nt_no_pseudo.fasta > ../output/alignments/all_sequences/nucleotide/trpB_nt_aln.fasta