#!/bin/bash
# Align marker sequences using mafft


# Markser sequences (Protein alignemnts)
#atpD was the only marker that had all 873 seuqnces (no pseudo genes) therefore will be saved to a diffrent folder
mafft ../output/protein/atpD_prot.fasta > ../output/alignments/all_sequences/protein/atpD_pro_aln.fasta
#rest of the markers will be saved to a diffrent folder
mafft ../output/protein/gyrB_prot.fasta > ../output/alignments/no_pseudo_genes/protein/gyrB_pro_aln.fasta

mafft ../output/protein/recA_prot.fasta > ../output/alignments/no_pseudo_genes/protein/recA_pro_aln.fasta

mafft ../output/protein/rpoB_prot.fasta > ../output/alignments/no_pseudo_genes/protein/rpoB_pro_aln.fasta

mafft ../output/protein/trpB_prot.fasta > ../output/alignments/no_pseudo_genes/protein/trpB_pro_aln.fasta

#Here, we also align 16S nulcoetide sequences

mafft ../output/nucleotide/16S_nucletide_found_STs.fasta > ../output/alignments/all_sequences/nucleotide/16S_nt_aln.fasta