#!/bin/bash
# download_genomes
# Download genomes from NCBI

# $1 - genera for which the genomes will be downloaded
# $2 - output file (path)


# $1  path to output directory
# MLST
ncbi-genome-download bacteria -F all -l all --genera $1 -o $2
