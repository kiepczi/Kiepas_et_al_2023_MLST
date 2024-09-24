#!/bin/bash
# download_genomes
# Download genomes from NCBI

# $1 - genera for which the genomes will be downloaded
# $2 - output file (path)


ncbi-genome-download bacteria -F all -l all -g $1 -o $2
