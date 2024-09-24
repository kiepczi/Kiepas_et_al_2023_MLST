#!/bin/bash
# mlst_v2.sh
# Running MLST v2 tool

# $1 - novel sequences file (path to file)
# $2 - input gbfiles (path to file)
# $3 - MLST scheme
# $4 - MLST output file


# MLST
mlst --novel $1 $2 --legacy --scheme $3 --mincov 80 >  $4


