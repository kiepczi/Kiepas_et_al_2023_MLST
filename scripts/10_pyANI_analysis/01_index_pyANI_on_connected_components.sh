#!/bin/bash
#Run pyANI v3


# for FILE in ~/Desktop/Kiepas_et_al_2023_MLST/data/pyani_analysis/connected_components/*; do
#     pyani index -i $FILE
#     done


for FILE in ~/Desktop/Kiepas_et_al_2023_MLST/data/pyani_analysis/shared_NCBI_names/*; do
    pyani index -i $FILE
    done