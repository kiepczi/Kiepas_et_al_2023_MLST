#!/bin/bash
#Run pyANI v3


for FILE in ../input/connected_components/*; do
    pyani index -i $FILE
    done


for FILE in ../input/shared_NCBI_names/*; do
    pyani index -i $FILE
    done