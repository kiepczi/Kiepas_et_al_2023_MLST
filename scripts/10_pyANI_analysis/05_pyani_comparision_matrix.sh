#!/bin/bash
#Get pyANI matrices comparisions


for i in {1..117}
do
   pyani report -o ~/Desktop/Kiepas_et_al_2023_MLST/data/pyani_analysis/connected_components_matrices --formats=stdout --run_matrices $i
done