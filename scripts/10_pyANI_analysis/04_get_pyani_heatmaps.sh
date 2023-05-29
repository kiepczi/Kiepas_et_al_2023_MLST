#!/bin/bash
#Get pyANI matrices comparisions



# for i in {1..117}
# do
#    pyani plot -o ~/Desktop/Kiepas_et_al_2023_MLST/figures/pyani_analysis_connected_components --run_id $i -v --formats pdf
# done


for i in {177..259}
do
   pyani plot -o ~/Desktop/Kiepas_et_al_2023_MLST/figures/pyani_analysis_NCBI_shared_names --run_id $i -v --formats pdf
done