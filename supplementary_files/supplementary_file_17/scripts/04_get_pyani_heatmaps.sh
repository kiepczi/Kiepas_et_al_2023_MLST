#!/bin/bash
#Get pyANI matrices comparisions



for i in {1..116}
do
   pyani plot -o ../output/pyani_heatmaps_connected_components --run_id $i -v --formats pdf
done


for i in {117..196}
do
   pyani plot -o ../output/pyani_heatmaps_shared_names --run_id $i -v --formats pdf
done