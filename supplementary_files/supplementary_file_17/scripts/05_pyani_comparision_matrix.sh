#!/bin/bash
#Get pyANI matrices comparisions


for i in {1..116}
do
   pyani report -o ../output/pyani_matrices_connected_components --formats=stdout --run_matrices $i
done

for i in {117..196}
do
   pyani report -o ../output/pyani_matrices_shared_names --formats=stdout --run_matrices $i
done