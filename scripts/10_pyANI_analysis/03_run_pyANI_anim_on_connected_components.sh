# for FILE in ~/Desktop/Kiepas_et_al_2023_MLST/data/pyani_analysis/connected_components/*; do

#     file="${FILE##}"
#     job_id=`echo "$FILE" | cut -d'/' -f10`
#     echo $job_id
    # pyani anim -i $file -o ${file}/output -v -l ${file}/output.log --name $job_id --labels ${file}/labels.txt --classes ${file}/custom_labels.txt
#     done

for FILE in ~/Desktop/Kiepas_et_al_2023_MLST/data/pyani_analysis/shared_NCBI_names/*; do

    file="${FILE##}"
    job_id=`echo "$FILE" | cut -d'/' -f9`
    echo $job_id
    pyani anim -i $file -o ${file}/output -v -l ${file}/output.log --name $job_id --labels ${file}/custom_labels.txt --classes ${file}/custom_classes.txt
    done