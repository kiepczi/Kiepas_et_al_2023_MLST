for FILE in ../input/connected_components/*; do

    file="${FILE##}"
    job_id=`echo "$FILE" | cut -d'/' -f4`
    echo $job_id 
    pyani anim -i $file -o ${file}/output -v -l ${file}/output.log --name $job_id --labels ${file}/custom_labels.txt --classes ${file}/custom_classes.txt
    done

for FILE in ../input/shared_NCBI_names/*; do

    file="${FILE##}"
    job_id=`echo "$FILE" | cut -d'/' -f4`
    pyani anim -i $file -o ${file}/output -v -l ${file}/output.log --name $job_id --labels ${file}/custom_labels.txt --classes ${file}/custom_classes.txt
    done