while read line; do
    sed "s/SAMPLE_NAME/${line}/g" inputs.template.json > inputs.json
    ~/bin/java -jar ~/bin/cromwell-29.jar run \
        -i inputs.json \
        epitope_prediction.wdl
done < $1;