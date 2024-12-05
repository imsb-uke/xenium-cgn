mkdir -p outputs

for csv_file in transcripts_codes/*.csv; do
    filename=$(basename "$csv_file" .csv)
    
    echo "Processing $filename"
    
    docker run \
        -v "$(pwd)/transcripts_codes:/data/input" \
        -v "$(pwd)/outputs:/data/output" \
        vpetukhov/baysor:latest \
        baysor run \
        -x x_location \
        -y y_location \
        -z z_location \
        -g feature_name \
        -m 10 \
        -p \
        --prior-segmentation-confidence 0.5 \
        --save-polygons=geojson \
        "/data/input/$(basename $csv_file)" \
        :cell_id_codes \
        -o "/data/output/output-$filename"
    
    if [ $? -eq 0 ]; then
        echo "Successfully processed $filename"
    else
        echo "Error processing $filename"
    fi
done
