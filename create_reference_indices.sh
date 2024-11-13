#!/bin/bash

mkdir -p /sp_library_refs
# Loop through all files in sp_libraries directory
for file in /sp_library_fasta/*; do
    # Check if it's a file (not a directory)
    if [ -f "$file" ]; then
        echo "Processing file: $file"
        # Extract just the filename without path and extension
        filename=$(basename "${file%.fasta}")
        bowtie2-build "$file" "/sp_library_refs/$filename"
    fi
done
