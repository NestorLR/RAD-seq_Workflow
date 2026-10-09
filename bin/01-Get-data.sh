#!/bin/bash

#This script downloads raw data from SRA database

ACCESSIONS_FILE="metadata/PRJEB91610_AccesionList.txt"
OUTPUT_DIR="data/raw/"

mkdir -p $OUTPUT_DIR

mapfile -t SRA_IDS < "$ACCESSIONS_FILE"

for SRA_ID in "${SRA_IDS[@]}"
do
    echo "Downloading $SRA_ID"

    # download .fastq file
    fasterq-dump $SRA_ID -O $OUTPUT_DIR

    # compress
    gzip $OUTPUT_DIR/*.fastq

    echo "Finished $SRA_ID"
done

