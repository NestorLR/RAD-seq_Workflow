#!/bin/bash

# Exit on error, undefined variables, or pipeline failures
set -euo pipefail

# ============================================
# Quality filtering of demultiplexed RAD-seq reads
# ============================================
# This is the first  step in Stacks. 
# Removes reads containing uncalled bases and permforms quality filtering.


# Input and output directories
RAW_DIR="data/raw"
STACKS_DIR="data/stacks/process_radtags"

mkdir -p "$STACKS_DIR"


echo "============================================"
echo "Running process_radtags"
echo "Restriction enzymes: EcoRI-HF + MspI"
echo "Paired-end data"
echo "============================================"

for R1 in "$RAW_DIR"/*_1.fastq.gz; do

    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="$RAW_DIR/${SAMPLE}_2.fastq.gz"

    echo "--------------------------------------------"
    echo "Processing sample: $SAMPLE"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "--------------------------------------------"

    if [[ ! -f "$R2" ]]; then
        echo "ERROR: R2 file not found for $SAMPLE"
        exit 1
    fi

    process_radtags \
      -1 "$R1" \
      -2 "$R2" \
      -o "$STACKS_DIR" \
      --renz-1 ecoRI \
      --renz-2 mspI \
      -c \
      -q \
      --basename "$SAMPLE" \
      -r

done

echo "============================================"
echo "process_radtags completed successfully."
echo "Output: $STACKS_DIR"
echo "============================================"



# ============================================
# Notes:
# ============================================
# - The input directory (-p) should contain demultiplexed FASTQ files.
# - La opción '-r' se refiere al rescate de barcodes y/o sitios de restricción cuando se detectan errores que pueden ser corregidos bajo los criterios del programa.
# - This step is critical for reducing sequencing errors prior to SNP calling.
# -q perform quality filtering using a sliding window
# -c remove reads with uncalled bases
