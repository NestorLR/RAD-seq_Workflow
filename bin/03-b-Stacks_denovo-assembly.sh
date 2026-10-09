#!/bin/bash

# Exit on error, undefined variables, or pipeline failures
set -euo pipefail

# ============================================
# De novo assembly of RAD-seq data (Stacks)
# ============================================
# This step runs the de novo pipeline to assemble loci and call SNPs without a reference genome.

#Inputs and Outputs directories
STACKS_DIR="data/stacks/process_radtags_fixed"
POPMAP="metadata/popMap.txt"
OUTPUT_DIR="data/stacks/denovo_M2_n2"

# Number of CPU threads
THREADS=2

# Create output directory
mkdir -p "$OUTPUT_DIR"


denovo_map.pl \
  --samples "$STACKS_DIR" \
  --popmap "$POPMAP" \
  -T "$THREADS" \
  -m 3 \
  -M 2 \
  -n 2 \
  --paired \
  --out-path "$OUTPUT_DIR"

# ============================================
# Notes:
# ============================================
# - This configuration represents a commonly used baseline parameter set for RAD-seq data.
# - The parameters (-M, -n) influence locus clustering stringency and can affect SNP recovery.
# -M número máximo de diferencias/mismatches permitidos entre stacks dentro de un individuo para agruparlos en un locus
# -n Controla cuántas diferencias se permiten al correlacionar los loci entre individuos para construir el catálogo
# -m 3' cobertura mínima de reads necesaria para formar un stack en ustacks
# - The option '--force-diff-len' allows inclusion of loci with indel variation,
#   which may be important for non-model organisms.
# - Parameter optimization (e.g., testing multiple M/n combinations) is recommended
#   to assess robustness of downstream population genomic inferences.
#  -X "ustacks: --force-diff-len" \  # allow loci of different lengths (useful for indel variation)
# --popmap data/stacks/popMap.txt archivo que indica a qué población corresponde la muestra 
