
#!/bin/bash

# Exit on error, undefined variables, or pipeline failures
set -euo pipefail

# ============================================
# SNP filtering and VCF export with Stacks
# ============================================

# Input paths
STACKS_DIR="data/stacks/denovo_M2_n2"
POPMAP="metadata/popMap.txt"

# Filtering threshold:
# Minimum proportion of individuals across populations
R=0.66

# Output directory
OUTDIR="data/stacks/populations_R70"

# Create output directory
mkdir -p "$OUTDIR"


echo "============================================"
echo "Running Stacks populations"
echo "Minimum overall sample proportion: $R"
echo "Output directory: $OUTDIR"
echo "============================================"

populations \
    -P "$STACKS_DIR" \
    -M "$POPMAP" \
    -R "$R" \
    --vcf \
    --fstats \
    --min-mac 3 \
    --min-maf 0.05 \
    --min-gt-depth 3 \
    --write-single-snp \
    -O "$OUTDIR"
    
echo "============================================"
echo "Analysis completed."
echo "Results: $OUTDIR"
echo "============================================"

