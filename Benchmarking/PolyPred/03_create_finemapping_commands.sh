#!/bin/bash
set -euo pipefail

# List of traits
TRAITS=("INI1003063" "INI20030780" "INI30120" "INI50030700")

# Paths
SUMSTAT_DIR="sumstat"
OUTDIR="output"
REGIONS_FILE="regions_file.tsv"
N=270920  # sample size

for TRAIT in "${TRAITS[@]}"; do
    echo "Creating jobs for $TRAIT..."

    python /home/lucytian/data/0_SOFTWARE/polyfun/create_finemapper_jobs.py \
        --sumstats ${SUMSTAT_DIR}/${TRAIT}_sumstats.txt.gz \
        --n ${N} \
        --method susie \
        --max-num-causal 5 \
        --out-prefix ${OUTDIR}/${TRAIT}/finemap \
        --jobs-file ${OUTDIR}/${TRAIT}_jobs.txt \
        --allow-missing \
        --regions-file ${REGIONS_FILE}
done
