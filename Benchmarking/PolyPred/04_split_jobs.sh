#!/bin/bash

TRAITS=("INI1003063" "INI20030780" "INI30120" "INI50030700")

for TRAIT in "${TRAITS[@]}"; do
    JOB_FILE="output/${TRAIT}_jobs.txt"
    CHUNK_DIR="output/${TRAIT}_job_chunks"
    mkdir -p "$CHUNK_DIR"

    # Remove any previous chunks just in case
    rm -f "$CHUNK_DIR"/part_*

    # Split into 500-line chunks, files named part_00, part_01, etc.
    split -l 200 "$JOB_FILE" "$CHUNK_DIR/part_"
    
    echo "Split $JOB_FILE into $(ls $CHUNK_DIR | wc -l) chunks at $CHUNK_DIR"
done
