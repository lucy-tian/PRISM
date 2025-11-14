#!/bin/bash
#SBATCH --job-name=finemap_array
#SBATCH --output=logs/finemap_%A_%a.out
#SBATCH --error=logs/finemap_%A_%a.err
#SBATCH -p kellis
#SBATCH --time=2-00:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-55

set -euo pipefail

# Optional: activate your conda env
# source activate polyfun

# Define traits
TRAITS=("INI1003063" "INI20030780" "INI30120" "INI50030700")

# Collect all part files in trait order
CHUNK_FILES=()
for TRAIT in "${TRAITS[@]}"; do
    for FILE in output/${TRAIT}_job_chunks/part_*; do
        CHUNK_FILES+=("$FILE")
    done
done

# Determine chunk file to run
CHUNK_FILE="${CHUNK_FILES[$SLURM_ARRAY_TASK_ID]}"
echo "Running chunk file: $CHUNK_FILE"

# Run each line (command) in the chunk file
while IFS= read -r CMD || [ -n "$CMD" ]; do
    # Extract LD block info from CMD
    CHR=$(echo "$CMD" | grep -oP '(?<=--chr )\d+')
    START=$(echo "$CMD" | grep -oP '(?<=--start )\d+')
    END=$(echo "$CMD" | grep -oP '(?<=--end )\d+')
    OUT=$(echo "$CMD" | grep -oP '(?<=--out )\S+')
    TRAIT=$(echo "$OUT" | awk -F'/' '{print $(NF-1)}')

    echo ">>> Running trait: $TRAIT | LD block: chr${CHR}:${START}-${END}"

    eval "$CMD"
done < "$CHUNK_FILE"
