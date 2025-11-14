# List of traits
TRAITS=("INI1003063" "INI20030780" "INI30120" "INI50030700")

# Paths
SUMSTAT_DIR="sumstat"
OUTDIR="output"

for TRAIT in "${TRAITS[@]}"; do
    echo "Aggregating results for $TRAIT..."


    python /home/lucytian/data/0_SOFTWARE/polyfun/aggregate_finemapper_results.py \
        --out-prefix ${OUTDIR}/${TRAIT}/finemap \
        --sumstats ${SUMSTAT_DIR}/${TRAIT}_sumstats.txt.gz \
        --out ${OUTDIR}/${TRAIT}/output.agg.txt.gz \
        --allow-missing-jobs \
        --adjust-beta-freq
        
done