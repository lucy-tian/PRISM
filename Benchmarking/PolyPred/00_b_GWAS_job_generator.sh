#!/bin/bash
set -euo pipefail

OUTDIR="./glm_job_lists"
mkdir -p "${OUTDIR}"

plink_script="/home/lucytian/data/4_polypred/00_a_GWAS_glm.sh" 
threads=2
memory=12000
suffix="QTs"
covar_str="pop_10PCs"

declare -A configs
configs[blood_biochemistry]="21 34"
configs[blood_count]="15"
configs[spirometry]="11"

keep_sets=("WB_train_val" "Afr_train_val")

for analysis_name in "${!configs[@]}"; do
    job_list="${OUTDIR}/${analysis_name}.jobs.sh"
    > "${job_list}"  # clear file if exists

    for col in ${configs[$analysis_name]}; do
        for keep in "${keep_sets[@]}"; do
            for idx in $(seq 1 69); do
                echo "bash ${plink_script} ${analysis_name} ${suffix} ${col} ${threads} ${memory} ${covar_str} ${keep} ${idx}" >> "${job_list}"
            done
        done
    done

    echo "✅ Generated: ${job_list}"
done
