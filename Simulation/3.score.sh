#!/bin/bash
set -beEuo pipefail

# 1. Input: Expects a directory path as the first argument
input_dir=$1
beta_f="${input_dir%/}/snpnet.BETAs.tsv.gz" # ${input_dir%/} removes trailing slash if present

ml load plink2

# 2. Define output prefix and trait name (the folder name)
out="${beta_f%.BETAs.tsv.gz}"
colname=$(basename "$(dirname "$beta_f")")

# 3. Execute PLINK2
# We use <() to feed the SNP list and the cleaned BETA file directly into memory
plink2 \
  --threads 6 \
  --memory 40000 \
  --pfile synthetic_v1_Afr_chr22only \
  --read-freq synthetic_freq.afreq \
  --extract <(zcat "$beta_f" | awk 'NR>1 {print $1}') \
  --score <(zcat "$beta_f" | awk -v FS='\t' 'NR==1 || $2 != ""' | sed "s/BETA/${colname}/g") \
  1 2 3 header-read zs \
  cols=maybefid,maybesid,phenos,nallele,dosagesum,scoreavgs,denom,scoresums \
  --out "$out"

# 4. Standardize the log name
if [ -f "${out}.log" ]; then
    mv "${out}.log" "${out}.sscore.log"
fi