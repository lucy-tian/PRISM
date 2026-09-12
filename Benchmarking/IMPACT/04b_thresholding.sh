#!/bin/bash
set -euo pipefail

ukb21942_d='/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942'
trait=${1}   # loop over your traits as needed
threshold=${2}

plink2 \
    --bfile genotype/Afr \
    --keep ${ukb21942_d}/sqc/population.20220316/Afr_train_val.keep \
    --score scoring/${1}_${2}.tsv 1 2 3 header ignore-dup-ids \
    --out prs/${1}_${2}