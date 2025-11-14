#!/bin/bash
set -euo pipefail

ukb21942_d='/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942'


plink2 \
    --bfile genotype/Afr \
    --keep ${ukb21942_d}/sqc/population.20220316/Afr_test.keep \
    --score scoring/INI1003063_0_0001.tsv 1 2 3 header ignore-dup-ids \
    --out prs/best_INI1003063
    
plink2 \
    --bfile genotype/Afr \
    --keep ${ukb21942_d}/sqc/population.20220316/Afr_test.keep \
    --score scoring/INI20030780_3e-05.tsv 1 2 3 header ignore-dup-ids \
    --out prs/best_INI20030780

plink2 \
    --bfile genotype/Afr_exclude_APOE \
    --keep ${ukb21942_d}/sqc/population.20220316/Afr_test.keep \
    --score scoring/INI20030780_3e-05.tsv 1 2 3 header ignore-dup-ids \
    --out prs/best_INI20030780_exclude_APOE


plink2 \
    --bfile genotype/Afr \
    --keep ${ukb21942_d}/sqc/population.20220316/Afr_test.keep \
    --score scoring/INI30120_0_1.tsv 1 2 3 header ignore-dup-ids \
    --out prs/best_INI30120