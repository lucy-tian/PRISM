#!/bin/bash

IMPACT_d='/home/lucytian/data/6_IMPACT/top_SNPs'
trait=$1

plink \
  --bfile genotype/Afr \
  --clump ${IMPACT_d}/${trait}.topsnps.tsv \
  --clump-p1 0.5 \
  --clump-p2 1 \
  --clump-r2 0.2 \
  --clump-snp-field 'ID' \
  --clump-field 'P' \
  --out clumping/${1}

