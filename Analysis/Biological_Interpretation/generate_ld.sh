#!/bin/bash
# Usage: bash run_ld.sh 15

CHR=$1

plink \
  --bfile /net/bmc-lab5/data/kellis/group/tanigawa/data/alkesgroup/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR} \
  --r2 \
  --ld-window 100000000 \
  --ld-window-kb 1000 \
  --ld-window-r2 0.1 \
  --out chr${CHR}_ld
