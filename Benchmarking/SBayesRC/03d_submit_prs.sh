#!/bin/bash

RSCRIPT="/net/bmc-lab5/data/kellis/users/lucytian/0_SOFTWARE/miniconda3/envs/sbayesrc/bin/Rscript"

mkdir -p logs
for pop in WB Afr; do
  for trait in INI30120 INI20030780 INI1003063 INI50030700; do
    sbatch -p kellis \
      -J "prs_${pop}_${trait}" \
      -o "logs/prs_${pop}_${trait}.%j.out" \
      -e "logs/prs_${pop}_${trait}.%j.err" \
      --wrap="Rscript 03c_PRS_pred.R ${pop} ${trait}"
  done
done
