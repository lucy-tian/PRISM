#!/bin/bash

RSCRIPT="/net/bmc-lab5/data/kellis/users/lucytian/0_SOFTWARE/miniconda3/envs/sbayesrc/bin/Rscript"

mkdir -p logs
for trait in INI30120 INI20030780 INI1003063 INI50030700; do
  sbatch -p kellis \
    -J "mix_${trait}" \
    -o "logs/mix_${trait}.%j.out" \
    -e "logs/mix_${trait}.%j.err" \
    --wrap="Rscript 04b_combine.R ${trait}"
done