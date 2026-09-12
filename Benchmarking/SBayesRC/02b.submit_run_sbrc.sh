#!/bin/bash

RSCRIPT="/net/bmc-lab5/data/kellis/users/lucytian/0_SOFTWARE/miniconda3/envs/sbayesrc/bin/Rscript"


mkdir -p logs
for pop in WB Afr; do
  for trait in INI30120 INI20030780 INI1003063 INI50030700; do
    sbatch -p kellis --time=1-00:00:00 --mem=96G --cpus-per-task=4 \
      -J "sbrc_${pop}_${trait}" \
      -o "logs/sbrc_${pop}_${trait}.%j.out" \
      -e "logs/sbrc_${pop}_${trait}.%j.err" \
      --wrap="Rscript 02a_run_sbrc.R ${pop} ${trait}"
  done
done

