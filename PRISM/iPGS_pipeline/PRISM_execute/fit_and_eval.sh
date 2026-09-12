#!/bin/bash
#SBATCH -o %j.out 
#SBATCH -e %j.err

set -beEuo pipefail

SRCNAME=$(readlink -f "${0}")
SRCDIR=$(dirname "${SRCNAME}")

source $(dirname ${SRCDIR})/helpers/functions.sh
source_paths "${SRCDIR}"

phenotype=$1
refit_str="fit_w_val"

p_factor_file="${p_factors}/${tissue}_cV2F.pfactor.rds"

sample_wts_file="None"
pop="all"                   # all WB
#genotype="genoHM3"          # geno genoHM3 cV2Fonly codingonly
#covariates_type="UKB_18PCs" # pop_10PCs UKB_18PCs hgdpKGP_10PCs
family="gaussian"
analysis_name=$(basename ${SRCDIR})

#save_d="${data_d}/${cohort}/${analysis_name}/${refit_str}/param_${param}/${phenotype}"

save_d="${data_d}/${analysis_name}/${refit_str}/${phenotype}"



echo ${save_d}
echo ${family}
echo ${analysis_name}


bash $(dirname ${SRCDIR})/helpers/fit_snpnet_eval_v3.sh \
    ${phenotype} \
    ${family} \
    ${refit_str} \
    ${analysis_name} \
    ${covariates_type} \
    ${p_factor_file} \
    ${sample_wts_file} \
    ${pop} \
    ${save_d} \

exit 0