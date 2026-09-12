#!/bin/bash
set -beEuo pipefail

# this script applies plink2 --sscore and compute
# PGS for each individual
echo "executed"
SRCNAME=$(readlink -f "${0}")
echo "$SRCNAME"
SRCDIR=$(dirname "${SRCNAME}")

source $(dirname ${SRCDIR})/helpers/functions.sh
source_paths "${SRCDIR}"


############################################################
# tmp dir
############################################################
tmp_dir="$(mktemp -p $(dirname $(mktemp -u)) -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

############################################################

phenotype=$1
family=$2
refit_str=$3       # "fit_w_val" "refit"
analysis_name=$4
covariates_type=$5 # "pop_10PCs" "UKB_18PCs" "hgdpKGP_10PCs"
p_factor_file=$6   # "None" or path to .rds file
sample_wts_file=$7 # "None" or path to .rds file
pop=$8             # "all", "w_others", "4pops", "WB"
save_d=${9}


############################################################
if [ ! -d "${save_d}" ] ; then mkdir -p "${save_d}" ; fi


echo ${phenotype_file}


refit_rdata_file="None"
echo ${keep_file}
echo ${sscore_afreq_ref}
echo ${save_d}
echo ${phenotype_file}
echo ${keep_file}
echo ${phenotype}
echo ${family}
echo ${genotype_pfile}
echo ${covariates_type}
echo ${p_factor_file}
echo ${sample_wts_file}
echo ${keep_file}


if [ ! -f "${save_d}/snpnet.BETAs.tsv.gz" ] ; then
    if [ -f "${save_d}/snpnet.RData" ] &&
       [ -f "${save_d}/snpnet.metrics.tsv.gz" ] ; then
        Rscript ${SRCDIR}/export_betas.R ${save_d}
    else
        rsync -az "${snpnet_repo_d}" "${save_d}/"

echo ${SRCDIR}

        OMP_NUM_THREADS=12 Rscript ${SRCDIR}/fit_snpnet_v3.R \
            ${save_d} \
            ${phenotype_file} \
            ${keep_file} \
            ${phenotype} \
            ${family} \
            ${refit_rdata_file} \
            ${genotype_pfile} \
            ${covariates_type} \
            ${p_factor_file} \
            ${sample_wts_file} 2>&1 \
        | tee -a ${save_d}/snpnet.fit.log
    fi
fi


if [ ! -f "${save_d}/snpnet.sscore.zst" ] ; then
    bash "${SRCDIR}/plink2_sscore.sh" "${save_d}/snpnet.BETAs.tsv.gz" "${sscore_afreq_ref}"
fi


if [ ! -d "${data_d}/406k_covars_${covariates_type}/${phenotype}" ] ; then
   Rscript ${SRCDIR}/eval_covars.R ${phenotype} ${family} ${covariates_type} all ${phenotype_file} ${pheno_names}
fi

if [ ! -f "${data_d}/406k_covars_${covariates_type}/covars_lead.score.tsv.gz" ] ; then
   Rscript ${SRCDIR}/eval_covars.R ${phenotype} ${family} ${covariates_type} all ${phenotype_file} ${pheno_names}
fi

if [ ! -f "${save_d}/snpnet.eval.tsv.gz" ] ; then
Rscript ${SRCDIR}/eval_PGS.R ${save_d} ${phenotype} ${family} ${covariates_type} snpnet ${phenotype_file} ${pheno_names}
fi
