#!/bin/bash
set -beEuo pipefail

if [ $# -lt 8 ] ; then
    echo "missing args!" >&2
    exit 1
fi


ukb21942_d='/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942'


analysis_name=$1
out_suffix=$2     # QTs
pheno_col_nums=$3 # 3-253
plink_threads=$4  # 2
plink_memory=$5   # 12000
covar_str=$6 # pop_10PCs pop_10PCs_no_phequantnorm UKB_18PCs
keep_set=$7  # WB
idx=$8       # 9
use_tmp_geno="FALSE"
overwrite="FALSE"

# var_split="array_both_9"
var_split=$(zcat ${ukb21942_d}/gwas_geno/ukb_geno.var_split.info.tsv.gz | grep -v '#' | awk -v idx=${idx} '(NR == idx){print $1}')


####################################################################
# output files
####################################################################

out=/home/lucytian/data/4_polypred/gwas_geno/${analysis_name}/${covar_str}/${keep_set}/var_split/gwas_${analysis_name}.${var_split}

if [ ! -d $(dirname "${out}") ] ; then mkdir -p $(dirname "${out}") ; fi

####################################################################
# input files
####################################################################

# please see the command below

####################################################################
# main
####################################################################

#if [ "${overwrite}" == "FALSE" ] ; then
#    # check if the output files exist
#    if [ -s "${out}.${out_suffix}.log" ] ; then
#        echo "log file already exists"
#        exit 0
#    fi
#fi

# include the array identity as a covariate if the variant is directly genotyped on both arrays
if [ $(echo ${var_split} | sed -e 's/[0-9]//g') == "array_both_" ] ; then
    var_split_covar=",array"
else
    var_split_covar=""
fi

if [ "${use_tmp_geno}" == "TRUE" ] ; then
    # copy genotype files into cache directory
    bash ${REPODIR}/helpers/cache.pgen.sh ukb_genoHM3

    # Update genotype_pfile variable with the path to the copy of the genotype file in the tmp directory
    genotype_pfile=${tmp_ukb21942_d}/geno/ukb_genoHM3/ukb_genoHM3
else
    genotype_pfile=${ukb21942_d}/geno/ukb_genoHM3/ukb_genoHM3
fi


plink_glm_wrapper () {
    zcat ${ukb21942_d}/gwas_geno/ukb_geno.var_split.tsv.gz \
    | awk -v var_split="${var_split}" '(var_split == $NF){ print $3 }' \
    | plink2 \
        --silent \
        --pfile ${genotype_pfile} vzs \
        --chr 1-22,X,XY,Y,MT \
        --extract /dev/stdin \
        --keep ${ukb21942_d}/sqc/population.20220316/${keep_set}.keep \
        --covar ${ukb21942_d}/sqc/sqc.20220316.tsv.gz \
        --glm zs omit-ref no-x-sex log10 hide-covar skip-invalid-pheno single-prec-cc cc-residualize firth-fallback \
        --covar-variance-standardize \
        $@
}

plink_glm_wrapper_wrapper () {
    # this script-specific phenotype files, phenotype column nums, etc.
    if [ ! -f "${out}.${out_suffix}.${pheno_col_nums}.log" ] ; then
        plink_glm_wrapper \
            --threads ${plink_threads} \
            --memory ${plink_memory} \
            --pheno ${ukb21942_d}/pheno/${analysis_name}.tsv.gz \
            --pheno-col-nums "${pheno_col_nums}" \
            --out "${out}.${out_suffix}.${pheno_col_nums}" \
            $@
    fi
}

if [ "${covar_str}" == "pop_10PCs_no_phequantnorm" ] ; then
    # Population-specific PCs PC1-PC10
    # without applying pheno-quantile-normalize
    # this is useful for effect size comparison beween GWAS and Elastic Net

    plink_glm_wrapper_wrapper \
        --covar-name age,sex,Townsend${var_split_covar},PC1-PC10

elif [ "${covar_str}" == "pop_10PCs" ] ; then
    # Population-specific PCs PC1-PC10

    plink_glm_wrapper_wrapper \
        --pheno-quantile-normalize \
        --covar-name age,sex,Townsend${var_split_covar},PC1-PC10

elif [ "${covar_str}" == "UKB_18PCs" ] ; then
    # UK Biobank PCs PC1-PC18

    plink_glm_wrapper_wrapper \
        --pheno-quantile-normalize \
        --covar-name age,sex,Townsend${var_split_covar},UKB_PC1-UKB_PC18
fi