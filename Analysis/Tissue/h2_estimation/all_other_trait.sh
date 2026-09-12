#!/bin/bash

ldsc_repo_d="/home/lucytian/data/0_SOFTWARE/ldsc"
weights_d="/home/lucytian/data/1_Single_Cell_PRS/ldsc-examples/1000G_Phase3_weights_hm3_no_MHC"
af_d="/home/lucytian/data/1_Single_Cell_PRS/ldsc-examples/1000G_Phase3_frq"
save_d="enrichment"
ld_d="ldscore"

if [[ ${1} == "INI30120" ]]; then
    munged_stat_d="/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/gwas_geno/blood_count/UKB_18PCs/WB/LDSC"
elif [[ ${1} == "INI1003063" ]]; then
    munged_stat_d="/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/gwas_geno/spirometry/UKB_18PCs/WB/LDSC"
else
    munged_stat_d="/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/gwas_geno/blood_biochemistry/UKB_18PCs/WB/LDSC"
fi

python ${ldsc_repo_d}/ldsc.py \
       --h2 ${munged_stat_d}/${1}.sumstats.gz \
       --ref-ld-chr ${ld_d}/1000G_EUR_geno. \
       --w-ld-chr ${weights_d}/weights.hm3_noMHC. \
       --overlap-annot \
       --frqfile-chr ${af_d}/1000G.EUR.QC. \
       --out ${save_d}/${1} \
       --print-coefficients \
       --print-delete-vals