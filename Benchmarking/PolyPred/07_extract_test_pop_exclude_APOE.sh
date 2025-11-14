ukb21942_d='/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942'
keep_set='Afr_test'


plink2 \
  --silent \
  --pfile /home/lucytian/group/data/geno/ukb_genoHM3_exclude_v2 vzs \
  --keep ${ukb21942_d}/sqc/population.20220316/${keep_set}.keep \
  --extract snp_extract.txt \
  --make-bed \
  --out APOE_pop/${keep_set}


plink \
  --bfile APOE_pop/${keep_set} \
  --update-name update_name.txt \
  --make-bed \
  --out APOE_pop/${keep_set}.rsid


plink \
  --bfile APOE_pop/${keep_set}.rsid \
  --a1-allele a1_effect.txt 1 2 \
  --make-bed \
  --out APOE_pop/${keep_set}.a1tag