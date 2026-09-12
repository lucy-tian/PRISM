
ukb21942_d='/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942'
genotype_pfile=${ukb21942_d}/geno/ukb_genoHM3/ukb_genoHM3
keep_set='Afr_test'



plink2 \
  --silent \
  --pfile ${genotype_pfile} vzs \
  --keep ${ukb21942_d}/sqc/population.20220316/${keep_set}.keep \
  --extract snp_extract.txt \
  --make-bed \
  --out test_pop/${keep_set}


plink \
  --bfile test_pop/${keep_set} \
  --update-name update_name.txt \
  --make-bed \
  --out test_pop/${keep_set}.rsid


plink \
  --bfile test_pop/${keep_set}.rsid \
  --a1-allele a1_effect.txt 1 2 \
  --make-bed \
  --out test_pop/${keep_set}.a1tag