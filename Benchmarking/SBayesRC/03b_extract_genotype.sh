ukb21942_d='/net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942'
genotype_pfile=${ukb21942_d}/geno/ukb_genoHM3/ukb_genoHM3


out_dir='genotype'
mkdir -p $out_dir 

  plink2 \
  --pfile ${genotype_pfile} vzs \
  --extract snp_extract.txt \
  --update-name update_name.txt \
  --keep Afr_val_test.keep \
  --make-bed \
  --out ${out_dir}/Afr_val_test
  
  
  plink2 \
  --pfile /home/lucytian/group/data/geno/ukb_genoHM3_exclude_v2 vzs \
  --extract snp_extract.txt \
  --update-name update_name.txt \
  --keep Afr_val_test.keep \
  --make-bed \
  --out ${out_dir}/Afr_val_test_exclude_APOE
  