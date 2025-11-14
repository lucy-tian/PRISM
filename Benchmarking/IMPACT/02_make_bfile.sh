plink2 --pfile /net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/geno/ukb_genoHM3/ukb_genoHM3 vzs \
       --keep /net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/sqc/population.20220316/Afr.keep \
       --make-bed --out genotype/Afr


plink2 --pfile /home/lucytian/group/data/geno/ukb_genoHM3_exclude_v2 vzs \
       --keep /net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/sqc/population.20220316/Afr.keep \
       --make-bed --out genotype/Afr_exclude_APOE