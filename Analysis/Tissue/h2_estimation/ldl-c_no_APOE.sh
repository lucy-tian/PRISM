ldsc_repo_d="/home/lucytian/data/0_SOFTWARE/ldsc"
save_d="/home/lucytian/data"
alkes_LDSCORE_d="/net/bmc-lab5/data/kellis/group/tanigawa/data/alkesgroup/LDSCORE"


python ${ldsc_repo_d}/ldsc.py \
       --h2 /home/lucytian/data/INI20030780_no_APOE.sumstats.gz \
       --ref-ld-chr ${alkes_LDSCORE_d}/baseline_v1.2/baseline. \
       --w-ld-chr ${alkes_LDSCORE_d}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
       --frqfile-chr ${alkes_LDSCORE_d}/1000G_Phase3_frq/1000G.EUR.QC.
       --overlap-annot \
       --out ${save_d}/INI20030780_no_APOE \
       --print-coefficients \
       --print-delete-vals