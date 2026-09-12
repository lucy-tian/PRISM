plink2 --pfile synthetic_v1_Afr_chr22only \
       --indep-pairwise 50 5 0.2 \
       --out pruned_data

plink2 --pfile synthetic_v1_Afr_chr22only \
       --extract pruned_data.prune.in \
       --pca approx 10 \
       --out afr_pca