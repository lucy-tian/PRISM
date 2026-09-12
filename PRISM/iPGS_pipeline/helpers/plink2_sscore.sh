#!/bin/bash
set -beEuo pipefail

# this script applies plink2 --sscore and compute
# PGS for each individual

SRCNAME=$(readlink -f "${0}")
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
beta_f=$(readlink -f $1)
gcount_file=$2

ml load plink2

out="${beta_f%.BETAs.tsv.gz}"
tmp_beta_f=${tmp_dir}/$(basename ${beta_f})
colname=$(basename $(dirname ${beta_f}))

# Drop covariate terms and replace the header name from BETA to trait ID
# and save it as a file in a temp directory
zcat ${beta_f} \
| awk -v FS='\t' 'NR==1 || $2 != ""' \
| sed -e "s/BETA/${colname}/g" > ${tmp_beta_f}

# Taking advantage of the sparsity of the PGS models,
# we specify the sets of variants with non-zero BETAs using
# plink2's --extract option.
# We extract variant IDs using awk and pass that information to
# standard input (/dev/stdin) via pipe
cat ${tmp_beta_f} | awk '(NR>1){print $1}' \
| plink2 \
--threads 6 \
--memory 40000 \
--pfile ${ukb21942_d}/geno/ukb_genoHM3 vzs \
--read-freq ${gcount_file} \
--extract /dev/stdin \
--out ${out} \
--score \
  ${tmp_beta_f} \
  1 2 3 header-read zs \
  cols=maybefid,maybesid,phenos,nallele,dosagesum,scoreavgs,denom,scoresums

mv ${out}.log ${out}.sscore.log
