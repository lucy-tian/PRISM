traits=(INI30120 INI20030780 INI1003063)
thresholds=('0_1' '0_03' '0_01' '0_003' '0_001' '0_0003' '0_0001' '3e-05' '1e-05')

mkdir -p log

for t in "${traits[@]}"; do
  for p in "${thresholds[@]}"; do
    sbatch -p kellis \
      -o log/${t}_${p}.%j.out \
      -e log/${t}_${p}.%j.err \
      --wrap="sh 04b_thresholding.sh $t $p"
  done
done
