for q in 0 0.2 0.4 0.6 0.8 1; do
    q_scaled=$(awk -v q="$q" 'BEGIN {printf "%.0f", q * 100}')
    out_dir="results/q${q_scaled}"

    mkdir -p "$out_dir"

    Rscript 2a.snpnet.R "$out_dir" "$q"
    Rscript 2b.result_curate.R "$out_dir"
done