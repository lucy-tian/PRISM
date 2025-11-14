TRAITS=("INI1003063" "INI20030780" "INI30120" "INI50030700")

# Paths
SUMSTAT_DIR="sumstat"
OUTDIR="output"



for TRAIT in "${TRAITS[@]}"; do
    echo "Merging effects for $TRAIT..."
    
    mkdir -p "PRS/${TRAIT}"
    LOGFILE="PRS/${TRAIT}/polypred.PRS.log"
    ERRFILE="PRS/${TRAIT}/polypred.PRS.err"
    
    python /home/lucytian/data/0_SOFTWARE/polyfun/polypred.py \
        --predict \
        --betas linear_input/${TRAIT}_sumstat.txt.gz,linear_input/${TRAIT}_finemapp.txt.gz \
        --mixweights-prefix PRS/${TRAIT} \
        --output-prefix PRS/${TRAIT}.pred \
        --plink-exe ~/data/0_SOFTWARE/plink \
        test_pop/Afr_test.a1tag.bed >"$LOGFILE" 2>"$ERRFILE"
done