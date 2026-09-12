TRAITS=("INI1003063" "INI20030780" "INI30120" "INI50030700")

# Paths
SUMSTAT_DIR="sumstat"
OUTDIR="output"


make -p 

for TRAIT in "${TRAITS[@]}"; do
    echo "Merging effects for $TRAIT..."
    
    mkdir -p "PRS/${TRAIT}"
    LOGFILE="PRS/${TRAIT}/polypred.log"
    ERRFILE="PRS/${TRAIT}/polypred.err"
    
    python /home/lucytian/data/0_SOFTWARE/polyfun/polypred.py \
        --estimate-mixweights \
        --betas linear_input/${TRAIT}_sumstat.txt.gz,linear_input/${TRAIT}_finemapp.txt.gz \
        --pheno test_pheno/${TRAIT}.pheno.txt \
        --output-prefix PRS/${TRAIT} \
        --plink-exe ~/data/0_SOFTWARE/plink \
        test_pop/Afr_test.a1tag.bed >"$LOGFILE" 2>"$ERRFILE"
done