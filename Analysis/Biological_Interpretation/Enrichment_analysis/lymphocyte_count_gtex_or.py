import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

data_d = '/home/lucytian/group/data/ENCODE4/chrmbpnet/features'
out_d = '/home/lucytian/data/10_Enrichment_analysis/binary_or_table'
os.makedirs(out_d, exist_ok=True)

trait = pd.read_csv('Lymphocyte Count.groupings.tsv', sep='\t')

df1 = pd.read_csv(f'{data_d}/baseline_annotation_chrmbpnet.1.txt.gz', sep='\t', compression='gzip')

cols = [
    i for i in df1.columns
    if (
        'blood' in i.lower() or
        'k562' in i.lower() or
        'gm12878' in i.lower() or
        'bld' in i.lower()
    )
]

binary = [
    i for i in cols
    if (
        i.startswith('cCREv4') or
        i.startswith('ABC_ENCODE') or
        i.startswith('GTEx_Flagship')
    )
]

def compute_or_ABC_vs_D(df, annot_col):
    signal = df['group'].isin(['Group C', 'Group A', 'Group B'])
    background = df['group'] == 'Group D'

    a = ((signal) & (df[annot_col] == 1)).sum()
    b = ((signal) & (df[annot_col] == 0)).sum()
    c = ((background) & (df[annot_col] == 1)).sum()
    d = ((background) & (df[annot_col] == 0)).sum()

    table = [[a, b], [c, d]]
    OR, p = fisher_exact(table)

    return pd.DataFrame([{
        'annotation': annot_col,
        'a': a, 'b': b, 'c': c, 'd': d,
        'odds_ratio': OR,
        'p_value': p
    }])

dfs = [df1[['CHR', 'BP', 'SNP', 'CM', 'base'] + binary]]

for i in range(2, 23):
    df = pd.read_csv(f'{data_d}/baseline_annotation_chrmbpnet.{i}.txt.gz', sep='\t', compression='gzip')
    df = df[['CHR', 'BP', 'SNP', 'CM', 'base'] + binary]
    dfs.append(df)

all_df = pd.concat(dfs, ignore_index=True)

merged = all_df.merge(trait, left_on='SNP', right_on='rsID')

results = []
for b in binary:
    res = compute_or_ABC_vs_D(merged, annot_col=b)
    results.append(res)

binary_or = pd.concat(results, ignore_index=True)
binary_or.to_csv(f'{out_d}/Lymphocyte_Count_table.tsv', sep='\t', index=False)