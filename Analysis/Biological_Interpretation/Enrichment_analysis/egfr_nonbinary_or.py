#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import pyreadr
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import os

data_d = '/home/lucytian/group/data/ENCODE4/chrmbpnet/features'
beta_d = '/home/lucytian/data/1_Single_Cell_PRS/2_cV2F/mvp_afr_0.9_0.01/20231221/406k_geno_v2_UKB_18PCs/fit_w_val'

df = pd.read_csv(f'{data_d}/baseline_annotation_chrmbpnet.1.txt.gz', sep='\t', compression='gzip')

trait = pd.read_csv('eGFR.groupings.tsv', sep='\t')

cols = [
    i for i in df.columns
    if (
        'kid' in i.lower()
    )
]

prefixes = ('cCREv4', 'ABC_ENCODE', 'GTEx_Flagship')

non_binary = [
    i for i in cols
    if not i.startswith(prefixes)
]

def compute_logistic_enrichment_ABC_vs_D(df, annot_col):
    df = df.copy()
    df = df.dropna(subset=[annot_col])
    
    df['is_signal'] = df['group'].isin(['Group A', 'Group B', 'Group C']).astype(int)
    
    x = (df[annot_col] - df[annot_col].mean()) / df[annot_col].std()
    
    X = sm.add_constant(x)
    y = df['is_signal']
    
    model = sm.Logit(y, X).fit(disp=0)
    
    beta = model.params.iloc[1]
    OR = np.exp(beta)
    p = model.pvalues.iloc[1]
    
    return pd.DataFrame([{
        'annotation': annot_col,
        'beta': beta,
        'odds_ratio_per_SD': OR,
        'p_value': p
    }])

dfs = [df[['CHR', 'BP', 'SNP', 'CM', 'base'] + non_binary]]
for i in np.arange(2, 23):
    df = pd.read_csv(f'{data_d}/baseline_annotation_chrmbpnet.{i}.txt.gz', sep='\t', compression='gzip')
    df = df[['CHR', 'BP', 'SNP', 'CM', 'base'] + non_binary]
    dfs.append(df)

all_df = pd.concat(dfs, ignore_index=True)

merged = all_df.merge(trait, left_on='SNP', right_on='rsID')

results = []
for nb in non_binary:
    res = compute_logistic_enrichment_ABC_vs_D(merged, annot_col=nb)
    results.append(res)

or_results = pd.concat(results, ignore_index=True)

out_d = '/home/lucytian/data/10_Enrichment_analysis/nonbinary_or_table'
os.makedirs(out_d, exist_ok=True)

or_results.to_csv(f'{out_d}/eGFR_table.tsv', sep='\t', index=False)