import pandas as pd
from fastparquet import ParquetFile
import snappy
def snappy_decompress(data, uncompressed_size):
    return snappy.decompress(data)
import h5py
import matplotlib.pyplot as plt
import numpy as np
import random
import seaborn as sns
import matplotlib as mpl

# 1. Transform continuous cv2f features to z-scores
cv2f_features_directory = '/net/bmc-lab5/data/kellis/group/lucytian/data/ENCODE4/chrmbpnet/features/'
cv2f_analysis_directory = 'cv2f_data/cv2f_features/'
snp_information_directory = 'cv2f_data/snp_information/'
chr_features = {}
n_features = []
n_chrom = 22
for i in range(1, n_chrom + 1):
    chr_features[i] = pd.read_csv(cv2f_features_directory + 
                                    f'baseline_annotation_chrmbpnet.{i}.txt.gz', 
                                    sep='\t')
chr_features_full = pd.concat(chr_features.values())
pd.read_csv(cv2f_features_directory + 
            f'baseline_annotation_chrmbpnet.2.txt.gz', 
            sep='\t', nrows = 1)
def z_score(df): return (df-df.mean())/df.std(ddof=0)
h5_zscore_file = f'{cv2f_analysis_directory}cv2f_zscores.h5'
n_prefix_features = 5
features = chr_features_full.columns[n_prefix_features:]
for feature in features:
    chr_feature_zscore = pd.concat([chr_features_full[['CHR', 'BP', 'SNP', 'CM']], 
               z_score(chr_features_full[feature])], axis=1)
    chr_feature_zscore.to_hdf(h5_zscore_file, key=feature, index=False)

# 2. Process SNP data for heatmaps
ipgs_snps = {}
ipgs_snps['ldl'] = pd.read_csv(f'{snp_information_directory}LDL-C_SNPs.tsv', sep='\t')
ipgs_snps['egfr'] = pd.read_csv(f'{snp_information_directory}eGFR_SNPs.tsv', sep='\t')
ipgs_snps['fev1'] = pd.read_csv(f'{snp_information_directory}FEV1_FVC ratio_SNPs.tsv', sep='\t')
ipgs_snps['lymphocyte'] = pd.read_csv(f'{snp_information_directory}Lymphocyte Count_SNPs.tsv', sep='\t')
selected_snps = {}
traits = ['ldl', 'egfr', 'fev1', 'lymphocyte']
for trait in traits:
    significant_snps = list(ipgs_snps[trait].loc[ipgs_snps[trait]['group'] != 'NOT_SIG']['rsID'])
    non_significant_snps = list(ipgs_snps[trait].loc[ipgs_snps[trait]['group'] == 'NOT_SIG']['rsID'])
    selected_snps_from_non_significant = random.sample(non_significant_snps, len(significant_snps))
    selected_snps[trait] = significant_snps + selected_snps_from_non_significant
for trait in traits:
    heatmap_dict = {}
    for feature in features:
        df = pd.read_hdf(h5_zscore_file, key=feature)
        heatmap_dict[feature] = df.loc[df['SNP'].isin(selected_snps[trait])][feature]
    heatmap_df = pd.DataFrame(heatmap_dict)
    overlap_snps = list(df.loc[df['SNP'].isin(selected_snps[trait])]['SNP'])
    heatmap_df.set_index([pd.Index(overlap_snps)], inplace=True)
    heatmap_df.to_hdf(f'{cv2f_analysis_directory}variant_feature_heatmap_data.h5', key=trait, index=False)
shapley = {}
shapley_file_directory = 'shapley_data/tissue/'
shapley['ldl'] = pd.read_csv(
    f'{shapley_file_directory}mvp_afr_0.9_0.01_all_annotation_chrmbpnet_LIVER_lava_ld_ld_maf.shap_values.csv')
shapley['egfr'] = pd.read_csv(
    f'{shapley_file_directory}mvp_afr_0.9_0.01_all_annotation_chrmbpnet_KIDNEY_lava_ld_ld_maf.shap_values.csv')
shapley['fev1'] = pd.read_csv(
    f'{shapley_file_directory}mvp_afr_0.9_0.01_all_annotation_chrmbpnet_LUNG_lava_ld_ld_maf.shap_values.csv')
shapley['lymphocyte'] = pd.read_csv(
    f'{shapley_file_directory}mvp_afr_0.9_0.01_all_annotation_chrmbpnet_BLOOD_lava_ld_ld_maf.shap_values.csv')
for trait in shapley:
    shapley[trait] = shapley[trait][['variable', 'mean_value']].drop_duplicates()

# 3. Generate heatmaps
## 3.1 fev1 main figure
trait = 'fev1'
relevance_set = set()
rsid_to_chrm_pos_map = {}
for i in range(ipgs_snps[trait].shape[0]):
    if ipgs_snps[trait].group[i] != 'NOT_SIG':
        relevance_set.add(ipgs_snps[trait].rsID[i])
        rsid_to_chrm_pos_map[ipgs_snps[trait].rsID[i]] = ipgs_snps[trait]['#ID'][i] + ' (' + ipgs_snps[trait].rsID[i] + ')'
trait_tissue_dict = {'lymphocyte': 'BLOOD', 'ldl': 'LIVER', 'egfr': 'KIDNEY', 'fev1': 'LUNG'}
heatmap_data = pd.read_hdf(f'{cv2f_analysis_directory}variant_feature_heatmap_data.h5', key=trait)
heatmap_data = heatmap_data[list(set(shapley[trait].variable))]
heatmap_data = heatmap_data.loc[list(relevance_set.intersection(heatmap_data.index.values))]
pastel = plt.get_cmap('tab10')
group_category_dict = {}
for variant in heatmap_data.index.values:
    if ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'beta_deviation':
        group_category_dict[variant] = '#FE6100'
    elif ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'cV2F_prioritized':
        group_category_dict[variant] = '#785EF0'
    elif ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'overlap':
        group_category_dict[variant] = '#DC267F'
    else:
        print('error')
bwr = plt.get_cmap('bwr')
prgn = plt.get_cmap('PRGn')
all_betas = []
for variant in heatmap_data.index.values:
    variant_beta_dev = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['cV2F_BETA - baseline_BETA']
    all_betas.append(variant_beta_dev)
min_beta_dev = np.min(all_betas)
max_beta_dev = np.max(all_betas)
beta_dev_dict = {}
beta_dev_scale = 0.004
for variant in heatmap_data.index.values:
    variant_beta_dev = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['cV2F_BETA - baseline_BETA']
    normalized_variant_beta_dev = (variant_beta_dev)/(2*beta_dev_scale) + 0.5
    beta_dev_dict[variant]  = bwr(normalized_variant_beta_dev)
beta_cv2f_dict = {}
variant_beta_scale = 0.007
for variant in heatmap_data.index.values:
    variant_beta = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['BETA_cV2F']
    normalized_variant_beta = variant_beta/(2*variant_beta_scale) + 0.5
    beta_cv2f_dict[variant]  = prgn(normalized_variant_beta)
pu_or = plt.get_cmap('Purples')
all_cv2f = []
for variant in heatmap_data.index.values:
    variant_cv2f = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant][f'cV2F_{trait_tissue_dict[trait]}']
    all_cv2f.append(variant_cv2f)
cv2f_dict = {}
for variant in heatmap_data.index.values:
    variant_cv2f = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant][f'cV2F_{trait_tissue_dict[trait]}']
    normalized_variant_cv2f = float((variant_cv2f))
    cv2f_dict[variant]  = (1-normalized_variant_cv2f/2, 1-normalized_variant_cv2f, 1)
pi_yg = plt.get_cmap('PiYG_r')
shapley_dict_init = {shapley[trait]['variable'].iloc[i]: shapley[trait]['mean_value'].iloc[i] 
                for i in range(len(shapley[trait]['mean_value']))}
shapley_scale = 0.2
shapley_dict = {key: pi_yg(shapley_dict_init[key]/(2*shapley_scale) + 0.5)
              for key in shapley_dict_init}
intersect_keys = [key for key in rsid_to_chrm_pos_map.keys() if key in heatmap_data.index]
intersect_keys_ordered = [key for key in rsid_to_chrm_pos_map.keys() if key in intersect_keys]
reordered_heatmap_data = heatmap_data.loc[intersect_keys_ordered]
labels0 = reordered_heatmap_data.index.values
lut0 = group_category_dict
row_colors0 = pd.Series(labels0).map(lut0)
row_colors0.rename('category', inplace=True)
labels1 = reordered_heatmap_data.index.values
lut1 = beta_dev_dict
row_colors1 = pd.Series(labels1).map(lut1)
row_colors1.rename('beta', inplace=True)
labels2 = reordered_heatmap_data.index.values
lut2 = cv2f_dict
row_colors2 = pd.Series(labels2).map(lut2)
row_colors2.rename('cV2F', inplace=True)
labels3 = reordered_heatmap_data.index.values
lut3 = beta_cv2f_dict
row_colors3 = pd.Series(labels3).map(lut3)
row_colors3.rename('beta_cv2f', inplace=True)
labels_col = reordered_heatmap_data.columns
lut_col = shapley_dict
col_colors = pd.Series(labels_col).map(lut_col)
g = sns.clustermap(reordered_heatmap_data.rename(index=rsid_to_chrm_pos_map), cmap="bwr", center=0, row_colors=[row_colors0, row_colors1, row_colors2, row_colors3], 
                   col_colors=[col_colors], 
                   cbar_kws=dict(orientation='horizontal'), xticklabels=1, yticklabels=1,
                   row_cluster=False,
                  )
g.ax_heatmap.set(xlabel='cV2F Features', ylabel='Variants')
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([x0-.15, .93, g.ax_row_dendrogram.get_position().width, 0.01])
g.ax_cbar.set_title('Functional Signature z-score')
g.fig.set_size_inches(8,18)
plt.savefig(f'heatmap_{trait}.pdf', bbox_inches = "tight")
fig, ax = plt.subplots(figsize=(1, .1))
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=prgn,
                                orientation='horizontal', boundaries=np.linspace(0, 1))
cb1.set_ticks([0,0.5, 1])
cb1.set_ticklabels([-variant_beta_scale, 0, variant_beta_scale])
plt.title(f'cV2F_Beta Colorbar {trait}')
plt.savefig(f'cv2f_beta_cbar_{trait}.pdf', bbox_inches = "tight")

## 3.2 Supplementary figures for other traits
### 3.2.1 eGFR supplementary figure
trait = 'egfr'
relevance_set = set()
rsid_to_chrm_pos_map = {}
for i in range(ipgs_snps[trait].shape[0]):
    if ipgs_snps[trait].group[i] != 'NOT_SIG':
        relevance_set.add(ipgs_snps[trait].rsID[i])
        rsid_to_chrm_pos_map[ipgs_snps[trait].rsID[i]] = ipgs_snps[trait]['#ID'][i] + ' (' + ipgs_snps[trait].rsID[i] + ')'
heatmap_data = pd.read_hdf(f'{cv2f_analysis_directory}variant_feature_heatmap_data.h5', key=trait)
heatmap_data = heatmap_data[list(set(shapley[trait].variable))]
heatmap_data = heatmap_data.loc[list(relevance_set.intersection(heatmap_data.index.values))]
group_category_dict = {}
for variant in heatmap_data.index.values:
    if ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'beta_deviation':
        group_category_dict[variant] = '#FE6100'
    elif ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'cV2F_prioritized':
        group_category_dict[variant] = '#785EF0'
    elif ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'overlap':
        group_category_dict[variant] = '#DC267F'
    else:
        print('error')
all_betas = []
for variant in heatmap_data.index.values:
    variant_beta_dev = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['cV2F_BETA - baseline_BETA']
    all_betas.append(variant_beta_dev)
beta_dev_dict = {}
beta_dev_scale = 0.006
for variant in heatmap_data.index.values:
    variant_beta_dev = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['cV2F_BETA - baseline_BETA']
    normalized_variant_beta_dev = (variant_beta_dev)/(2*beta_dev_scale) + 0.5
    beta_dev_dict[variant]  = bwr(normalized_variant_beta_dev)
varaint_beta_scale = 0.01
beta_cv2f_dict = {}
for variant in heatmap_data.index.values:
    variant_beta = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['BETA_cV2F']
    normalized_variant_beta = (variant_beta)/(2*varaint_beta_scale) + 0.5
    beta_cv2f_dict[variant]  = prgn(normalized_variant_beta)
all_cv2f = []
for variant in heatmap_data.index.values:
    variant_cv2f = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant][f'cV2F_{trait_tissue_dict[trait]}']
    all_cv2f.append(variant_cv2f)
cv2f_dict = {}
for variant in heatmap_data.index.values:
    variant_cv2f = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant][f'cV2F_{trait_tissue_dict[trait]}']
    normalized_variant_cv2f = float((variant_cv2f))
    cv2f_dict[variant]  = (1-normalized_variant_cv2f/2, 1-normalized_variant_cv2f, 1)
shapley_dict_init = {shapley[trait]['variable'].iloc[i]: shapley[trait]['mean_value'].iloc[i] 
                for i in range(len(shapley[trait]['mean_value']))}
shapley_scale = 0.3
shapley_dict = {key: pi_yg(shapley_dict_init[key]/(2*shapley_scale) + 0.5 )
              for key in shapley_dict_init}
intersect_keys = [key for key in rsid_to_chrm_pos_map.keys() if key in heatmap_data.index]
intersect_keys_ordered = [key for key in rsid_to_chrm_pos_map.keys() if key in intersect_keys]
reordered_heatmap_data = heatmap_data.loc[intersect_keys_ordered]
labels0 = reordered_heatmap_data.index.values
lut0 = group_category_dict
row_colors0 = pd.Series(labels0).map(lut0)
row_colors0.rename('category', inplace=True)
labels1 = reordered_heatmap_data.index.values
lut1 = beta_dev_dict
row_colors1 = pd.Series(labels1).map(lut1)
row_colors1.rename('beta', inplace=True)
labels2 = reordered_heatmap_data.index.values
lut2 = cv2f_dict
row_colors2 = pd.Series(labels2).map(lut2)
row_colors2.rename('cV2F', inplace=True)
labels3 = reordered_heatmap_data.index.values
lut3 = beta_cv2f_dict
row_colors3 = pd.Series(labels3).map(lut3)
row_colors3.rename('beta_cv2f', inplace=True)
labels_col = reordered_heatmap_data.columns
lut_col = shapley_dict
col_colors = pd.Series(labels_col).map(lut_col)
g = sns.clustermap(reordered_heatmap_data.rename(index=rsid_to_chrm_pos_map), cmap="bwr", center=0, row_colors=[row_colors0, row_colors1, row_colors2, row_colors3], 
                   col_colors=[col_colors], 
                   cbar_kws=dict(orientation='horizontal'), xticklabels=1, yticklabels=1,
                   row_cluster=False,
                  )
g.ax_heatmap.set(xlabel='cV2F Features', ylabel='Variants')
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([x0-.15, .93, g.ax_row_dendrogram.get_position().width, 0.01])
g.ax_cbar.set_title('Functional Signature z-score')
g.fig.set_size_inches(5,30)
plt.savefig(f'heatmap_{trait}.pdf', bbox_inches = "tight")
fig, ax = plt.subplots(figsize=(1, .1))
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=prgn,
                                orientation='horizontal', boundaries=np.linspace(0, 1))
cb1.set_ticks([0,0.5, 1])
cb1.set_ticklabels([-variant_beta_scale, 0, variant_beta_scale])
plt.title(f'cV2F_Beta Colorbar {trait}')
plt.savefig(f'cv2f_beta_cbar_{trait}.pdf', bbox_inches = "tight")

### 3.2.2 lymphocyte supplementary figure
trait = 'lymphocyte'
relevance_set = set()
rsid_to_chrm_pos_map = {}
for i in range(ipgs_snps[trait].shape[0]):
    if ipgs_snps[trait].group[i] != 'NOT_SIG':
        relevance_set.add(ipgs_snps[trait].rsID[i])
        rsid_to_chrm_pos_map[ipgs_snps[trait].rsID[i]] = ipgs_snps[trait]['#ID'][i] + ' (' + ipgs_snps[trait].rsID[i] + ')'
heatmap_data = pd.read_hdf(f'{cv2f_analysis_directory}variant_feature_heatmap_data.h5', key=trait)
heatmap_data = heatmap_data[list(set(shapley[trait].variable))]
heatmap_data = heatmap_data.loc[list(relevance_set.intersection(heatmap_data.index.values))]
group_category_dict = {}
for variant in heatmap_data.index.values:
    if ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'beta_deviation':
        group_category_dict[variant] = '#FE6100'
    elif ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'cV2F_prioritized':
        group_category_dict[variant] = '#785EF0'
    elif ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'overlap':
        group_category_dict[variant] = '#DC267F'
    else:
        print('error')
all_betas = []
for variant in heatmap_data.index.values:
    variant_beta_dev = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['cV2F_BETA - baseline_BETA']
    all_betas.append(variant_beta_dev)
beta_dev_scale = 0.01
beta_dev_dict = {}
for variant in heatmap_data.index.values:
    variant_beta_dev = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['cV2F_BETA - baseline_BETA']
    normalized_variant_beta_dev = (variant_beta_dev)/(2*beta_dev_scale) + 0.5
    beta_dev_dict[variant]  = bwr(normalized_variant_beta_dev)
variant_beta_scale = 0.03
beta_cv2f_dict = {}
for variant in heatmap_data.index.values:
    variant_beta = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['BETA_cV2F']
    normalized_variant_beta = (variant_beta)/(2*variant_beta_scale) + 0.5
    beta_cv2f_dict[variant]  = prgn(normalized_variant_beta)
all_cv2f = []
for variant in heatmap_data.index.values:
    variant_cv2f = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant][f'cV2F_{trait_tissue_dict[trait]}']
    all_cv2f.append(variant_cv2f)
cv2f_dict = {}
for variant in heatmap_data.index.values:
    variant_cv2f = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant][f'cV2F_{trait_tissue_dict[trait]}']
    normalized_variant_cv2f = float((variant_cv2f))
    cv2f_dict[variant]  = (1-normalized_variant_cv2f/2, 1-normalized_variant_cv2f, 1)
shapley_dict_init = {shapley[trait]['variable'].iloc[i]: shapley[trait]['mean_value'].iloc[i] 
                for i in range(len(shapley[trait]['mean_value']))}
shapley_scale = 0.3
shapley_dict = {key: pi_yg(shapley_dict_init[key]/(2*shapley_scale) + 0.5 )
              for key in shapley_dict_init}
intersect_keys = [key for key in rsid_to_chrm_pos_map.keys() if key in heatmap_data.index]
intersect_keys_ordered = [key for key in rsid_to_chrm_pos_map.keys() if key in intersect_keys]
reordered_heatmap_data = heatmap_data.loc[intersect_keys_ordered]
labels0 = reordered_heatmap_data.index.values
lut0 = group_category_dict
row_colors0 = pd.Series(labels0).map(lut0)
row_colors0.rename('category', inplace=True)
labels1 = reordered_heatmap_data.index.values
lut1 = beta_dev_dict
row_colors1 = pd.Series(labels1).map(lut1)
row_colors1.rename('beta', inplace=True)
labels2 = reordered_heatmap_data.index.values
lut2 = cv2f_dict
row_colors2 = pd.Series(labels2).map(lut2)
row_colors2.rename('cV2F', inplace=True)
labels3 = reordered_heatmap_data.index.values
lut3 = beta_cv2f_dict
row_colors3 = pd.Series(labels3).map(lut3)
row_colors3.rename('beta_cv2f', inplace=True)
labels_col = reordered_heatmap_data.columns
lut_col = shapley_dict
col_colors = pd.Series(labels_col).map(lut_col)
g = sns.clustermap(reordered_heatmap_data.rename(index=rsid_to_chrm_pos_map), cmap="bwr", center=0, row_colors=[row_colors0, row_colors1, row_colors2, row_colors3], 
                   col_colors=[col_colors], 
                   cbar_kws=dict(orientation='horizontal'), xticklabels=1, yticklabels=1,
                   row_cluster=False,
                  )
g.ax_heatmap.set(xlabel='cV2F Features', ylabel='Variants')
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([x0-.15, .93, g.ax_row_dendrogram.get_position().width, 0.01])
g.ax_cbar.set_title('Functional Signature z-score')
g.fig.set_size_inches(15,23)
plt.savefig(f'heatmap_{trait}.pdf', bbox_inches = "tight")
fig, ax = plt.subplots(figsize=(1, .1))
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=prgn,
                                orientation='horizontal', boundaries=np.linspace(0, 1))
cb1.set_ticks([0,0.5, 1])
cb1.set_ticklabels([-variant_beta_scale, 0, variant_beta_scale])
plt.title('cV2F_Beta Colorbar lymphocyte')
plt.savefig('cv2f_beta_cbar_lymphocyte.pdf', bbox_inches = "tight")

## 3.2.3 LDL supplementary figure
trait = 'ldl'
relevance_set = set()
rsid_to_chrm_pos_map = {}
for i in range(ipgs_snps[trait].shape[0]):
    if ipgs_snps[trait].group[i] != 'NOT_SIG':
        relevance_set.add(ipgs_snps[trait].rsID[i])
        rsid_to_chrm_pos_map[ipgs_snps[trait].rsID[i]] = ipgs_snps[trait]['#ID'][i] + ' (' + ipgs_snps[trait].rsID[i] + ')'
heatmap_data = pd.read_hdf(f'{cv2f_analysis_directory}variant_feature_heatmap_data.h5', key=trait)
heatmap_data = heatmap_data[list(set(shapley[trait].variable))]
heatmap_data = heatmap_data.loc[list(relevance_set.intersection(heatmap_data.index.values))]
group_category_dict = {}
for variant in heatmap_data.index.values:
    if ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'beta_deviation':
        group_category_dict[variant] = '#FE6100'
    elif ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'cV2F_prioritized':
        group_category_dict[variant] = '#785EF0'
    elif ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant].group.values[0] == 'overlap':
        group_category_dict[variant] = '#DC267F'
    else:
        print('error')
all_betas = []
for variant in heatmap_data.index.values:
    variant_beta_dev = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['cV2F_BETA - baseline_BETA']
    all_betas.append(variant_beta_dev)
beta_dev_scale = 0.02
beta_dev_dict = {}
for variant in heatmap_data.index.values:
    variant_beta_dev = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['cV2F_BETA - baseline_BETA']
    normalized_variant_beta_dev = (variant_beta_dev)/(2*beta_dev_scale) + 0.5
    beta_dev_dict[variant]  = bwr(normalized_variant_beta_dev)
variant_beta_scale = 0.03
beta_cv2f_dict = {}
for variant in heatmap_data.index.values:
    variant_beta = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant]['BETA_cV2F']
    normalized_variant_beta = (variant_beta)/(2*variant_beta_scale) + 0.5
    beta_cv2f_dict[variant]  = prgn(normalized_variant_beta)
all_cv2f = []
for variant in heatmap_data.index.values:
    variant_cv2f = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant][f'cV2F_{trait_tissue_dict[trait]}']
    all_cv2f.append(variant_cv2f)
cv2f_dict = {}
for variant in heatmap_data.index.values:
    variant_cv2f = ipgs_snps[trait].loc[ipgs_snps[trait]['rsID'] == variant][f'cV2F_{trait_tissue_dict[trait]}']
    normalized_variant_cv2f = float((variant_cv2f))
    cv2f_dict[variant]  = (1-normalized_variant_cv2f/2, 1-normalized_variant_cv2f, 1)
shapley_dict_init = {shapley[trait]['variable'].iloc[i]: shapley[trait]['mean_value'].iloc[i] 
                for i in range(len(shapley[trait]['mean_value']))}
shapley_scale = 0.3
shapley_dict = {key: pi_yg(shapley_dict_init[key]/(2*shapley_scale) + 0.5 )
              for key in shapley_dict_init}
intersect_keys = [key for key in rsid_to_chrm_pos_map.keys() if key in heatmap_data.index]
intersect_keys_ordered = [key for key in rsid_to_chrm_pos_map.keys() if key in intersect_keys]
reordered_heatmap_data = heatmap_data.loc[intersect_keys_ordered]
labels0 = reordered_heatmap_data.index.values
lut0 = group_category_dict
row_colors0 = pd.Series(labels0).map(lut0)
row_colors0.rename('category', inplace=True)
labels1 = reordered_heatmap_data.index.values
lut1 = beta_dev_dict
row_colors1 = pd.Series(labels1).map(lut1)
row_colors1.rename('beta', inplace=True)
labels2 = reordered_heatmap_data.index.values
lut2 = cv2f_dict
row_colors2 = pd.Series(labels2).map(lut2)
row_colors2.rename('cV2F', inplace=True)
labels3 = reordered_heatmap_data.index.values
lut3 = beta_cv2f_dict
row_colors3 = pd.Series(labels3).map(lut3)
row_colors3.rename('beta_cv2f', inplace=True)
labels_col = reordered_heatmap_data.columns
lut_col = shapley_dict
col_colors = pd.Series(labels_col).map(lut_col)
g = sns.clustermap(reordered_heatmap_data.rename(index=rsid_to_chrm_pos_map), cmap="bwr", center=0, row_colors=[row_colors0, row_colors1, row_colors2, row_colors3], 
                   col_colors=[col_colors], 
                   cbar_kws=dict(orientation='horizontal'), xticklabels=1, yticklabels=1,
                   row_cluster=False,
                  )
g.ax_heatmap.set(xlabel='cV2F Features', ylabel='Variants')
x0, _y0, _w, _h = g.cbar_pos
g.ax_cbar.set_position([x0-.15, .93, g.ax_row_dendrogram.get_position().width, 0.01])
g.ax_cbar.set_title('Functional Signature z-score')
g.fig.set_size_inches(10,15)
plt.savefig(f'heatmap_{trait}.pdf', bbox_inches = "tight")
fig, ax = plt.subplots(figsize=(1, .1))
cb1 = mpl.colorbar.ColorbarBase(ax, cmap=prgn,
                                orientation='horizontal', boundaries=np.linspace(0, 1))
cb1.set_ticks([0,0.5, 1])
cb1.set_ticklabels([-variant_beta_scale, 0, variant_beta_scale])
plt.title('cV2F_Beta Colorbar LDL')
plt.savefig('cv2f_beta_cbar_ldl.pdf', bbox_inches = "tight")