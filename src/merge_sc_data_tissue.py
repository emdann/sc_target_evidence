'''Merge single datasets downloaded from CxG into one obsect per tested'''
import os,sys
import pandas as pd
import numpy as np
import scanpy as sc
import anndata

import argparse
# parse arguments t, data_dir, target_table_file
parser = argparse.ArgumentParser()
parser.add_argument('t', type=str, help='Tissue name')
parser.add_argument('--data_dir', default='/oak/stanford/groups/pritch/users/emma/bin/sc_target_evidence/data/', type=str, help='Directory where data is stored')
parser.add_argument('--target_table_file', default = 'filtered_nelson_disease_relevant_tissues_06172024.csv', type=str, help='Name of target table csv file (stored in data_dir)')
args = parser.parse_args()

t = args.t
data_dir = args.data_dir
target_table_file = f'{data_dir}/{args.target_table_file}'

# Read tissue2dataset mapping
tissue2dataset = pd.read_csv(f'{data_dir}/tissue2dataset_id.csv', index_col=0)
all_tissues = tissue2dataset.tissue_general.unique()

## Read clinical targets
clin_order = ['Preclinical', 'Phase I', 'Phase II', 'Phase III', 'Launched']
target_table = pd.read_csv(target_table_file, index_col=0)
target_table['combined_max_phase'] = target_table['combined_max_phase'].astype('category').cat.reorder_categories(clin_order, ordered=True)
try:
    tissue_targets_df = target_table[~target_table[f'dtr_{t.replace("-", "_")}'].isna()].copy()
except KeyError:
    if t == 'eye':
        tissue_targets_df = target_table[~target_table[f'dtr_retina'].isna()].copy()
    else:
        raise ValueError(f'No target data for {t}')

# Read all anndata objects 
tissue_dataset_ids = tissue2dataset[tissue2dataset.tissue_general == t].tissue_dataset.tolist()
file_paths = [f'../data/cellxgene_targets_{id}.pbulk_all_genes.h5ad' for id in tissue_dataset_ids]
adata_ls = [sc.read_h5ad(f) for f in file_paths]

# Make common var table
n_datasets_per_var = pd.Series(sum([a.var_names.tolist() for a in adata_ls], [])).value_counts()
n_datasets_per_var.name = 'n_datasets'
all_var = pd.concat([a.var[['feature_id', 'feature_name']] for a in adata_ls]).drop_duplicates()
all_var = pd.concat([all_var, n_datasets_per_var], axis=1)

# Merge on all measured genes
adata = anndata.concat(adata_ls, join='outer', merge='same')
assert sum([a.n_obs for a in adata_ls]) == adata.n_obs
adata.var = all_var.loc[adata.var_names].copy()

# Count non-zero entries per gene across datasets 
nonzero_counts = np.diff(adata.X.tocsc().indptr)
adata.var['nnz'] = nonzero_counts

# Filter genes expressed in at least 1% of samples (nnz >= 0.1 * n_obs) and at least 25% of the datasets 
tot_datasets = len(adata_ls)
adata.var['is_detected'] = (adata.var['n_datasets'] >= (0.25 * tot_datasets)) & (adata.var['nnz'] >= (0.01 * adata.n_obs))
adata = adata[:, adata.var['is_detected']].copy()

# Add targets from Nelson table
target_max_phase = tissue_targets_df.groupby('target')['combined_max_phase'].max().to_frame().reset_index().rename({'target':'feature_name'}, axis=1)
adata.var['is_tissue_target'] = adata.var['feature_name'].isin(tissue_targets_df.target.unique())
adata.var = pd.merge(adata.var, target_max_phase, how='left')

# Save adata as h5ad
adata.write(f'{data_dir}/cellxgene_targets_{t}.pbulk_all_genes.h5ad')
