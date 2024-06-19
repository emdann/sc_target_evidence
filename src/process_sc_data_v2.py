## Get sc data for disease-relevant targets and pseudo-bulk
import pandas as pd
import json
import cellxgene_census
import numpy as np
from sc_target_evidence_utils import preprocessing_utils, cellontology_utils

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("tissue_dataset_id",
                    type=str,
                    help="ID of dataset to download")
parser.add_argument("--census_version",
                    default="2023-12-15",
                    type=str,
                    help="cellxgene census version for reproducibility")
parser.add_argument("--filters_table",
                    default='/oak/stanford/groups/pritch/users/emma/bin/sc_target_evidence/data/cellxgene_filters.ct_marker_evidence.csv',
                    type=str,
                    help="Path to table matching tissue_dataset_id to cellxgene census filters - computed in cellxgene_census_stats_v3.ipynb")
parser.add_argument("--data_dir",
                    default='/oak/stanford/groups/pritch/users/emma/bin/sc_target_evidence/data/',
                    type=str,
                    help="Directory to save pseudo-bulked expression objects.")
args = parser.parse_args()


def clean_disease(pbulk_adata):
    '''Uniform disease naming.'''        
    DISEASE_RENAME = {'cardiomyopathy':[
    'arrhythmogenic right ventricular cardiomyopathy',
     'dilated cardiomyopathy',
     'non-compaction cardiomyopathy'],
  'renal cell carcinoma':['chromophobe renal cell carcinoma', 'clear cell renal carcinoma', 'nonpapillary renal cell carcinoma',
                         'kidney oncocytoma', ],
  'colorectal cancer': ['colorectal cancer', 'colorectal neoplasm'],
    'non-small cell lung carcinoma':['lung large cell carcinoma', 'non-small cell lung carcinoma']
    }

    disease_rename_rev = {x:k for k,v in DISEASE_RENAME.items() for x in v }

    pbulk_adata.obs['disease_name_original'] = pbulk_adata.obs['disease'].copy()
    pbulk_adata.obs['disease'] = [disease_rename_rev[x] if x in disease_rename_rev.keys() else x for x in pbulk_adata.obs.disease]


# Parse args
tissue_dataset_id = args.tissue_dataset_id
data_dir = args.data_dir
census_version = args.census_version
filters_table_path = args.filters_table

## Read cxg filter for obs values
filters_df = pd.read_csv(filters_table_path, index_col=0)
obs_filter_str = filters_df.loc[tissue_dataset_id, 'filter_string']

# ## Load target data
# targets = pd.read_csv(data_dir + 'TargetDiseasePairs_OpenTargets_cellXgeneID_12072023.csv', index_col=0)
# ## Get target genes
# if "_" in disease_ontology_id:
#     disease_ontology_id = disease_ontology_id.replace("_", ":")
# target_genes = targets[targets['cxg_matched_id'] == disease_ontology_id.replace(":", "_")].targetId.tolist()
# ## Get list of all possible genes you might need 
# json_file =   data_dir + 'target_universe_dict.json'
# with open(json_file, "r") as json_file:
#     universe_dict = json.load(json_file)
# all_genes_ls = sum([x for x in universe_dict.values()], [])
# all_genes_ls.extend(target_genes)
# all_genes_ls = list(set(all_genes_ls))
# if keep_all_genes:
#     keep_genes = all_genes_ls
# else:
#     keep_genes = target_genes

## Query cellxgene census
OBS_COLS = [
            "assay", 
            "tissue_general", 
            "suspension_type", 
            "disease", 
            'donor_id', 
            'cell_type_ontology_term_id',
            'cell_type',
            'is_primary_data',
            'raw_sum',
            'nnz'
    ]
var_filter_str = None

with cellxgene_census.open_soma(census_version=census_version) as census:
    print("Getting anndata")
    # Get expression of target genes as anndata
    adata = cellxgene_census.get_anndata(
            census = census,
            organism = "Homo sapiens",
            var_value_filter = var_filter_str,
            obs_value_filter = obs_filter_str,
            column_names = {"obs": OBS_COLS}
        )

adata.obs.rename({'raw_sum':'total_counts'}, axis=1, inplace=True)
adata.var_names = adata.var['feature_id'].values
adata.obs['tissue_dataset_id'] = tissue_dataset_id
assert adata.obs['tissue_general'].nunique() == 1, 'More tissues than expected'
tissue = adata.obs['tissue_general'].unique().tolist()[0]

# Some checks
assay_blacklist = [
    'BD Rhapsody Targeted mRNA',
    'STRT-seq',
    'inDrop'
    ]

assert all(adata.obs['is_primary_data']), 'Primary data is not filtered'
assert not any(adata.obs['assay'].isin(assay_blacklist)), 'Blacklisted assays are not filtered'

## Clean cell ontologies
print("Cleaning cell ontologies...")
graph = cellontology_utils.get_cellontology_graph(data_dir)
with open(f'{data_dir}/celltype_harmonization_dict.json', 'r') as json_file:
    ct_harmonize_dict = json.load(json_file)
ct_rename_dict = ct_harmonize_dict[tissue]

# Exclude cell types with less than 5 cells
adata.obs["high_level_cell_type_ontology_term_id"] = [
    ct_rename_dict[x] if x in ct_rename_dict.keys() else 'low_quality_annotation' for x in adata.obs["cell_type_ontology_term_id"] 
    ]

## Pseudo-bulk by cell ontology and sample
print("Pseudo-bulking...")
pbulk_adata = preprocessing_utils.anndata2pseudobulk(adata, 
                                 group_by=['high_level_cell_type_ontology_term_id', 'assay', 'suspension_type', 'disease','donor_id'],
                                 min_ncells=5       
                                )

clean_disease(pbulk_adata)
pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'] != 'low_quality_annotation'].copy()
pbulk_adata.obs['high_level_cell_type'] = [f'{cellontology_utils.ontology2name(x, graph)} ({x})' for x in pbulk_adata.obs['high_level_cell_type_ontology_term_id'].tolist()]
pbulk_adata.obs['sample_id'] = ['-'.join(x[1:]) for x in pbulk_adata.obs_names.str.split("-")]
pbulk_adata.obs['tissue_dataset_id'] = tissue_dataset_id

# Filter out genes that are not expressed in dataset
pbulk_adata.var['sum_counts'] = np.array(pbulk_adata.X.sum(0)).flatten()
pbulk_adata = pbulk_adata[:, pbulk_adata.var['sum_counts'] > 0].copy()

## Save 
print("Saving...")
pbulk_adata.write_h5ad(data_dir + f'cellxgene_targets_{tissue_dataset_id.replace(":","_")}.pbulk_all_genes.h5ad')
