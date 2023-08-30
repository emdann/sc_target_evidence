## Get sc data for disease-relevant targets and pseudo-bulk
import pandas as pd
import anndata
import json
from sc_target_evidence_utils import cxg_utils, plotting_utils, preprocessing_utils, cellontology_utils

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("mondo_id",
                    type=str,
                    help="Mondo ID of disease of interest (format should be MONDO:00000000)")
parser.add_argument("--keep_all_genes",
                    type=bool,
                    help="Whether to keep all the genes (within possible universes)")
parser.add_argument("--output_dir",
                    default='/nfs/team205/ed6/bin/sc_target_evidence/data/',
                    type=str,
                    help="Directory to save pseudo-bulked expression objects.")
args = parser.parse_args()


def clean_disease(pbulk_adata):
    '''Uniform disease naming.'''        
    DISEASE_RENAME = {'cardiomyopathy':[
        'arrhythmogenic right ventricular cardiomyopathy',
        'dilated cardiomyopathy',
        'non-compaction cardiomyopathy'],
    'renal cell carcinoma':['chromophobe renal cell carcinoma', 'clear cell renal carcinoma'],
    'colorectal cancer': ['colorectal cancer', 'colorectal neoplasm'],
        'non-small cell lung carcinoma':['lung large cell carcinoma', 'non-small cell lung carcinoma']
        }

    disease_rename_rev = {x:k for k,v in DISEASE_RENAME.items() for x in v }

    pbulk_adata.obs['disease_name_original'] = pbulk_adata.obs['disease'].copy()
    pbulk_adata.obs['disease'] = [disease_rename_rev[x] if x in disease_rename_rev.keys() else x for x in pbulk_adata.obs.disease]

def _download_dataset(dataset_id, tissue_ids, disease_ids, ct_rename_dict):
    '''Download and preprocess one dataset'''
    adata = cxg_utils.get_disease_targets_dataset(
        dataset_ids = [dataset_id], tissue_ids = tissue_ids, disease_ids= disease_ids,
        target_genes=keep_genes, 
        keep_all_genes = False
        )
    adata.var_names = adata.var['feature_id'].values
    adata.obs['disease_ontology_id'] = disease_ontology_id

    # Exclude blacklisted assays
    assay_blacklist = [
        'BD Rhapsody Targeted mRNA',
        'STRT-seq',
        'inDrop'
        ]

    adata = adata[~adata.obs['assay'].isin(assay_blacklist)].copy()
    adata = adata[adata.obs['is_primary_data']].copy()

    adata.obs["high_level_cell_type_ontology_term_id"] = [
        ct_rename_dict[x] if x in ct_rename_dict.keys() else 'low_quality_annotation' for x in adata.obs["cell_type_ontology_term_id"] 
        ]

    ## Pseudo-bulk by cell ontology and sample
    print("Pseudo-bulking...")
    pbulk_adata = preprocessing_utils.anndata2pseudobulk(adata, 
                                     group_by=['high_level_cell_type_ontology_term_id', 'assay', 'suspension_type', 'disease', 'donor_id'],
                                     min_ncells=5       
                                    )
    return(pbulk_adata)

# Parse args
disease_ontology_id = args.mondo_id
output_dir = args.output_dir
keep_all_genes = args.keep_all_genes

## Read cxg cleaned metadata
cxg_metadata = pd.read_csv(output_dir + 'cellxgene_hsapiens_donor_metadata.disease_relevant_annotation.csv', index_col=0)

## Load target data
targets = pd.read_csv(output_dir + 'TargetDiseasePairs_OpenTargets_cellXgeneID_12072023.csv', index_col=0)

## Get target genes
if "_" in disease_ontology_id:
    disease_ontology_id = disease_ontology_id.replace("_", ":")
target_genes = targets[targets['cxg_matched_id'] == disease_ontology_id.replace(":", "_")].targetId.tolist()

## Get list of all possible genes you might need 
json_file =   output_dir + 'target_universe_dict.json'
with open(json_file, "r") as json_file:
    universe_dict = json.load(json_file)
    
all_genes_ls = sum([x for x in universe_dict.values()], [])
all_genes_ls.extend(target_genes)
all_genes_ls = list(set(all_genes_ls))

if keep_all_genes:
    keep_genes = all_genes_ls
else:
    keep_genes = target_genes

print(f'# genes: {len(keep_genes)}')

## Get datasets of interest
sc_metadata = cxg_utils._get_cxg_datasets(cxg_metadata, disease_ontology_id, keep_disease_relevant_tissue=True, keep_normal=True)
dataset_ids = sc_metadata['dataset_id'].unique()
tissue_ids = sc_metadata['tissue_general_original'].unique()
disease_names = sc_metadata['disease'].unique()
disease_ids = sc_metadata['disease_ontology_id_original'].unique()

## Get annotations from cell metadata
cxg_cell_metadata = pd.read_csv('/nfs/team205/ed6/data/cellxgene_hsapiens_cell_metadata.csv')
cxg_cell_metadata = cxg_cell_metadata[(cxg_cell_metadata.tissue_general.isin(tissue_ids)) & (cxg_cell_metadata.disease.isin(disease_names)) & (cxg_cell_metadata.dataset_id.isin(dataset_ids))].copy()
cxg_cell_metadata['disease_ontology_id'] = disease_ontology_id

# Clean cell ontologies
print("Cleaning cell ontologies...")
graph = cellontology_utils.get_cellontology_graph(output_dir)

# Exclude cell types with less than 5 cells
ct_counts = cxg_cell_metadata["cell_type_ontology_term_id"].value_counts()
ontology_terms = ct_counts[ct_counts > 5].index.tolist()
ct_rename_dict = cellontology_utils.rename_cts_to_high_level(ontology_terms, graph)
cxg_cell_metadata["high_level_cell_type_ontology_term_id"] = [
    ct_rename_dict[x] if x in ct_rename_dict.keys() else 'low_quality_annotation' for x in cxg_cell_metadata["cell_type_ontology_term_id"] 
    ]
plotting_utils.plot_celltype_rename(cxg_cell_metadata, disease_ontology_id, graph, savedir=output_dir + '/plots/')

## Get target expression in disease-relevant tissue
pbulk_ls = []
for dataset in dataset_ids:
    print(f"Processing dataset {dataset}...")
    pbulk_dataset = _download_dataset(dataset, tissue_ids, disease_ids, ct_rename_dict)
    pbulk_ls.append(pbulk_dataset)
    
pbulk_adata = anndata.concat(pbulk_ls)

clean_disease(pbulk_adata)
pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'] != 'low_quality_annotation'].copy()
pbulk_adata.obs['high_level_cell_type'] = [f'{cellontology_utils.ontology2name(x, graph)} ({x})' for x in pbulk_adata.obs['high_level_cell_type_ontology_term_id'].tolist()]
pbulk_adata.obs['sample_id'] = ['-'.join(x[1:]) for x in pbulk_adata.obs_names.str.split("-")]
pbulk_adata.obs['disease_ontology_id'] = disease_ontology_id

## Save 
print("Saving...")
if not keep_all_genes:
    pbulk_adata.write_h5ad(output_dir + f'cellxgene_targets_{disease_ontology_id.replace(":","_")}.pbulk_all_OT_targets.h5ad')
else:
    pbulk_adata.write_h5ad(output_dir + f'cellxgene_targets_{disease_ontology_id.replace(":","_")}.pbulk_all_genes.h5ad')
