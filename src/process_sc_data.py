## Get sc data for disease-relevant targets and pseudo-bulk
import pandas as pd
from sc_target_evidence_utils import cxg_utils, plotting_utils, preprocessing_utils, cellontology_utils

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("mondo_id",
                    type=str,
                    help="Mondo ID of disease of interest (format should be MONDO:00000000)")
parser.add_argument("--output_dir",
                    default='/nfs/team205/ed6/bin/sc_target_evidence/data/',
                    type=str,
                    help="Directory to save pseudo-bulked expression objects.")
args = parser.parse_args()

# Parse args
disease_ontology_id = args.mondo_id
output_dir = args.output_dir

## Read cxg cleaned metadata
cxg_metadata = pd.read_csv(output_dir + 'cellxgene_hsapiens_donor_metadata.disease_relevant_annotation.csv', index_col=0)

## Load target data
targets = pd.read_csv(output_dir + 'TargetDiseasePairs_OpenTargets_cellXgeneID_12072023.csv', index_col=0)

## Get target genes
if "_" in disease_ontology_id:
    disease_ontology_id = disease_ontology_id.replace("_", ":")
target_genes = targets[targets['cxg_matched_id'] == disease_ontology_id.replace(":", "_")].targetId.tolist()

## Get target expression in disease-relevant tissue
print("Loading data from cellxgene...")
adata = cxg_utils.get_disease_targets_sc_data(
    target_genes=target_genes, 
    disease_ontology_id=disease_ontology_id,
    cxg_metadata=cxg_metadata,
    keep_disease_relevant_tissue=True,
    keep_normal=True,
    keep_all_genes = False
    )
adata.var_names = adata.var['feature_id'].values
adata.obs['disease_ontology_id'] = disease_ontology_id

## Clean cell ontologies
print("Cleaning cell ontologies...")
graph = cellontology_utils.get_cellontology_graph(output_dir)
ontology_terms = adata.obs["cell_type_ontology_term_id"].unique().tolist()
ct_rename_dict = cellontology_utils.rename_cts_to_high_level(ontology_terms, graph)
adata.obs["high_level_cell_type_ontology_term_id"] = [
    ct_rename_dict[x] if x in ct_rename_dict.keys() else 'low_quality_annotation' for x in adata.obs["cell_type_ontology_term_id"] 
    ]

plotting_utils.plot_celltype_rename(adata.obs, graph, savedir=output_dir + '/plots/')


## Pseudo-bulk by cell ontology and sample
print("Pseudo-bulking...")
pbulk_adata = preprocessing_utils.anndata2pseudobulk(adata, 
                                 group_by=['high_level_cell_type_ontology_term_id', 'assay', 'suspension_type', 'disease', 'donor_id'],
                                 min_ncells=50
                                )

# ## Exclude cell types with data from less than 3 replicates
# n_replicates = pbulk_adata.obs.value_counts(['high_level_cell_type_ontology_term_id', 'disease']).reset_index()
# ct_labels = n_replicates[n_replicates[0] >= 3]['high_level_cell_type_ontology_term_id'].unique()
# pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'].isin(ct_labels)].copy()

## Save 
print("Saving...")
pbulk_adata.write_h5ad(output_dir + f'cellxgene_targets_{disease_ontology_id.replace(":","_")}.pbulk_all_OT_targets.h5ad')
