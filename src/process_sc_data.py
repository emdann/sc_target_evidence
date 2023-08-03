## Get sc data for disease-relevant targets and pseudo-bulk
import cellxgene_census
import tiledbsoma as soma
from tiledbsoma.experiment_query import X_as_series
import pandas as pd
import numpy as np
import obonet

import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import scipy
import anndata

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("mondo_id",
                    type=str,
                    help="Mondo ID of disease of interest (format should be MONDO:00000000)")
parser.add_argument("--output_dir",
                    default='/nfs/team205/ed6/bin/sc_target_evidence/data/',
                    type=str,
                    help="Directory to save pseudo-bulked expression objects.")
# parser.add_argument("--all_genes",
#                     default=False,
#                     type=bool,
#                     help="Boolean indicating whether all genes should be pulled")
args = parser.parse_args()

# Parse args
disease_ontology_id = args.mondo_id
output_dir = args.output_dir

### ---- Cellxgene census utils ----- ###

def _get_cxg_datasets(cxg_metadata, disease_ontology_id, keep_disease_relevant_tissue, keep_normal):
    filt_metadata = cxg_metadata[cxg_metadata['disease_ontology_id'] == disease_ontology_id]
    
    if keep_disease_relevant_tissue:
        filt_metadata = filt_metadata[filt_metadata['tissue_general'] == filt_metadata['disease_relevant_tissue']].copy()

    if keep_normal:
        normal_metadata = cxg_metadata[ (cxg_metadata['disease_ontology_id'] == 'PATO:0000461') & (cxg_metadata['tissue_general'] == filt_metadata['disease_relevant_tissue'].iloc[0]) ]
        filt_metadata = pd.concat([normal_metadata, filt_metadata], axis=0)
    
    return filt_metadata


def get_disease_targets_sc_data(
    disease_ontology_id,
    target_genes, 
    cxg_metadata,
    keep_disease_relevant_tissue = True,
    keep_normal = True,
    keep_all_genes = False
    ):
    '''
    Params:
    -------
    
    - disease_ontology_id: disease of interest (MONDO ID) - downloads data for disease and healthy tissue for disease-relevant tissue
    - target_genes: list of target genes for the disease
    - keep_disease_relevant_tissue: whether to filter only samples from the disease-relevant tissue
    - keep_normal: whether to include data from normal tissue
    - keep_all_genes: boolean indicating whether to filter target genes only or not
    
    Returns:
    --------
    adata: AnnData object for targets of interest
    ''' 
    
    OBS_COLS = [
            "assay", 
            "tissue_general", 
            "suspension_type", 
            "disease", 
            'donor_id', 
            'cell_type_ontology_term_id',
            'cell_type'
    ]
    
    # Make string of features 
    target_genes_str = '[' + ', '.join(f"'{item}'" for item in target_genes) + ']'
    if keep_all_genes:
        var_filter_str = None
    else:
        var_filter_str = f"feature_id in {target_genes_str}"
    
    # Subset to relevant datasets
    sc_metadata = _get_cxg_datasets(cxg_metadata, disease_ontology_id, keep_disease_relevant_tissue, keep_normal)
    dataset_ids = sc_metadata['dataset_id'].unique()
    tissue_ids = sc_metadata['tissue_general_original'].unique()
    disease_ids = sc_metadata['disease_ontology_id_original'].unique()
    tissue_ids_str = '[' + ', '.join(f"'{item}'" for item in tissue_ids) + ']'
    disease_ids_str = '[' + ', '.join(f"'{item}'" for item in disease_ids) + ']'
    dataset_ids_str = '[' + ', '.join(f"'{item}'" for item in dataset_ids) + ']'
    obs_filter_str = f"dataset_id in {dataset_ids_str} and disease_ontology_term_id in {disease_ids_str} and tissue_general in {tissue_ids_str}"
    
    with cellxgene_census.open_soma(census_version="2023-05-15") as census:
        print("Getting anndata")
        # Get expression of target genes as anndata
        adata = cellxgene_census.get_anndata(
                census = census,
                organism = "Homo sapiens",
                var_value_filter = var_filter_str,
                obs_value_filter = obs_filter_str,
                column_names = {"obs": OBS_COLS}
            )
        
        # Store total counts per cell (library size for DE analysis) 
        print("Getting total counts x cell")
        hsapiens = census["census_data"]["homo_sapiens"]
        if not keep_all_genes:
            with hsapiens.axis_query(
                measurement_name="RNA",
                obs_query=soma.AxisQuery(value_filter=obs_filter_str),
            ) as query:
                obs_df = query.obs().concat().to_pandas().set_index("soma_joinid")
                n_obs = len(obs_df)

                raw_sum = np.zeros((n_obs,), dtype=np.float64)  # accumulate the sum of expression

                # query.X() returns an iterator of pyarrow.Table, with X data in COO format.
                # You can request an indexer from the query that will map it to positional indices
                indexer = query.indexer
                for arrow_tbl in query.X("raw").tables():
                    obs_dim = indexer.by_obs(arrow_tbl["soma_dim_0"])
                    data = arrow_tbl["soma_data"]
                    np.add.at(raw_sum, obs_dim, data)
        else:
            raw_sum = np.array(adata.X.sum(1)).flatten()

        adata.obs['total_counts'] = raw_sum
        
    return(adata)


### --- Utils for cell ontology --- ###

obo_file = output_dir + "cl.obo"  # downloaded from http://obofoundry.org/ontology/cl.html


# Load the ontology from the OBO file
graph = obonet.read_obo(obo_file)

def get_ancestors(term_id, graph, levels=3):
    ancestors = set()
    stack = [(term_id, 0)]  # Using a stack to keep track of terms and their levels

    while stack:
        current_term, current_level = stack.pop()
        ancestors.add(current_term)

        if current_level < levels:
            parents = [parent for parent in graph.successors(current_term) if parent.startswith("CL:")]
            stack.extend([(parent, current_level + 1) for parent in parents])

    return ancestors

def get_CL_ancestors(term_id, graph):
    ancestors = get_ancestors(term_id, graph)
    return list(ancestors)

def rename_cts_to_high_level(adata, graph):
    original_terms = adata.obs['cell_type_ontology_term_id'].unique().tolist()
    # Exclude VERY high level annotations (e.g. native cell)
    keep_terms = []
    for o in original_terms:
        if len(get_ancestors(o, graph)) > 5:
            keep_terms.append(o)
    
    new_terms = original_terms.copy()
    for i,term in enumerate(original_terms):
        # Check if cell_type ontology term is the most high level in dataset  
        upstream_CL = get_CL_ancestors(term, graph)
        parent_terms = np.intersect1d(upstream_CL, keep_terms)
        if len(parent_terms) > 0:
            new_terms[i] = parent_terms[0]
    ct_rename_dict = dict(zip(original_terms, new_terms))
    
    # Rename to higher level if it exists
    while any([ct_rename_dict[v] != v for v in ct_rename_dict.values()]):
        for k,v in ct_rename_dict.items():
            if ct_rename_dict[v] != v:
                ct_rename_dict[k] = ct_rename_dict[v]
                
    adata.obs['high_level_cell_type_ontology_term_id'] = [ct_rename_dict[x] for x in adata.obs['cell_type_ontology_term_id']]
    return(ct_rename_dict)

### --- Pseudobulking utils --- ###

def anndata2pseudobulk(adata, group_by, 
                       agg="s", 
                       min_ncells = 3,
                       use_layer = None
                      ):
    '''
    Do pseudo-bulking of raw counts in anndata object
    
    Params:
    ------
    - adata: the anndata object
    - group_by: list of obs columns to use for aggregation
    - agg: "s" for sum (if adata.X are counts), "m" for mean (if adata.X are log-counts)
    - min_ncells: minimum number of cells to keep pseudobulk sample (default: 10)
    - use_layer: which layer to use for pseudobulking (default: None, use adata.X)
    
    Returns:
    -------
    AnnData object of same n_vars as adata, but pseudobulked obs
    '''    
    if use_layer is None:
        X = adata.X.copy()
    else:
        X = adata.layers[use_layer].copy()
        
    ## Make obs for pseudobulk
    pseudobulk_obs = adata.obs[group_by].drop_duplicates()
    pseudobulk_obs = pseudobulk_obs[group_by].astype("str")
    pseudobulk_obs.index = pseudobulk_obs[group_by].agg("-".join, axis=1)
    
    ## Add column to obs assigning cells to pseudobulk samples
    adata.obs[group_by] = adata.obs[group_by].astype("str")
    adata.obs["pseudobulk_sample"] = adata.obs[group_by].agg("-".join, axis=1)
    
    ## Sum counts from same sample
    sample_dummies = pd.get_dummies(adata.obs["pseudobulk_sample"])[pseudobulk_obs.index].values
    sample_dummies = scipy.sparse.csr_matrix(sample_dummies)
    pseudobulk_X = X.T.dot(sample_dummies)
    
    ## Compute library sizes (for DE)
    total_counts = adata.obs['total_counts'].values
    size_factors = scipy.sparse.csr_matrix(total_counts).dot(sample_dummies)
    pseudobulk_obs['size_factors'] = size_factors.toarray().flatten()
    
    ## Make new anndata object
    pseudobulk_adata = anndata.AnnData(pseudobulk_X.T, obs=pseudobulk_obs, var=adata.var)
    
    ## Add number of cells to obs 
    n_cells = adata.obs.groupby('pseudobulk_sample').count().iloc[:,0]
    n_cells.name = "n_cells"
    pseudobulk_adata.obs = pd.concat([pseudobulk_adata.obs, n_cells], axis=1)
    
    ## Filter obs by number of cells threshold
    pseudobulk_adata = pseudobulk_adata[pseudobulk_adata.obs['n_cells'] >= min_ncells].copy()
    return(pseudobulk_adata)



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
adata = get_disease_targets_sc_data(
    target_genes=target_genes, 
    disease_ontology_id=disease_ontology_id,
    cxg_metadata=cxg_metadata,
    keep_disease_relevant_tissue=True,
    keep_normal=True,
    keep_all_genes = False
    )
adata.var_names = adata.var['feature_id'].values

## Clean cell ontologies
print("Cleaning cell ontologies...")
ct_rename_dict = rename_cts_to_high_level(adata, graph)

## Pseudo-bulk by cell ontology and sample
print("Pseudo-bulking...")
pbulk_adata = anndata2pseudobulk(adata, 
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
