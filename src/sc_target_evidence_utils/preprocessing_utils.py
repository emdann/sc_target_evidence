from typing import List, Union, Tuple, Dict, Optional

import pandas as pd
import scipy
import anndata

# from sc_target_evidence_utils import cellontology_utils


def anndata2pseudobulk(
        adata: anndata.AnnData, 
        group_by: List[str], 
        min_ncells: int = 3,
        use_layer: Union[str, None] = None
        ):
    '''
    Pseudo-bulk raw counts in anndata object.
    
    Params:
    ------
    - adata: the anndata object
    - group_by: list of obs columns to use for aggregation
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

def clean_disease_name(disease_name):
    return "_".join(disease_name.split(' ')).replace('/', "_").replace('-', "_")


# def prep_pbulk_adata(pbulk_adata, data_dir):
#         graph = cellontology_utils.get_cellontology_graph(data_dir)
#         pbulk_adata.obs['high_level_cell_type'] = [f'{cellontology_utils.ontology2name(x, graph)} ({x})' for x in pbulk_adata.obs['high_level_cell_type_ontology_term_id'].tolist()]
#         pbulk_adata.obs['sample_id'] = ['-'.join(x[1:]) for x in pbulk_adata.obs_names.str.split("-")]

#         ## Run PCA
#         pbulk_adata.layers['counts'] = pbulk_adata.X.copy()
#         cpms = scipy.sparse.csr_matrix(pbulk_adata.X.T / pbulk_adata.obs['size_factors'].values.flatten()) * 1000000
#         pbulk_adata.X  = np.log1p(cpms).T
#         sc.tl.pca(pbulk_adata, n_comps=20, use_highly_variable=False, random_state=42)

#         ## Add dendrogram
#         sc.tl.dendrogram(pbulk_adata, groupby='high_level_cell_type', use_rep='X_pca')
#         pbulk_adata.X  = pbulk_adata.layers['counts'].copy()