import os
import anndata
import pandas as pd

import rpy2.robjects.pandas2ri
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import STAP


def _run_glmGamPoi_DE(pbulk_adata, 
                     design = '~ disease',
                    ref_level = 'normal',
                    contrast = 'disease'
                    ):
    '''
    Run R code for DE analysis with glmGamPoi
    
    Params:
    ------
    - pbulk_adata: anndata object with pseudobulked data
    
    '''
    
    ## Define R function
    glmgampoi_str = f'''
        library(SingleCellExperiment)
        library(glmGamPoi)
        library(scran)
        library(scater)

        run_de <- function(args){{
            pbulk_sdata_X <- args[[1]]
            pbulk_sdata_obs <- args[[2]]
            pbulk_sdata_var <- args[[3]]

            sce <- SingleCellExperiment(assays = list(counts = t(pbulk_sdata_X)), colData = pbulk_sdata_obs)

            ## Fit
            fit <- glm_gp(sce, design = {design}, reference_level = '{ref_level}', size_factors = colData(sce)[['size_factors']])

            ## Test 
            de_res <- test_de(fit, contrast = '{contrast}')    
            de_res[,'gene_name'] <- pbulk_sdata_var['feature_name']
            de_res[,'gene_id'] <- pbulk_sdata_var['feature_id']
            return(de_res)
        }}'''
    
    ## Get anndata components
    pbulk_sdata_X = pbulk_adata.X.toarray().copy()
    pbulk_sdata_obs = pbulk_adata.obs.copy()
    pbulk_sdata_var = pbulk_adata.var.copy()
    
    r_pkg = STAP(glmgampoi_str, "r_pkg")
    # this was needed for the code to run on jhub
    # if you have a different version of rpy2 you may not need these two lines
    rpy2.robjects.pandas2ri.activate()
    rpy2.robjects.numpy2ri.activate()
    
    # PASS OBJECTS INTO FUNCTION
    args = [pbulk_sdata_X, pbulk_sdata_obs, pbulk_sdata_var]
    de_res_df = r_pkg.run_de(args).to_csvfile('./DE_results.csv')
    de_res_df = pd.read_csv('./DE_results.csv', index_col=0)
    de_res_df.index = de_res_df['gene_name']
    de_res_df.drop('name', 1, inplace=True)
    os.remove('./DE_results.csv')
    return(de_res_df)


# def prep_adata_for_de(pbulk_adata):
#     '''Prepare pseudobulk anndata object for DE analysis by removing cell types with insufficient replicates.'''
#     n_replicates = pbulk_adata.obs.value_counts(['high_level_cell_type_ontology_term_id', 'disease']).reset_index()
#     ct_labels = n_replicates[n_replicates[0] >= 3]['high_level_cell_type_ontology_term_id'].unique()
#     pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'].isin(ct_labels)].copy()
#     return(pbulk_adata)


def celltype_marker_targets(pbulk_adata: anndata.AnnData, min_replicates: int =3):
    # Keep only normal tissue
    pbulk_adata_test = pbulk_adata[pbulk_adata.obs['disease'] == 'normal'].copy()

    # Exclude cell types with insufficient replicates
    n_replicates = pbulk_adata.obs.value_counts(['high_level_cell_type_ontology_term_id', 'disease']).reset_index()
    ct_labels = n_replicates[n_replicates[0] >= min_replicates]['high_level_cell_type_ontology_term_id'].unique()
    pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'].isin(ct_labels)].copy()

    pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'] = pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].str.replace(":", "_").astype('category')
    ct_categories = pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].cat.categories.tolist()

    celltype_de_results = pd.DataFrame()
    for ct_term in pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].unique():
        # Reorder cell type categories to include in design matrix
        ct_categories.remove(ct_term)
        ct_categories.append(ct_term)
        pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'] = pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].cat.reorder_categories(ct_categories)

        de_results = _run_glmGamPoi_DE(
            pbulk_adata_test,
            design = '~ assay + suspension_type + high_level_cell_type_ontology_term_id', # + n_cells
            ref_level = 'normal',
            contrast = f'high_level_cell_type_ontology_term_id{ct_term}'
            )
        de_results['high_level_cell_type_ontology_term_id'] = ct_term
        celltype_de_results = pd.concat([celltype_de_results, de_results], axis=0)


