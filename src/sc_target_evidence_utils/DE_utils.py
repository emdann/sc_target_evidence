from typing import List
import warnings
import os
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
from datetime import datetime
from sc_target_evidence_utils import preprocessing_utils

import rpy2.robjects.pandas2ri
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import STAP


def _run_glmGamPoi_DE(pbulk_adata, 
                     design = '~ disease',
                    ref_level = None,
                    contrast = 'disease',
                    n_hvgs = None
                    ) -> pd.DataFrame:
    '''
    Run R code for DE analysis with glmGamPoi
    
    Params:
    ------
    - pbulk_adata: anndata object with pseudobulked data
    - design: R formula for design matrix
    - ref_level: reference level for design matrix (default: None, one vs all test)
    - contrast: contrast to test (default: 'disease')

    Returns:
    -------
    - de_res_df: dataframe with DE results for each tested gene
    '''
    
    ## Define R function
    if ref_level is None:
        ref_level = 'NULL'
    else:
        ref_level = f'"{ref_level}"'
        
    if n_hvgs is None:
        n_hvgs = 'NULL'
    else:
        n_hvgs = f'{n_hvgs}'

    glmgampoi_str = f'''
        library(SingleCellExperiment)
        library(glmGamPoi)
        library(scran)
        library(scuttle)

        run_de <- function(args){{
            pbulk_sdata_X <- args[[1]]
            pbulk_sdata_obs <- args[[2]]
            pbulk_sdata_var <- args[[3]]

            sce <- SingleCellExperiment(assays = list(counts = t(pbulk_sdata_X)), colData = pbulk_sdata_obs)
            
            if (!is.null({n_hvgs})){{
                sce <- logNormCounts(sce)

                ## Feature selection w scran (just on test cell types)
                dec <- modelGeneVar(sce)
                hvgs <- getTopHVGs(dec, n = {n_hvgs})
                sce <- sce[hvgs,]
            }}

            ## Fit
            fit <- glm_gp(sce, design = {design}, reference_level = {ref_level}, size_factors = colData(sce)[['size_factors']])

            ## Test 
            de_res <- test_de(fit, contrast = '{contrast}')
            if (!is.null({n_hvgs})){{
                de_res[,'gene_name'] <- pbulk_sdata_var[hvgs,]['feature_name']
                de_res[,'gene_id'] <- pbulk_sdata_var[hvgs,]['feature_id']
            }} else {{
                de_res[,'gene_name'] <- pbulk_sdata_var['feature_name']
                de_res[,'gene_id'] <- pbulk_sdata_var['feature_id']
            }}
            return(de_res)
        }}'''
    
    sc.pp.filter_genes(pbulk_adata, min_cells = 3) # exclude genes expressed in < 3 pseudobulk samples
    
    ## Get anndata components
    pbulk_sdata_X = pbulk_adata.X.toarray().copy()
    pbulk_sdata_obs = pbulk_adata.obs.copy()
    pbulk_sdata_var = pbulk_adata.var.copy()
    
    r_pkg = STAP(glmgampoi_str, "r_pkg")
    # this was needed for the code to run on jhub
    # if you have a different version of rpy2 you may not need these two lines
    rpy2.robjects.pandas2ri.activate()
    rpy2.robjects.numpy2ri.activate()
    
    # Make timestamp to save intermediate files
    dt = datetime.now()
    ts = datetime.timestamp(dt)
    out_filename = f'./DE_results_{ts}.csv'

    # PASS OBJECTS INTO FUNCTION
    args = [pbulk_sdata_X, pbulk_sdata_obs, pbulk_sdata_var]
    de_res_df = r_pkg.run_de(args).to_csvfile(out_filename)
    de_res_df = pd.read_csv(out_filename, index_col=0)
    de_res_df.index = de_res_df['gene_name']
    de_res_df.drop('name', axis=1, inplace=True)
    os.remove(out_filename)
    return(de_res_df)


# def prep_adata_for_de(pbulk_adata):
#     '''Prepare pseudobulk anndata object for DE analysis by removing cell types with insufficient replicates.'''
#     n_replicates = pbulk_adata.obs.value_counts(['high_level_cell_type_ontology_term_id', 'disease']).reset_index()
#     ct_labels = n_replicates[n_replicates[0] >= 3]['high_level_cell_type_ontology_term_id'].unique()
#     pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'].isin(ct_labels)].copy()
#     return(pbulk_adata)


def celltype_marker_targets(pbulk_adata: anndata.AnnData, 
                            min_replicates: int = 3,
                            confounders: List[str] = ['suspension_type', 'assay'],
                            remove_intercept = False,
                            n_hvgs = None
                            ) -> pd.DataFrame:
    '''
    Run DE analysis to identify cell type marker genes.

    Params:
    ------
    - pbulk_adata: anndata object with pseudobulked data
    - min_replicates: minimum number of replicates per cell type (default: 3)
    - confounders: list of confounding variables to include in design matrix (default: ['suspension_type', 'assay'])
    - remove_intercept: indicates whether intercept should be removed (i.e. assuming the reference level in = 0)
    - n_hvgs: number of HVGs to select for DE testing (default: None, no HVG selection)

    Returns:
    -------
    - celltype_de_results: dataframe with DE results for each cell type and tested gene

    '''
    # Keep only normal tissue & Exclude low quality annotation
    pbulk_adata_test = pbulk_adata[
        (pbulk_adata.obs['disease'] == 'normal') & (pbulk_adata.obs['high_level_cell_type_ontology_term_id'] != 'low_quality_annotation')
        ].copy()

    # Exclude cell types with insufficient replicates
    n_replicates = pbulk_adata_test.obs.value_counts(['high_level_cell_type_ontology_term_id', 'disease']).reset_index()
    if 'count' not in n_replicates.columns:
        n_replicates.columns = ['high_level_cell_type_ontology_term_id', 'disease', 'count']
    ct_labels = n_replicates[n_replicates['count'] >= min_replicates]['high_level_cell_type_ontology_term_id'].unique()
    pbulk_adata_test = pbulk_adata_test[pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].isin(ct_labels)].copy()

    pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'] = pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].str.replace(":", "_").astype('category')
    ct_categories = pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].cat.categories.tolist()

    ## Select confounders with multiple values
    for c in confounders:
        if pbulk_adata_test.obs[c].nunique() == 1:
            confounders.remove(c)
    confounders_str = '+'.join(confounders)

    celltype_de_results = pd.DataFrame()
    for ct_term in pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].unique():
        # Reorder cell type categories to include in design matrix
        ct_categories.remove(ct_term)
        ct_categories.append(ct_term)
        pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'] = pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].cat.reorder_categories(ct_categories)
        
        if remove_intercept:
            DE_design = f'~ 0 + {confounders_str} + n_cells + high_level_cell_type_ontology_term_id'
        else:
            DE_design = f'~ {confounders_str} + n_cells + high_level_cell_type_ontology_term_id'
            
        de_results = _run_glmGamPoi_DE(
            pbulk_adata_test,
                design = DE_design,
                ref_level = None,
                contrast = f'high_level_cell_type_ontology_term_id{ct_term}',
                n_hvgs = n_hvgs
                )

        de_results['high_level_cell_type_ontology_term_id'] = ct_term
        celltype_de_results = pd.concat([celltype_de_results, de_results], axis=0)

    return(celltype_de_results)


def disease_marker_targets(pbulk_adata: anndata.AnnData, 
                            min_replicates: int = 3,
                            confounders: List[str] = ['suspension_type', 'assay'],
                            n_hvgs = None
                            ) -> pd.DataFrame:
    '''
    Run DE analysis to identify disease cell type marker genes.

    Params:
    ------
    - pbulk_adata: anndata object with pseudobulked data
    - min_replicates: minimum number of replicates per cell type (default: 3)
    - confounders: list of confounding variables to include in design matrix (default: ['suspension_type', 'assay'])
    - n_hvgs: number of HVGs to select for DE testing (default: None, no HVG selection)

    Returns:
    -------
    - disease_de_results: dataframe with DE results for each cell type and tested gene
    - bulk_de_results: dataframe with DE results for each cell type and tested gene in bulk
    '''
    # Get clean disease label
    pbulk_adata_test = pbulk_adata[
        pbulk_adata.obs['high_level_cell_type_ontology_term_id'] != 'low_quality_annotation'
        ].copy()
    pbulk_adata_test.obs['disease'] = ["_".join(x) for x in pbulk_adata_test.obs['disease'].str.split(" ")]
    pbulk_adata_test.obs['disease'] = [x.replace('/', "_") for x in pbulk_adata_test.obs['disease']]
    pbulk_adata_test.obs['disease'] = [x.replace('-', "_") for x in pbulk_adata_test.obs['disease']]
    disease_term = pbulk_adata_test.obs['disease'][pbulk_adata_test.obs['disease'] !='normal'].unique()[0]

    # Exclude cell types with insufficient replicates in healthy and disease
    n_replicates = pbulk_adata_test.obs.value_counts(['high_level_cell_type_ontology_term_id', 'disease']).reset_index()
    n_replicates = n_replicates[n_replicates['count'] >= min_replicates].groupby(['high_level_cell_type_ontology_term_id']).count()['count']
    ct_labels = n_replicates[n_replicates > 1].index.tolist()
    pbulk_adata_test = pbulk_adata_test[pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].isin(ct_labels)].copy()

    # Reorder disease label categories
    pbulk_adata_test.obs['disease'] = pbulk_adata_test.obs['disease'].astype('category').cat.reorder_categories(['normal', disease_term])

    # Test DE between disease and healthy in cell types
    disease_de_results = pd.DataFrame()
    for ct_term in pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'].unique():
        ad_test = pbulk_adata_test[pbulk_adata_test.obs['high_level_cell_type_ontology_term_id'] == ct_term].copy()
        if ad_test.obs['disease'].unique().shape[0] >= 2:
            ## Select confounders with multiple values
            DE_design = '~' 
            confounding_warning = 0
            for c in confounders:
                if ad_test.obs[c].nunique() > 1:
                    if not _test_matrix_rank(ad_test, 'disease', c):
                        confounding_warning += 1
                        warnings.warn(f'Confounder {c} is collinear with disease label')
                    else:
                        DE_design = DE_design + c + " + "   

            DE_design = DE_design + 'n_cells + disease'

            de_results = _run_glmGamPoi_DE(
                ad_test,
                design = DE_design,
                ref_level = 'normal',
                contrast = f'disease{disease_term}',
                n_hvgs = n_hvgs
                )
            de_results['high_level_cell_type_ontology_term_id'] = ct_term
#             if n_hvgs is not None:
#                 assert de_results.shape[0] == n_hvgs
            de_results['confounder_warning'] = confounding_warning
            disease_de_results = pd.concat([disease_de_results, de_results], axis=0)
            

    # Test DE in "simulated bulk" - aggregating all cell types
    pbulk_adata_test.obs['total_counts'] = pbulk_adata_test.obs['size_factors'].values
    bulk_adata_test = preprocessing_utils.anndata2pseudobulk(
        pbulk_adata_test, 
        group_by=['donor_id', 'disease', 'assay', 'suspension_type'],
        min_ncells=0
        )

    bulk_de_results = _run_glmGamPoi_DE(
                bulk_adata_test,
                design = '~ disease',
                ref_level = 'normal',
                contrast = f'disease{disease_term}',
                n_hvgs = n_hvgs
                )
    return disease_de_results, bulk_de_results


def _test_matrix_rank(pbulk_adata, cov1, cov2):
    model_mat = pd.concat([pd.get_dummies(pbulk_adata.obs[cov1]), pd.get_dummies(pbulk_adata.obs[cov2])], axis=1).astype('int')
    return(np.linalg.matrix_rank(model_mat) >= model_mat.shape[1])
