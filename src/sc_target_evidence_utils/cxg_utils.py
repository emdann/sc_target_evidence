from typing import List
from anndata import AnnData

import cellxgene_census
import tiledbsoma as soma
import pandas as pd
import numpy as np

def _get_cxg_datasets(
        cxg_metadata: pd.DataFrame, 
        disease_ontology_id: str, 
        keep_disease_relevant_tissue: bool, 
        keep_normal: bool) -> AnnData:
    '''
    Get dataset_ids from curated metadata table for a given disease_ontology_id.
    '''
    filt_metadata = cxg_metadata[cxg_metadata['disease_ontology_id'] == disease_ontology_id]
    
    if keep_disease_relevant_tissue:
        filt_metadata = filt_metadata[filt_metadata['tissue_general'] == filt_metadata['disease_relevant_tissue']].copy()

    if keep_normal:
        normal_metadata = cxg_metadata[ (cxg_metadata['disease_ontology_id'] == 'PATO:0000461') & (cxg_metadata['tissue_general'] == filt_metadata['disease_relevant_tissue'].iloc[0]) ]
        filt_metadata = pd.concat([normal_metadata, filt_metadata], axis=0)
    
    return filt_metadata


def get_disease_targets_sc_data(
    disease_ontology_id: str,
    target_genes: List[str], 
    cxg_metadata: pd.DataFrame,
    keep_disease_relevant_tissue: bool = True,
    keep_normal: bool = True,
    keep_all_genes: bool = False
    ):
    '''
    Download data for disease of interest from cellxgene database, using datasets from curated metadata.

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
            'cell_type',
            'is_primary_data'
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
    obs_filter_str = f"dataset_id in {dataset_ids_str} and disease_ontology_term_id in {disease_ids_str} and tissue_general in {tissue_ids_str} and is_primary_data == True"
    
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

def get_disease_targets_dataset(
    dataset_ids,
    tissue_ids,
    disease_ids,
    target_genes: List[str], 
    keep_all_genes: bool = False
    ):
    '''
    Download data for disease of interest from cellxgene database.
    
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
            'cell_type',
            'is_primary_data'
    ]
    
    # Make string of features 
    target_genes_str = '[' + ', '.join(f"'{item}'" for item in target_genes) + ']'
    if keep_all_genes:
        var_filter_str = None
    else:
        var_filter_str = f"feature_id in {target_genes_str}"
    
    # Make query string
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