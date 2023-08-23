import pandas as pd
import numpy as np

def DE2evidence_celltype(
    ct_res: pd.DataFrame,
    signif_thresh: float = 0.1,
    lfc_thresh = 1,
    ct_fraction_thresh = 0.5
    ) -> dict:
    '''Get targets with cell type marker evidence from DE analysis results.
    
    Parameters
    ----------
    ct_res : pd.DataFrame
        DE results from cell type marker analysis (output of `DE_utils.celltype_marker_targets`).
    signif_thresh : float, optional
        Significance threshold for FDR (default: 0.1).
    lfc_thresh : float, optional
        Log fold change threshold (default: 1).
    ct_fraction_thresh : float, optional
        Max fraction of cell types in which a gene must be significantly over-expressed to be considered a marker (default: 0.5).

    Returns
    -------
    evidence_dict : dict
        Dictionary containing list of target genes with cell type marker evidence.
    '''
    evidence_dict = {}

    tot_celltypes = ct_res['high_level_cell_type_ontology_term_id'].nunique()
    # Are the target genes cell type markers (significantly over-expressed in a subset of cell types)
    ct_marker_evidence_df = ct_res[
        (ct_res['adj_pval'] < signif_thresh) & (ct_res['lfc'] > lfc_thresh)
        ].groupby('gene_id').size()/tot_celltypes
    evidence_dict['ct_marker_evidence'] = ct_marker_evidence_df.index[ct_marker_evidence_df < ct_fraction_thresh].values
    return(evidence_dict)

def DE2evidence_disease(
    disease_res: pd.DataFrame,
    bulk_res: pd.DataFrame,
    signif_thresh: float = 0.1,
    lfc_thresh = 1,
    lfc_group = 'all'
    ) -> dict:
    '''Get targets with sc evidence from DE analysis results.
    
    Parameters
    ----------
    disease_res : pd.DataFrame
        DE results from disease analysis (output of `DE_utils.disease_marker_targets`).
    disease_res : pd.DataFrame
        DE results from disease analysis on bulk tissue (output of `DE_utils.disease_marker_targets`).
    signif_thresh : float, optional
        Significance threshold for FDR (default: 0.1).
    lfc_thresh : float, optional
        Log fold change threshold (default: 1).
    lfc_group : bool, optional
        whether to consider both up- and down-regulated genes (one of 'all', 'positive', 'negative', default: all).
        
    Returns
    -------
    evidence_dict : dict
        Dictionary containing list of target genes with disease marker evidence.
        - disease_evidence: genes that are differentially expressed between disease and control in a subset of cell types
        - disease_ct_evidence: genes that are differentially expressed between disease and control in a subset of cell types and not in bulk tissue
        - bulk_disease_evidence: genes that are differentially expressed between disease and control in bulk tissue
    '''
    
    evidence_dict = {}

    # Are the target genes differentially expressed between disease and control in a subset of cell types (and not detectable in bulk)?
    evidence_dict['bulk_disease_evidence'] = bulk_res[
        _filter_res(bulk_res, signif_thresh, lfc_thresh, lfc_group)
        ].gene_id.unique()

    evidence_dict['disease_evidence'] = disease_res[
        _filter_res(disease_res, signif_thresh, lfc_thresh, lfc_group)
        ].gene_id.unique()
    
    evidence_dict['disease_ct_evidence'] = np.setdiff1d(evidence_dict['disease_evidence'], evidence_dict['bulk_disease_evidence'])
    return(evidence_dict)

def _filter_res(de_res, signif_thresh, lfc_thresh, lfc_group):
    '''util for making filter'''
    if lfc_group == 'all':
        return (de_res['adj_pval'] < signif_thresh) & (np.abs(de_res['lfc']) > lfc_thresh)
    elif lfc_group == 'positive':
        return (de_res['adj_pval'] < signif_thresh) & (de_res['lfc'] > lfc_thresh)        
    elif lfc_group == 'negative':
        return (de_res['adj_pval'] < signif_thresh) & (de_res['lfc'] < -lfc_thresh)
    else:
        raise ValueError('lfc_group must be one of "all", "positive", "negative"')