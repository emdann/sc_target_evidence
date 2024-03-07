import pandas as pd
from scipy.stats.contingency import odds_ratio
from scipy.stats import fisher_exact

def get_OR(
    targets_disease_df, 
    evidence_col,
    clinical_status_col,
    gene_universe = None
    ):
    '''
    Compute Odds-ratio for association between omic evidence and target clinical status.

    Parameters
    ----------
    targets_disease_df : pd.DataFrame
        DataFrame of target-disease pairs containing information on omic evidence and clinical status in different columns.
    evidence_col : str, optional
        Name of column containing boolean or integer for omic evidence (default: 'ct_marker_evidence').
    clinical_status_col : str, optional
        Name of column containing boolean or integer for clinical status (default: 'is_effective').
    gene_universe : list, optional
        List of genes to consider in the universe. If None (default), use all the genes in the input DataFrame.
    
    Returns
    -------
    or_df : pd.DataFrame
        DataFrame containing odds-ratio and other statistics for each target-disease pair.
    '''
    if gene_universe is not None:
        ## Filter genes in universe
        targets_disease_df = targets_disease_df[targets_disease_df['gene_id'].isin(gene_universe)].copy()

    if sum(targets_disease_df[evidence_col]) == 0:
        raise ValueError(f'no gene with {evidence_col} evidence')
    if sum(targets_disease_df[clinical_status_col]) == 0:
        raise ValueError(f'no gene in status {clinical_status_col}')
    
    # Get contingency table
    cont = pd.crosstab(
        targets_disease_df[evidence_col], 
        targets_disease_df[clinical_status_col]
    ).loc[[1,0],[1,0]].T
    
    if gene_universe is not None:
        ## Add missing unsuccessful genes for each disease_ontology_id
        n_targets_by_disease = targets_disease_df.groupby('disease_ontology_id').size() 
        missing_genes_by_disease = len(gene_universe) - n_targets_by_disease
        cont.loc[0,0] = cont.loc[0,0] + missing_genes_by_disease.sum()
        assert cont.sum().sum() == (targets_disease_df['disease_ontology_id'].nunique() * len(gene_universe))
        
    # Compute Odds ratios
    or_res = odds_ratio(cont.loc[[1,0],[1, 0]].T)
    or_enrichment = or_res.statistic
    or_ci = or_res.confidence_interval(confidence_level=0.95)
    pval = fisher_exact(cont, alternative='greater').pvalue

    # Get number of targets with evidence and clinical status
    n_targets_evidence = cont.loc[1,1]
    n_success = targets_disease_df[clinical_status_col].sum()
    if gene_universe is not None:
        n_insuccess = cont.sum().sum() - n_success 
    else:
        n_insuccess = targets_disease_df[clinical_status_col].shape[0] - n_success 
    n_supported = targets_disease_df[evidence_col].sum()
    
    or_df = pd.DataFrame([or_enrichment, or_ci.low, or_ci.high, pval, n_success,n_insuccess, n_targets_evidence, n_supported], index=['odds_ratio', 'ci_low', 'ci_high', 'pval', 'n_success', 'n_insuccess', 'n_supported_approved','n_supported']).T
    or_df['evidence'] = evidence_col
    or_df['clinical_status'] = clinical_status_col
    return(or_df)

def compute_grouped_OR(
    targets_evidence_all, 
    group_by,
    evidence_cols = ['disease_ct_evidence', 'ct_marker_evidence', 'disease_evidence','bulk_disease_evidence', 'has_genetic_support'],
    clinical_status_cols = ['is_druggable', 'is_safe', 'is_effective', 'is_approved'],
    gene_universe = None
    ):
    '''
    Compute OddsRatios stratifying target-disease pairs by a given column.

    Parameters
    ----------
    targets_evidence_all : pd.DataFrame
        DataFrame of target-disease pairs containing information on omic evidence and clinical status in different columns.
    group_by : str, optional
        Name of column to stratify target-disease pairs (e.g. 'disease_ontology_id').
    evidence_cols : list, optional
        List of columns containing boolean or integer for omic evidence (default: ['disease_ct_evidence', 'ct_marker_evidence', 'disease_evidence','bulk_disease_evidence', 'has_genetic_support']).
    clinical_status_cols : list, optional
        List of columns containing boolean or integer for clinical status (default: ['is_druggable', 'is_safe', 'is_effective', 'is_approved']).
    
    Returns
    -------
    or_df_all : pd.DataFrame
        DataFrame containing odds-ratio and other statistics for each target-disease pair for each group.
    '''
    or_df_all = pd.DataFrame()
    for g in targets_evidence_all[group_by].unique():
        targets_evidence_df = targets_evidence_all[targets_evidence_all[group_by] == g].copy()
        for ev in evidence_cols:
            for status in clinical_status_cols:
                if targets_evidence_df[status].sum() > 0:
                    try:
                        or_df = get_OR(targets_evidence_df, ev, status, gene_universe)
                        or_df[group_by] = g
                        or_df_all = pd.concat([or_df_all, or_df], axis=0)
                    except ValueError:
                        continue
                else:
                    continue
    return(or_df_all)