## Data and output pointers

Data folder: 
 - In repo: `sc_target_evidence/data/`
 - On server (includes git-ignored large files): `/nfs/team205/ed6/bin/sc_target_evidence/data/`


### Figshare data 
- **suppl_table_cxg_sample_metadata.post_qc.csv** - sample-level metadata for data from CellXGene collections used in study 
- **suppl_table_diseases.csv** - table of diseases available in CellXGene considered in study, with disease-relevant tissue annotations and reasons to exclude
- **suppl_table_disease_target_evidence.csv** - table of target-disease pairs with indication of clinical success and omic support  
- **suppl_table_odds_ratios.all.csv** - Results of association analysis between omic support (`evidence`) and clinical success (`clinical status`) across diseases
- **suppl_table_odds_ratios.disease.csv** - Results of association analysis between omic support (`evidence`) and clinical success (`clinical status`) by disease

### Metadata

- **cellxgene_hsapiens_donor_metadata.csv** - donor-level metadata for all adult studies in CellxGene collections (see `cellxgene_census_stats.ipynb` for download)
- **cellxgene_hsapiens_donor_metadata.disease_relevant_annotation.csv** - curated donor-level metadata for all adult studies in CxG collection (including grouping of diseases, grouping of tissues, annotation of disease-relevant tissues)
- [on server] **TargetDiseasePairs_OpenTargets_cellXgeneID_12072023.clean.csv** - target-disease pairs in OpenTargets, with clinical status annotation based on OT drug evidence score.
- **target_universe_dict.json** - dictionary of all the tested gene target universes used throughout analysis (built in [`analysis/make_target_universe.py`](https://github.com/emdann/sc_target_evidence/blob/master/analysis/make_target_universe.py))
- **all_CellOntologies.csv** - table of all cell ontology terms considered across diseases and tissues
- **all_diseaseIDs.txt** - list of all disease IDs considered for analysis (*note: not all diseases end in final result*)

### scRNA-seq data

- [on server] **cellxgene_targets_{MONDO_ID}.pbulk_all_genes.h5ad** - AnnData object with aggregated expression of all available genes in disease-relevant tissue and pseudo-bulk.

#### Diagnostic plots

Plot folders: `sc_target_evidence/data/plots/{disease_id}_{disease_name}`

- **cellxgene_{disease_id}.celltype_harmonization.*** - confusion table of original cell ontology annotations and uniformed ontology annotations. Heatmap color and number in cells denotes the number of cells for each category.
- **cellxgene_targets_{disease_id}.n_cells_boxplot.*** - boxplot of numbers of cells per sample and cell type in healthy and disease tissue
- **cellxgene_targets_{disease_id}.target_expression.*** - heatmap of log-normalized expression of a sample of drug targets for the disease
- **cellxgene_targets_{disease_id}.celltype_distribution.*** - confusion table of assignment of uniformed cell type ontology to each donor, to check differences in cell type distribution across donors/datasets/diseases

#### Differential expression analysis outputs

DE analysis results: 
- [on server] **DE_celltype_{disease_id}.hvgs.csv** - results of differential expression analysis between cell types, ran on top 7500 highly variable genes (output of [`sc_target_evidence_utils.DE_utils.celltype_marker_targets`](https://github.com/emdann/sc_target_evidence/blob/9e9658d9443f6f1ca642f008ffc18e847982c476/src/sc_target_evidence_utils/DE_utils.py#L155)). Results for all tested cell types are shown. 
- [on server] **DE_diseasecelltype_{disease_id}.hvgs.csv** - results of differential expression analysis between disease and healthy cells in each cell types, ran on top 7500 highly variable genes (output of [`sc_target_evidence_utils.DE_utils.disease_marker_targets`](https://github.com/emdann/sc_target_evidence/blob/9e9658d9443f6f1ca642f008ffc18e847982c476/src/sc_target_evidence_utils/DE_utils.py#L224)). Results for all tested cell types are shown. 
- [on server] **DE_diseasebulk_{disease_id}.hvgs.csv** - results of differential expression analysis between disease and healthy cells in bulked samples (ignoring cell type), ran on top 7500 highly variable genes (output of [`sc_target_evidence_utils.DE_utils.disease_marker_targets`](https://github.com/emdann/sc_target_evidence/blob/9e9658d9443f6f1ca642f008ffc18e847982c476/src/sc_target_evidence_utils/DE_utils.py#L224)).

Single-cell evidence outputs: 

- **suppl_table_disease_target_evidence.csv** - merged table of target-disease pairs with clinical status from OpenTargets, genetic evidence from OpenTargets and single-cell evidence from DE analysis, for all tested diseases. 

