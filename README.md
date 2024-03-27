# Drug target evidence in single-cell data
Meta-analysis of drug target evidence in single-cell data

## Contents

- `analysis` - Notebooks and scripts for analysis
- `data` - Metadata and output files (see [DATA_INFO](https://github.com/emdann/sc_target_evidence/blob/master/DATA_INFO.md) for details on metadata and output files)
- `src` - Main workflow scripts and functions
    - `sc_target_evidence_utils` - Python package with utility functions
- `tests` - Unit tests for utility functions

## Set-up

```bash
# Make conda env (see also sc-target-evidence-env.yml)
conda create --name sc-target-evidence-env
conda activate sc-target-evidence-env

# Install R dependencies (for DE analysis)
conda install conda-forge::r-base==4.0.5 
ENVPATH=$(conda info --envs | grep sc-target-evidence-env | cut -d " " -f 5) # get path to conda environment
Rscript --vanilla -e "install.packages(c('BiocManager'), repos='http://cran.us.r-project.org', lib='${ENVPATH}/lib/R/library'); library('BiocManager'); BiocManager::install('glmGamPoi', lib='${ENVPATH}/lib/R/library')"
Rscript --vanilla -e "install.packages(c('tidyverse'), repos='http://cran.us.r-project.org', lib='${ENVPATH}/lib/R/library')"
Rscript --vanilla -e "library('BiocManager'); BiocManager::install('scater', lib='${ENVPATH}/lib/R/library')"

# install utils package
pip install .
```

## Data pointers

### Metadata

- [suppl_table_diseases.csv](data/suppl_table_diseases.csv) - table of diseases available in CZ CellxGene database considered for study (supplementary table 1)
- [cellxgene_hsapiens_donor_metadata.disease_relevant_annotation.csv](data/cellxgene_hsapiens_donor_metadata.disease_relevant_annotation.csv) - curated donor-level metadata for all adult studies in CxG collection, including grouping of diseases, grouping of tissues, annotation of disease-relevant tissues (supplementary table 2).
- [target_universe_dict.json](data/target_universe_dict.json) - dictionary of all the tested gene target universes used throughout analysis (built in [`analysis/make_target_universe.py`](https://github.com/emdann/sc_target_evidence/blob/master/analysis/make_target_universe.py))

### scRNA-seq data

- [on server] **cellxgene_targets_{MONDO_ID}.pbulk_all_genes.h5ad** - AnnData object with aggregated expression of all available genes in disease-relevant tissue and pseudo-bulk.

### Diagnostic plots

Plot folders: `sc_target_evidence/data/plots/{disease_id}_{disease_name}`

- **cellxgene_{disease_id}.celltype_harmonization.*** - confusion table of original cell ontology annotations and uniformed ontology annotations. Heatmap color and number in cells denotes the number of cells for each category.
- **cellxgene_targets_{disease_id}.n_cells_boxplot.*** - boxplot of numbers of cells per sample and cell type in healthy and disease tissue
- **cellxgene_targets_{disease_id}.target_expression.*** - heatmap of log-normalized expression of a sample of drug targets for the disease
- **cellxgene_targets_{disease_id}.celltype_distribution.*** - confusion table of assignment of uniformed cell type ontology to each donor, to check differences in cell type distribution across donors/datasets/diseases

### Analysis outputs

- **suppl_table_disease_target_evidence.csv** - merged table of target-disease pairs with clinical status from OpenTargets, genetic evidence from OpenTargets and single-cell evidence from DE analysis, for all tested diseases. 


#### Differential expression analysis outputs

DE analysis results: 
- [on server] **DE_celltype_{disease_id}.hvgs.csv** - results of differential expression analysis between cell types, ran on top 7500 highly variable genes (output of [`sc_target_evidence_utils.DE_utils.celltype_marker_targets`](https://github.com/emdann/sc_target_evidence/blob/9e9658d9443f6f1ca642f008ffc18e847982c476/src/sc_target_evidence_utils/DE_utils.py#L155)). Results for all tested cell types are shown. 
- [on server] **DE_diseasecelltype_{disease_id}.hvgs.csv** - results of differential expression analysis between disease and healthy cells in each cell types, ran on top 7500 highly variable genes (output of [`sc_target_evidence_utils.DE_utils.disease_marker_targets`](https://github.com/emdann/sc_target_evidence/blob/9e9658d9443f6f1ca642f008ffc18e847982c476/src/sc_target_evidence_utils/DE_utils.py#L224)). Results for all tested cell types are shown. 
- [on server] **DE_diseasebulk_{disease_id}.hvgs.csv** - results of differential expression analysis between disease and healthy cells in bulked samples (ignoring cell type), ran on top 7500 highly variable genes (output of [`sc_target_evidence_utils.DE_utils.disease_marker_targets`](https://github.com/emdann/sc_target_evidence/blob/9e9658d9443f6f1ca642f008ffc18e847982c476/src/sc_target_evidence_utils/DE_utils.py#L224)).




