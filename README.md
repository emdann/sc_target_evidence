# Drug target evidence in single-cell data
Meta-analysis of drug target evidence in single-cell data

## Contents

- `analysis` - Notebooks and scripts for analysis
- `data` - Metadata and output files
- `src` - Main worflow scripts and functions
    - `sc_target_evidence_utils` - Python package with utility functions
- `tests` - Unit tests for utility functions

## Set-up

```bash
# Make conda env (see also sc-target-evidence-env.yml)
conda create --name sc-target-evidence-env
conda activate sc-target-evidence-env

# Install R dependencies (for DE analysis)
conda install conda-forge::r-base==4.0.5 
Rscript --vanilla -e "install.packages(c('BiocManager'), repos='http://cran.us.r-project.org', lib='/nfs/team205/ed6/miniconda3/envs/sc-target-evidence-env/lib/R/library'); library('BiocManager'); BiocManager::install('glmGamPoi', lib='/nfs/team205/ed6/miniconda3/envs/sc-target-evidence-env/lib/R/library')"
Rscript --vanilla -e "install.packages(c('tidyverse'), repos='http://cran.us.r-project.org', lib='/nfs/team205/ed6/miniconda3/envs/sc-target-evidence-env/lib/R/library')"
Rscript --vanilla -e "library('BiocManager'); BiocManager::install('scater', lib='/nfs/team205/ed6/miniconda3/envs/sc-target-evidence-env/lib/R/library')"

# install utils package
pip install .
```

