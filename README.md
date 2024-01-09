# Drug target evidence in single-cell data
Meta-analysis of drug target evidence in single-cell data

## Contents

- `analysis` - Notebooks and scripts for analysis
- `data` - Metadata and output files (see [DATA_INFO](https://github.com/emdann/sc_target_evidence/blob/master/DATA_INFO.md) for details on metadata and output files)
- `src` - Main workflow scripts and functions
    - `sc_target_evidence_utils` - Python package with utility functions
- `tests` - Unit tests for utility functions

See Teams folder for [Main results figures](https://sanofi.sharepoint.com/:p:/r/sites/singlecellevidencemeta-analysisproject/Shared%20Documents/General/outline%20documents/Figures.pptx?d=w3e15f2572539439c864a7aa7d1e2ddaf&csf=1&web=1&e=MdlZfy) and [paper outline](https://sanofi.sharepoint.com/:w:/r/sites/singlecellevidencemeta-analysisproject/Shared%20Documents/General/outline%20documents/scTargetID%20-%20outline%20-%20et%20comments.docx?d=w9a958626a4304a79987ed593a231e41b&csf=1&web=1&e=efvuJm)

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



