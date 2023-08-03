# sc_target_evidence
Meta-analysis of drug target evidence in single-cell data

## Set-up

```
# Make conda env
conda create --name sc-target-evidence-env
conda activate sc-target-evidence-env

# Install R dependencies (for DE analysis)
conda install conda-forge::r-base==4.0.5 bioconda::bioconductor-edger==3.32.1 rpy2==3.4.2 conda-forge::r-statmod==1.4.37

# install utils package
pip install .
```

