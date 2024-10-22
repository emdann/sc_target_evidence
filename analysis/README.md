Contents
--------

Preprocessing scripts

- `clean_ot_table.py` - script to extract clinical status assignment from OpenTargets target-disease table
- `make_target_universe.py` - script to save gene target universes used throughout the analysis (see [`data/target_universe_dict.json`](https://github.com/emdann/sc_target_evidence/blob/master/data/target_universe_dict.json))

Cellxgene metadata filtering and visualization:
- `cellxgene_census_stats.ipynb` - notebook to download and exploratory analysis of cellxgene database metadata 
- `OT_scDisease_overlap.ipynb` - notebook for cxg metadata filtering and crossing with OpenTargets tables

Association analysis
- `disease_specific_analysis_revision.ipynb` - notebook for odds ratio association analysis between single-cell evidence and clinical status on 25 diseases with patient data in CxG
- `celltype_specific_analysis_revision.ipynb` - notebook for odds ratio association analysis between cell type specificity and clinical status on >200 diseases

Old notebooks and scratch analyses are saved in `_misc/`