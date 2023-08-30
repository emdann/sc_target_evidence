Contents
--------

Preprocessing scripts

- `clean_ot_table.py` - script to extract clinical status assignment from OpenTargets target-disease table
- `make_target_universe.py` - script to save gene target universes used throughout the analysis (see [`data/target_universe_dict.json`](https://github.com/emdann/sc_target_evidence/blob/master/data/target_universe_dict.json))

Cellxgene metadata filtering and visualization:
- `cellxgene_census_stats.ipynb` - notebook to download and exploratory analysis of cellxgene database metadata
- `OT_scDisease_overlap.ipynb` - notebook for cxg metadata filtering and crossing with OpenTargets tables

Association analysis
- `OR_association_analysis.ipynb` - notebook for odds ratio association analysis between single-cell evidence and clinical status

Old notebooks and scratch analyses are saved in `_misc/`