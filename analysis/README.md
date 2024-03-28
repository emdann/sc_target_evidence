Contents
--------

Preprocessing scripts

- `clean_ot_table.py` - script to extract clinical status assignment from OpenTargets target-disease table
- `make_target_universe.py` - script to save gene target universes used throughout the analysis (see [`data/target_universe_dict.json`](https://github.com/emdann/sc_target_evidence/blob/master/data/target_universe_dict.json))

Cellxgene metadata filtering and visualization:
- `cellxgene_census_stats_v2.ipynb` - notebook to download and exploratory analysis of cellxgene database metadata
- `OT_scDisease_overlap.ipynb` - notebook for cxg metadata filtering and crossing with OpenTargets tables

Association analysis
- `OR_association_analysis.ipynb` - notebook for odds ratio association analysis between single-cell evidence and clinical status
- `fine_vs_coarse.ipynb` - comparison of coarse and fine annotations on lung diseases
- `disease_followup.ipynb` - follow-up target analysis on selected diseases

Open Targets data loading, figures and analysis scripts
- OTplatform_drug_disease_evidence_summary_03062024: This file loads Open Targets direct association evidence from platform json format files and matches Open Targets Disease Terms with Cell x Gene mapped disease terms.
- plots_final_03282024_share_version.Rmd: This notebook loads files containing drug information for the studied diseases and their targets. Code generates plots by year of approval, counts of drug mechanism, and number of approved and exploratory indications per drug. This file also includes multiple linear regression analysis predicting number of approved or investigational indications from year of first approval, presence single cell evidence, and presence of genetic evidence. These associations are merged with features generated from single cell data for target-disease pairs.
- drugs_OTplatform.ipynb: This notebook combines drug information tables downloaded from Open Targets for the 30 disease terms included in this analysis. We include in our analysis only drugs with the indication matching the disease term and not indirectly associated indications for specificity.
- merge_sc_disease_evidence_split_OT_03142024.ipynb: This notebook matches gene-disease evidence using cell x gene disease terms back with gene-disease terms in open targets and loads Open Targets data for drugs and drug mechanisms of action. Drugs can have multiple targets and multiple investigational or approved indications. This notebook searches the target-disease associations at the drug level and also creates a file with mechanism of action information aggregated in Open Targets. 


Additional utility code is saved in `_misc/`
