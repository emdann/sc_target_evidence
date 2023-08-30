import pandas as pd

data_dir  = '/nfs/team205/ed6/bin/sc_target_evidence/data/'

## Read OT targets evidence
OT_targets_df = pd.read_csv(data_dir + 'TargetDiseasePairs_OpenTargets_cellXgeneID_12072023.csv', index_col=0)
OT_targets_df = OT_targets_df[['targetId','cxg_matched_id', 'genetic_association', 'known_drug']].copy()
OT_targets_df.columns = ['gene_id', 'disease_ontology_id', 'genetic_association', 'known_drug']
drug_score_class = {
    'druggable': 0,
    'safe': 0.1,
    'effective': 0.2,
    'approved': 0.7
}

for k,v in drug_score_class.items():
    OT_targets_df["is_" + k] = OT_targets_df['known_drug'].apply(lambda x: 1 if x > v else 0)

OT_targets_df.to_csv(data_dir + 'TargetDiseasePairs_OpenTargets_cellXgeneID_12072023.clean.csv')