import os
import scanpy as sc
import pandas as pd

output_dir = '/nfs/team205/ed6/bin/sc_target_evidence/data/'

cxg_metadata = pd.read_csv(output_dir + 'cellxgene_hsapiens_donor_metadata.disease_relevant_annotation.csv', index_col=0)
disease_info_df = cxg_metadata[['disease_relevant_tissue', 'disease_ontology_id', 'disease']].drop_duplicates()
disease_name_mapper = dict(zip(disease_info_df['disease_ontology_id'], disease_info_df['disease']))
disease_tissue_mapper = dict(zip(disease_info_df['disease_ontology_id'], disease_info_df['disease_relevant_tissue']))

h5ad_files = [x for x in os.listdir(output_dir) if x.endswith('all_genes.h5ad')]

all_ct_obs_df = pd.DataFrame()
for f in h5ad_files:
    obs_df = sc.read_h5ad(output_dir + f, backed=True).obs[
        ['high_level_cell_type', 'high_level_cell_type_ontology_term_id', 'disease_ontology_id']
        ]
    obs_df = obs_df.drop_duplicates().reset_index()[['high_level_cell_type', 'high_level_cell_type_ontology_term_id', 'disease_ontology_id']]
    all_ct_obs_df = pd.concat([all_ct_obs_df, obs_df])

all_ct_obs_df = all_ct_obs_df[['high_level_cell_type', 'high_level_cell_type_ontology_term_id']].drop_duplicates()
all_ct_obs_df.to_csv(f'{output_dir}/all_CellOntologies.csv')