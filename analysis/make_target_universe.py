### Prepare targets universes ### 
import os
import json
import numpy as np
import pandas as pd

## Get table of all genes
import genomic_features as gf
ensdb = gf.ensembl.annotation(species="Hsapiens", version="108")
genes = ensdb.genes()

data_dir = '/nfs/team205/ed6/bin/sc_target_evidence/data/'
universe_dict = {}

## Protein-coding genes
protein_coding_genes = genes[
    (genes['gene_biotype'] == 'protein_coding') & \
    (genes['description'] != 'novel protein') & \
    (genes['seq_name'].isin([str(x) for x in np.arange(1,23)] + ['X', 'Y']) )
]
universe_dict['protein_coding_targets'] = protein_coding_genes.gene_id.tolist()

## Small molecule / antibody tractable targets
nelson_anno_dataset = data_dir + 'genetic_support/data/gene_lists/' ## Downloaded from https://github.com/ericminikel/genetic_support
gene_annotations = os.listdir(nelson_anno_dataset)

annos = {}
for g in gene_annotations:
    anno_name = g.split(".tsv")[0]
    annos[anno_name] = pd.read_table(nelson_anno_dataset + g, header=None).values.flatten()

ens_ids_to_name = genes[['gene_id', 'gene_name']].copy()

for g, g_ls in annos.items():
    ens_ids_to_name[g] = ens_ids_to_name['gene_name'].isin(g_ls).astype('int')
    
universe_dict['sm_tractable_targets'] = ens_ids_to_name[ens_ids_to_name['sm_tractable'] == 1 ].gene_id.tolist()
universe_dict['ab_tractable_targets'] = ens_ids_to_name[ens_ids_to_name['ab_tractable'] == 1 ].gene_id.tolist()

json_file = data_dir + '/target_universe_dict.json'
with open(json_file, "w") as json_file:
    json.dump(universe_dict, json_file)