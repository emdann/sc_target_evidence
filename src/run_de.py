import scanpy as sc
from sc_target_evidence_utils import DE_utils, cellontology_utils
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("mondo_id",
                    type=str,
                    help="Mondo ID of disease of interest (format should be MONDO:00000000)")
parser.add_argument("--data_dir",
                    default='/nfs/team205/ed6/bin/sc_target_evidence/data/',
                    type=str,
                    help="Directory storing pseudo-bulked expression objects.")
args = parser.parse_args()

# Parse args
disease_ontology_id = args.mondo_id
data_dir = args.data_dir

# Load pseudobulk data
pbulk_adata = sc.read_h5ad(data_dir + f'cellxgene_targets_{disease_ontology_id.replace(":", "_")}.pbulk_all_OT_targets.h5ad')

# clean pseudobulk data
clean_disease(pbulk_adata)
graph = cellontology_utils.get_cellontology_graph(data_dir)
pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'] != 'low_quality_annotation'].copy()
pbulk_adata.obs['high_level_cell_type'] = [f'{cellontology_utils.ontology2name(x, graph)} ({x})' for x in pbulk_adata.obs['high_level_cell_type_ontology_term_id'].tolist()]
pbulk_adata.obs['sample_id'] = ['-'.join(x[1:]) for x in pbulk_adata.obs_names.str.split("-")]
pbulk_adata.obs['disease_ontology_id'] = disease_ontology_id

# Exclude blacklisted assays
assay_blacklist = [
    'BD Rhapsody Targeted mRNA',
    'STRT-seq',
    'inDrop'
    ]

pbulk_adata = pbulk_adata[~pbulk_adata.obs['assay'].isin(assay_blacklist)].copy()

# Run DE analysis
ct_res = DE_utils.celltype_marker_targets(pbulk_adata)
ct_res.to_csv(f'{data_dir}/DE_celltype_{disease_ontology_id.replace(":","_")}.all_targets.csv')

disease_res, bulk_res = DE_utils.disease_marker_targets(pbulk_adata)
disease_res.to_csv(f'{data_dir}/DE_diseasecelltype_{disease_ontology_id.replace(":","_")}.all_targets.csv')
bulk_res.to_csv(f'{data_dir}/DE_diseasebulk_{disease_ontology_id.replace(":","_")}.all_targets.csv')