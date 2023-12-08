import os
import scanpy as sc
import genomic_features as gf
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
parser.add_argument("--n_hvgs",
                    default=5000,
                    type=int,
                    help="Number of highly variable genes to use in DE analysis.")
args = parser.parse_args()

# Parse args
disease_ontology_id = args.mondo_id
data_dir = args.data_dir
n_hvgs = args.n_hvgs

# Load pseudobulk data
pbulk_adata = sc.read_h5ad(data_dir + f'cellxgene_targets_{disease_ontology_id.replace(":", "_")}.pbulk_all_genes.h5ad')
sc.pp.filter_genes(pbulk_adata, min_cells = 3) # exclude genes expressed in < 3 pseudobulk samples

if not 'feature_name' in pbulk_adata.var:
    ensdb = gf.ensembl.annotation(species="Hsapiens", version="108")
    genes = ensdb.genes()
    pbulk_adata.var['feature_name'] = genes.set_index('gene_id').loc[pbulk_adata.var['feature_id']].gene_name.values

# Run DE analysis
ct_res = DE_utils.celltype_marker_targets(pbulk_adata, n_hvgs = n_hvgs)
ct_res.to_csv(f'{data_dir}/DE_celltype_{disease_ontology_id.replace(":","_")}.hvgs.csv')