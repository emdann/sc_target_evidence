import os
import scanpy as sc
import genomic_features as gf
from sc_target_evidence_utils import DE_utils, cellontology_utils
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("tissue_id",
                    type=str,
                    help="Tissue identifier")
parser.add_argument("--data_dir",
                    default='/nfs/team205/ed6/bin/sc_target_evidence/data/',
                    type=str,
                    help="Directory storing pseudo-bulked expression objects.")
parser.add_argument("--mode",
                    default=None,
                    type=int,
                    help="Either 'hvgs' or 'targets'. If none, both are run If 'hvgs', DE analysis is performed on highly variable genes. If 'targets', DE analysis is performed on preclinical target genes for the tissue of interest.")
parser.add_argument("--n_hvgs",
                    default=5000,
                    type=int,
                    help="Number of highly variable genes to use in DE analysis (only if mode=='hvgs').")
parser.add_argument('--target_table_file', 
                    default = 'filtered_nelson_disease_relevant_tissues_06172024.csv', 
                    type=str, 
                    help="Name of target table csv file (stored in data_dir) (used only if mode == 'targets')")

args = parser.parse_args()

# Parse args
tissue_id = args.tissue_id
data_dir = args.data_dir
mode = args.mode
n_hvgs = args.n_hvgs

if mode is None:
    mode = ['hvgs', 'targets']
else:
    mode = [mode]

# Load pseudobulk data
pbulk_adata = sc.read_h5ad(data_dir + f'cellxgene_targets_{tissue_id.replace(" ", "-")}.pbulk_all_genes.h5ad')
assert 'feature_name' in pbulk_adata.var, "Please provide gene names in the feature_name column of the var slot."

# Run DE analysis
if 'hvgs' in mode:
    ct_res_hvgs = DE_utils.celltype_marker_targets(pbulk_adata, n_hvgs = n_hvgs)
    ct_res_hvgs.to_csv(f'{data_dir}/DE_celltype_{tissue_id.replace(" ", "-")}.hvgs.csv')
if 'targets' in mode:
    try:
        pbulk_adata = pbulk_adata[:, pbulk_adata.var['is_tissue_target']].copy()
    except KeyError:
        # Read target table
        target_table = pd.read_csv(f'{data_dir}/{args.target_table_file}', index_col=0)
        try:
            tissue_targets_df = target_table[~target_table[f'dtr_{tissue_id.replace("-", "_")}'].isna()].copy()
        except KeyError:
            if tissue_id == 'eye':
                tissue_targets_df = target_table[~target_table[f'dtr_retina'].isna()].copy()
            else:
                raise ValueError(f'No target data for {tissue_id}')
        adata.var['is_tissue_target'] = adata.var['feature_name'].isin(tissue_targets_df.target.unique())
        pbulk_adata = pbulk_adata[:, pbulk_adata.var['is_tissue_target']].copy()
    
    ct_res_targets = DE_utils.celltype_marker_targets(pbulk_adata, n_hvgs = None)
    ct_res_targets.to_csv(f'{data_dir}/DE_celltype_{tissue_id.replace(" ", "-")}.targets.csv')