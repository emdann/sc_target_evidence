import pytest
import scanpy as sc
from sc_target_evidence_utils.DE_utils import celltype_marker_targets, disease_marker_targets

@pytest.fixture
def pbulk_adata():
    pbulk_adata = sc.read_h5ad('./tests/test_data/test_pbulk.h5ad')
    return (pbulk_adata)

def test_celltype_marker_targets(pbulk_adata):
    ct_de_genes = celltype_marker_targets(pbulk_adata)
    assert 'gene_name' in ct_de_genes.columns, 'Missing gene_name column'
    assert not 'low_quality_annotation' in ct_de_genes['high_level_cell_type_ontology_term_id'], 'Low quality cell type annotations included'
    assert all(ct_de_genes['pval'] <= ct_de_genes['adj_pval']), 'Incorrect p-value adjustment'
    assert ct_de_genes['gene_id'].nunique() == pbulk_adata.n_vars, 'Incorrect number of genes in DE results'

def test_disease_marker_targets(pbulk_adata):
    disease_de_genes, bulk_de_genes = disease_marker_targets(pbulk_adata)
    assert 'gene_name' in disease_de_genes.columns, 'Missing gene_name column'
    assert not 'low_quality_annotation' in disease_de_genes['high_level_cell_type_ontology_term_id'], 'Low quality cell type annotations included'
    assert disease_de_genes['gene_id'].nunique() == pbulk_adata.n_vars, 'Incorrect number of genes in DE results'
    assert bulk_de_genes['gene_id'].shape[0] == pbulk_adata.n_vars, 'Incorrect number of genes in DE results'

