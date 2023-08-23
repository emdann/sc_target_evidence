import pytest
import scanpy as sc
from sc_target_evidence_utils.DE_utils import celltype_marker_targets, disease_marker_targets
from sc_target_evidence_utils.sc_evidence_utils import DE2evidence_celltype, DE2evidence_disease

@pytest.fixture
def pbulk_adata():
    pbulk_adata = sc.read_h5ad('./tests/test_data/test_pbulk.h5ad')
    return (pbulk_adata)

def test_DE2evidence_celltype(pbulk_adata):
    ct_de_genes = celltype_marker_targets(pbulk_adata)
    evidence_dict = DE2evidence_celltype(ct_de_genes)
    assert 'ct_marker_evidence' in evidence_dict.keys(), 'Missing ct_marker_evidence key'
    
def test_DE2evidence_disease(pbulk_adata):
    disease_de_genes, bulk_de_genes = disease_marker_targets(pbulk_adata)
    evidence_dict = DE2evidence_disease(disease_de_genes, bulk_de_genes)
    assert 'disease_evidence' in evidence_dict.keys(), 'Missing disease_evidence key'