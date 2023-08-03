import pandas as pd

from sc_target_evidence_utils.cellontology_utils import (
    get_ancestors,
    get_cellontology_graph,
    rename_cts_to_high_level,
)


def test_get_ancestors():
    """Test that the ancestors are correct."""
    graph = get_cellontology_graph("./data/")
    term_id = "CL:0000985"
    levels = 3

    expected_ancestors = [
        term_id,
        "CL:0000974",  # long-lived plasma cell
        "CL:0000786",  # plasma cell
        "CL:0000946",  # antibody secreting cell
    ]

    ancestors = get_ancestors(term_id, graph, levels)
    assert set(ancestors) == set(expected_ancestors), "Ancestors don't match"
    assert all([term.startswith("CL:") for term in ancestors]), "Ancestors must be cell ontology terms"

    # Test that the order is right
    assert ancestors[0] == term_id
    assert ancestors[-1] == "CL:0000946"


def test_rename_cts_to_high_level():
    """Test that the cell types are renamed correctly."""
    graph = get_cellontology_graph("./data/")
    adata_obs = pd.read_csv("./tests/test_data/test_obs.csv")
    ct_rename_dict = rename_cts_to_high_level(adata_obs["cell_type_ontology_term_id"].unique().tolist(), graph)
    expected_renamed_n = 13
    assert sum([k == v for k, v in ct_rename_dict.items()]) == expected_renamed_n
