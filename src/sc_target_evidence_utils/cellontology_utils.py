import os
from typing import List,Union

import networkx as nx
import numpy as np
import obonet
import requests
from networkx import MultiDiGraph


# --- Utils for cell ontology --- #
def get_cellontology_graph(dirpath):
    """
    Load cell ontology graph from file (or download it).

    Parameters
    ----------
    dirpath : str
        Path to directory where ontology file is stored

    Returns
    -------
    graph : networkx.MultiDiGraph
        Cell ontology graph
    """
    CELLONTOLOGY_OBO_URL = "http://obofoundry.org/ontology/cl.html/cl.obo"
    obo_file = dirpath + "/cl.obo"

    # Download and load cell ontology graph
    if not os.path.exists(obo_file):
        print("Downloading the ontology file...")
        response = requests.get(CELLONTOLOGY_OBO_URL)
        if response.status_code == 200:
            with open(obo_file, "wb") as file:
                file.write(response.content)
        else:
            print(f"Failed to download the Cell ontology file. Status code: {response.status_code}")
            exit(1)

    # Load the ontology from the OBO file
    graph = obonet.read_obo(obo_file)
    return graph


def ontology2name(term_id: str, graph: MultiDiGraph) -> str:
    """Get name of ontology term from its ID."""
    return graph.nodes.get(term_id)["name"]


def get_edge_type(term1: str, term2: str, graph: MultiDiGraph) -> str:
    """Get edge type between two terms in the ontology graph"""
    assert term1.startswith("CL:"), "term1 must be a cell ontology term"
    assert term2.startswith("CL:"), "term2 must be a cell ontology term"
    edge_type = [k for k in graph.get_edge_data(term1, term2).keys()][0]
    return edge_type


def get_ancestors(term_id: str, graph: MultiDiGraph, levels: int = 3):
    """Get ancestors of a term in the ontology graph"""
    ancestors = set()
    stack = [(term_id, 0)]  # Using a stack to keep track of terms and their levels

    while stack:
        current_term, current_level = stack.pop()
        ancestors.add(current_term)

        if current_level < levels:
            parents = [
                parent
                for parent in graph.successors(current_term)
                if parent.startswith("CL:") and get_edge_type(current_term, parent, graph) == "is_a"
            ]
            stack.extend([(parent, current_level + 1) for parent in parents])

    ## Order ancestors by shortest path to term_id
    dist2term = [nx.shortest_path_length(graph, source=term_id, target=x) for x in list(ancestors)]
    sorted_combined = sorted(zip(list(ancestors), dist2term), key=lambda x: x[1])
    sorted_ancestors = [item[0] for item in sorted_combined]

    return sorted_ancestors


def rename_cts_to_high_level(
    original_terms: List[str], graph: MultiDiGraph, term_blacklist: Union[List[str], None] = None, depth_factor: int = 5
):
    """
    Rename cell types to common higher level annotation when present.

    Works as follows:
    - for each term search for ancestor terms up to n levels up the ontology graph
    - store original and ancestor terms in the same list
    - remove blacklisted terms (e.g. native cell, animal cell)
    - for each term, replace with the highest level ancestor term in the list

    Parameters
    ----------
    original_terms: List[str]
        list containing cell_type_ontology_terms to convert to higher level.
    graph : networkx.MultiDiGraph
        Cell ontology graph
    terms_blacklist: List[str] or None
        List of terms to exclude from the conversion. Default: None, use predefined set of low quality annotations

    Returns
    -------
    ct_rename_dict : dict
        Dictionary mapping original cell type ontology terms to new high level cell type ontology term.
    """

    ## Exclude blacklist terms
    if term_blacklist is None:
        term_blacklist = [
            "CL:0000000",  # cell
            "CL:0000003",  # native cell
            "CL:0000548",  # animal cell
            "CL:0000081",  # blood cell
            "CL:0000034",  # stem cell
            "CL:0000145",  # professional antigen-presenting cell
        ]

    for t in original_terms:
        if t in term_blacklist:
            original_terms.remove(t)
    ## Get number of ancestors for each term to determine depth level
    ancestors_depth = {}
    for t in original_terms:
        ancestors_depth[t] = np.floor(len(get_ancestors(t, graph, levels=1000)) / depth_factor)

    new_terms = original_terms.copy()
    for i, term in enumerate(original_terms):
        # Check if cell_type ontology term is the most high level in dataset
        upstream_CL = get_ancestors(term, graph, levels=ancestors_depth[term])
        parent_terms = np.intersect1d(upstream_CL, original_terms)
        if len(parent_terms) > 0:
            new_terms[i] = parent_terms[0]

    ct_rename_dict = dict(zip(original_terms, new_terms))
    return ct_rename_dict
