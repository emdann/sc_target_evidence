import os

import numpy as np
import obonet
import requests


### --- Utils for cell ontology --- ###
def get_cellontology_graph(dirpath):
    """
    Load cell ontology graph from file (or download it)

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


def get_ancestors(term_id: str, graph, levels: int = 3):
    """Get ancestors of a term in the ontology graph"""
    ancestors = set()
    stack = [(term_id, 0)]  # Using a stack to keep track of terms and their levels

    while stack:
        current_term, current_level = stack.pop()
        ancestors.add(current_term)

        if current_level < levels:
            parents = [parent for parent in graph.successors(current_term) if parent.startswith("CL:")]
            stack.extend([(parent, current_level + 1) for parent in parents])

    return ancestors


def get_CL_ancestors(term_id, graph):
    ancestors = get_ancestors(term_id, graph)
    return list(ancestors)


def rename_cts_to_high_level(adata, graph):
    original_terms = adata.obs["cell_type_ontology_term_id"].unique().tolist()
    # Exclude VERY high level annotations (e.g. native cell)
    keep_terms = []
    for o in original_terms:
        if len(get_ancestors(o, graph)) > 5:
            keep_terms.append(o)

    new_terms = original_terms.copy()
    for i, term in enumerate(original_terms):
        # Check if cell_type ontology term is the most high level in dataset
        upstream_CL = get_CL_ancestors(term, graph)
        parent_terms = np.intersect1d(upstream_CL, keep_terms)
        if len(parent_terms) > 0:
            new_terms[i] = parent_terms[0]
    ct_rename_dict = dict(zip(original_terms, new_terms))

    # Rename to higher level if it exists
    while any([ct_rename_dict[v] != v for v in ct_rename_dict.values()]):
        for k, v in ct_rename_dict.items():
            if ct_rename_dict[v] != v:
                ct_rename_dict[k] = ct_rename_dict[v]

    adata.obs["high_level_cell_type_ontology_term_id"] = [
        ct_rename_dict[x] for x in adata.obs["cell_type_ontology_term_id"]
    ]
    return ct_rename_dict
