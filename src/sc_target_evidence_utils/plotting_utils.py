import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from networkx import MultiDiGraph
from scanpy.metrics import confusion_matrix

from .cellontology_utils import ontology2name


def plot_celltype_rename(adata_obs, dataset_id, graph: MultiDiGraph, savedir=None):
    """
    Plot original vs aggregated cell types (using output to rename_cts_to_high_level)
    
    Parameters
    ----------
    adata_obs
        obs of anndata with original and high level cell types.
    graph
        Cell ontology graph
    savedir
        Directory to save plots.
    """

    # Make confusion table
    confmat = confusion_matrix(
        "cell_type_ontology_term_id", "high_level_cell_type_ontology_term_id", data=adata_obs, normalize=False
    )
    confmat_nonzero = confmat.replace(0, np.nan)

    # Sort by number of grouped annotations
    df = adata_obs[["cell_type_ontology_term_id", "high_level_cell_type_ontology_term_id"]].drop_duplicates()
    df.columns = ["original", "new"]
    df = pd.merge(df, df.value_counts("new").reset_index())
    df = df.sort_values(["count", "new"], ascending=False)
    confmat_nonzero = confmat_nonzero.loc[df["original"], df["new"].unique()]

    # Rename
    confmat_nonzero.columns = [
        f"{ontology2name(x, graph)} ({x})" if x != "low_quality_annotation" else x for x in confmat_nonzero.columns
    ]
    confmat_nonzero.index = [f"{ontology2name(x, graph)} ({x})" for x in confmat_nonzero.index]

    with plt.rc_context({"figure.figsize": (df.shape[0] / 3, df.shape[0] / 3)}):
        # sns.scatterplot(data=df, x='new_name', y='original_name');
        sns.heatmap(confmat_nonzero, annot=True, cmap="Blues", fmt="g")
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)

        plt.xlabel("High level label", fontsize=14)  # Increase the font size of the x-axis label
        plt.ylabel("Original label", fontsize=14)  # Increase the font size of the y-axis label

        plt.xticks(rotation=90)

        if savedir is not None:
            plt.savefig(f'{savedir}/cellxgene_{dataset_id.replace(":","_")}.celltype_harmonization.pdf', bbox_inches='tight')
            plt.savefig(f'{savedir}/cellxgene_{dataset_id.replace(":","_")}.celltype_harmonization.png', bbox_inches='tight')
   
