### Make diagnostic plots ###
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import scipy
import anndata

import obonet

obo_file = "../data/cl.obo"  # downloaded from http://obofoundry.org/ontology/cl.html

# Load the ontology from the OBO file
graph = obonet.read_obo(obo_file)

def ontology2name(o):
    return(graph.nodes.get(o)['name'])


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("mondo_id",
                    type=str,
                    help="Mondo ID of disease of interest (format should be MONDO:00000000)")
parser.add_argument("--data_dir",
                    default='/nfs/team205/ed6/bin/sc_target_evidence/data/',
                    type=str,
                    help="Directory of pseudo-bulked expression objects.")
parser.add_argument("--plot_dir",
                    default='/home/jovyan/mount/gdrive/sc_targetID/plots/',
                    type=str,
                    help="Directory to save plots.")
args = parser.parse_args()

# Parse args
disease_ontology_id = args.mondo_id
data_dir = args.data_dir
plot_dir = args.plot_dir

def plot_ncells(pbulk_adata, savedir=None):
    '''Plot number of cells x donor and cell type in pseudobulk dataset'''
    pl_df = pbulk_adata.obs[['high_level_cell_type', 'sample_id', 'disease', 'n_cells']].copy()
    pl_df.loc[:,'log10_n_cells'] = np.log10(pl_df['n_cells'])

    sorted_cts = pl_df.value_counts('high_level_cell_type').index.tolist()
    pl_df['high_level_cell_type'] 

    fig_height = len(sorted_cts)

    with plt.rc_context({'figure.figsize':[7, fig_height]}):
        sns.set_context('poster')
        sns.boxplot(data=pl_df, y='high_level_cell_type', x='log10_n_cells', 
                    hue='disease', dodge=True, palette='Set1', showfliers=False)
        sns.stripplot(data=pl_df, y='high_level_cell_type', x='log10_n_cells', 
                      hue='disease', dodge=True, color='black', s=3);
        plt.legend(
          bbox_to_anchor=(1.05, 1), # relative position on x and y axis (> 1 indicates outside of axis)
          loc='upper left', # equiv to hjust/vjust in ggplot
          borderaxespad=0, # The pad between the axes and legend border, in font-size units.
          frameon=False,
          title='Disease')
        plt.xlabel("log10(# cells)");
        plt.ylabel("High level cell type");
        plt.title(f'{disease_name} ({disease_ontology_id}) - {disease_relevant_tissue}');
        if savedir is not None:
            plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.n_cells_boxplot.pdf', bbox_inches='tight')

def plot_targets_dotplot(
    pbulk_adata: anndata.AnnData,
    targets: pd.DataFrame,
    disease_ontology_id: str,
    disease_name: str,
    include_other_targets: bool = True,
    split_drug_targets_by_phase: bool = True,
    max_genes: int = 20,
    savedir=None
    ):
    '''
    Plot expression of OT target genes in pseudobulk datasets.
    '''

    ## Genes to plot ## 
    if "_" in disease_ontology_id:
        disease_ontology_id = disease_ontology_id.replace("_", ":")

    # Drug targets
    disease_targets =targets[(targets['targetId'].isin(pbulk_adata.var_names))].copy()
    disease_targets['target_name'] = pbulk_adata.var.loc[disease_targets.targetId.tolist()]['feature_name'].tolist()
    disease_drug_targets = disease_targets[disease_targets['known_drug'] > 0].copy()
    drug_target_genes = disease_targets.target_name.tolist()
    if split_drug_targets_by_phase:
        disease_drug_targets['phase_class'] = pd.cut(disease_drug_targets.known_drug, [0, 0.1, 0.2, 0.7, 1.1], include_lowest=False, labels = ['druggable','safe', 'effective', 'approved'])
        pl_target_dict = disease_drug_targets.groupby('phase_class')['target_name'].apply(list).to_dict()
    else:
        pl_target_dict = {}
        pl_target_dict['Drug targets'] = drug_target_genes

    # Other targets
    if include_other_targets:
        other_target_genes = disease_targets[disease_targets['known_drug'] == 0].target_name.tolist()
        other_target_genes = np.random.choice(other_target_genes, size=len(drug_target_genes))
        pl_target_dict['Other targets'] = other_target_genes

    if max_genes is not None:
        for cl,target_ls in pl_target_dict.items():
            if len(target_ls) > max_genes:
                pl_target_dict[cl] = np.random.choice(target_ls, max_genes)
    
    ## Plot target expression ##
    sc.set_figure_params(scanpy=True, fontsize=20)
    dpl = sc.pl.MatrixPlot(
            pbulk_adata, 
            pl_target_dict, 'high_level_cell_type', 
            gene_symbols='feature_name', 
            layer='logcounts',
        
            )

    dpl.style(cmap='magma', edge_lw=0)
    dpl.legend(title =  'Mean\nlog-normalized\nexpression')
    _ = dpl.get_axes()['mainplot_ax'].set_xlabel(f'Target genes - {disease_name} ({disease_ontology_id})',fontsize=20)

    if savedir is not None:
        plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.target_expression.pdf', bbox_inches='tight')
        plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.target_expression.png', bbox_inches='tight')
    
    plt.show()

        
cxg_metadata = pd.read_csv(data_dir + 'cellxgene_hsapiens_donor_metadata.disease_relevant_annotation.csv', index_col=0)

## Load target data
targets = pd.read_csv(data_dir + 'TargetDiseasePairs_OpenTargets_cellXgeneID_12072023.csv', index_col=0)

## Load pseudobulk_object
pbulk_adata = sc.read_h5ad(data_dir + f'cellxgene_targets_{disease_ontology_id.replace(":", "_")}.pbulk_all_OT_targets.h5ad')

## Exclude low quality cells
pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'] != 'low_quality_annotation'].copy()

## Preprocess expression
cpms = scipy.sparse.csr_matrix(pbulk_adata.X.T / pbulk_adata.obs['size_factors'].values.flatten()) * 1000000
pbulk_adata.layers['logcounts'] = np.log1p(cpms).T
pbulk_adata.obs['high_level_cell_type'] = [f'{ontology2name(x)} ({x})' for x in pbulk_adata.obs['high_level_cell_type_ontology_term_id'].tolist()]
pbulk_adata.obs['sample_id'] = ['-'.join(x[1:]) for x in pbulk_adata.obs_names.str.split("-")]

disease_relevant_tissue = cxg_metadata[cxg_metadata['disease_ontology_id'] == disease_ontology_id].disease_relevant_tissue.unique()[0]
disease_name = cxg_metadata[cxg_metadata['disease_ontology_id'] == disease_ontology_id].disease.unique()[0]

plot_ncells(pbulk_adata, savedir = plot_dir)
plot_targets_dotplot(pbulk_adata, targets, disease_ontology_id, disease_name, savedir = plot_dir)