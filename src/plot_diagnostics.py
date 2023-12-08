### Make diagnostic plots ###
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import scipy
import anndata
import genomic_features as gf
import matplotlib

import obonet
from sc_target_evidence_utils import preprocessing_utils

obo_file = "../data/cl.obo"  # downloaded from http://obofoundry.org/ontology/cl.html

# Load the ontology from the OBO file
graph = obonet.read_obo(obo_file)

plt.rcParams['axes.grid'] = False

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
    
    sns.set_context('poster')
    fig = plt.figure(figsize=[7, fig_height])
        
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
        fig.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.n_cells_boxplot.pdf', bbox_inches='tight')
        fig.show()

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
            layer='logcounts')

    dpl.style(cmap='magma', edge_lw=0)
    dpl.legend(title =  'Mean\nlog-normalized\nexpression')
    _ = dpl.get_axes()['mainplot_ax'].set_xlabel(f'Target genes - {disease_name} ({disease_ontology_id})',fontsize=20)

    if savedir is not None:
        plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.target_expression.pdf', bbox_inches='tight')
        plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.target_expression.png', bbox_inches='tight')
        plt.show()
    else:
        plt.show()
    return()
    
    
def plot_celltype_distribution(pbulk_adata, savedir=None):  
    '''
    Plot clustered heatmap of cell type matches per donor.
    '''
    df = pbulk_adata.obs.copy()
    df['donor_id'] = df['donor_id'].astype('str')
    conf_mat = sc.metrics.confusion_matrix('donor_id', 'high_level_cell_type', df, normalize=False)

    disease_palette = sns.color_palette("Set1", n_colors=len(df['disease'].cat.categories))
    disease_color_mapping = dict(zip(df['disease'].cat.categories, disease_palette))

    fig_height, fig_width = (x / 5 for x in conf_mat.shape)

    sns.set_context('paper')

    # Clustermap for the confusion matrix
    disease_df = df[['donor_id', 'disease']].drop_duplicates().set_index('donor_id')
    col_colors = pd.Series([disease_color_mapping[x] for x in disease_df['disease']])
    col_colors.index = disease_df.index
    sns.clustermap(
        conf_mat.T, 
        figsize = (10,10), 
        col_colors = col_colors,
        robust=True,
        xticklabels=True,
        yticklabels=True,
        linewidths=0
    )
    legend_handles = [plt.Line2D([0], [0], marker='o', color='w', label=cat, markerfacecolor=disease_color_mapping[cat], markersize=10) for cat in df['disease'].cat.categories]
    plt.legend(
        handles=legend_handles, 
        title='Disease', frameon=False, 
        loc='upper right', bbox_to_anchor=(0, 1.0))

    if savedir is not None:
        plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.celltype_distribution.pdf', bbox_inches='tight')
        plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.celltype_distribution.png', bbox_inches='tight')
        plt.show()
    else:
        plt.show()


    
def plot_var_stats(pbulk_adata, savedir=None):
    '''Plot mean expression vs number of samples expressing.'''
    pbulk_adata.layers['expr'] = pbulk_adata.X.copy()
    pbulk_adata.layers['expr'][pbulk_adata.layers['expr'].nonzero()] = 1

    pbulk_adata.var['mean_counts'] = np.array(pbulk_adata.X.mean(0)).flatten()
    pbulk_adata.var['mean_logcounts'] = np.array(pbulk_adata.layers['logcounts'].mean(0)).flatten()
    pbulk_adata.var['n_expressing'] = np.array(pbulk_adata.layers['expr'].sum(0)).flatten()

    plt.hist2d(pbulk_adata.var['n_expressing'], pbulk_adata.var['mean_logcounts'], bins=100, norm=matplotlib.colors.LogNorm());
    plt.xlabel('# pseudo-bulks expressing');
    plt.ylabel('Mean log-normalized counts');
    plt.axhline(pbulk_adata.var['mean_logcounts'].quantile(0.90), color='red')
    
    if savedir is not None:
        plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.var_stats.pdf', bbox_inches='tight')
        plt.savefig(f'{savedir}/cellxgene_targets_{disease_ontology_id.replace(":","_")}.var_stats.png', bbox_inches='tight')
        plt.show()
    else:
        plt.show()
    return(None)
    
        
cxg_metadata = pd.read_csv(data_dir + 'cellxgene_hsapiens_donor_metadata.disease_relevant_annotation.csv', index_col=0)

## Load target data
targets = pd.read_csv(data_dir + 'TargetDiseasePairs_OpenTargets_cellXgeneID_12072023.csv', index_col=0)

## Load pseudobulk_object
pbulk_adata = sc.read_h5ad(data_dir + f'cellxgene_targets_{disease_ontology_id.replace(":", "_")}.pbulk_all_genes.h5ad')

## Exclude low quality cells
pbulk_adata = pbulk_adata[pbulk_adata.obs['high_level_cell_type_ontology_term_id'] != 'low_quality_annotation'].copy()

## Preprocess expression
cpms = scipy.sparse.csr_matrix(pbulk_adata.X.T / pbulk_adata.obs['size_factors'].values.flatten()) * 1000000
pbulk_adata.layers['logcounts'] = np.log1p(cpms).T
pbulk_adata.obs['high_level_cell_type'] = [f'{ontology2name(x)} ({x})' for x in pbulk_adata.obs['high_level_cell_type_ontology_term_id'].tolist()]
pbulk_adata.obs['sample_id'] = ['-'.join(x[1:]) for x in pbulk_adata.obs_names.str.split("-")]

if not 'feature_name' in pbulk_adata.var:
    ensdb = gf.ensembl.annotation(species="Hsapiens", version="108")
    genes = ensdb.genes()
    pbulk_adata.var['feature_name'] = genes.set_index('gene_id').loc[pbulk_adata.var['feature_id']].gene_name.values

disease_relevant_tissue = cxg_metadata[cxg_metadata['disease_ontology_id'] == disease_ontology_id].disease_relevant_tissue.unique()[0]
disease_name = cxg_metadata[cxg_metadata['disease_ontology_id'] == disease_ontology_id].disease.unique()[0]

plot_dir = plot_dir + f'/{disease_ontology_id.replace(":", "_")}_{preprocessing_utils.clean_disease_name(disease_name)}'
if not os.path.exists(plot_dir):
    os.mkdir(plot_dir)

plot_var_stats(pbulk_adata, savedir = plot_dir)
plt.clf()
plot_ncells(pbulk_adata, savedir = plot_dir)
plt.clf()
plot_targets_dotplot(pbulk_adata, targets, disease_ontology_id, disease_name, savedir = plot_dir)
plt.clf()
plot_celltype_distribution(pbulk_adata, savedir = plot_dir)
plt.clf()