import numpy as np
# import muon as mu
import scanpy as sc
from matplotlib import pyplot as plt

def split_umap(adata, split_by, color_by=None, ncols = 4, basis = 'X_umap', **kwargs):
    """
        This can work on both adata and anndata objects
    """
    if color_by is None:
        color_by = split_by
    color_col = adata.obs[color_by]
    split_col = adata.obs[split_by]
    split_cats = split_col.cat.categories if split_col.dtype.name == 'category' else sorted(split_col.unique())

    # figsize = kwargs.pop('figsize', plt.rcParams.get('figure.figsize'))
    
    nrows = int(np.ceil(len(split_cats) / ncols))
    figsize = kwargs.pop('figsize', (ncols*5, nrows*5))
    

    fig, axs = plt.subplots(nrows, ncols, figsize = figsize)
    axs = axs.flatten()
    
    temp_col_name = 'TEMP_COL_NAME_123'
    for ax, cat in zip(axs, split_cats):
        adata.obs[temp_col_name] = color_col
        adata.obs.loc[split_col!=cat, temp_col_name] = np.nan
        sc.pl.embedding(adata, color = temp_col_name, na_color = '#F3F3F3', ax=ax, show=False, title=cat, na_in_legend=False, basis = basis, **kwargs)
    del adata.obs[temp_col_name]
    
    for ax in axs[len(split_cats):]:
        ax.axis('off')
    fig.tight_layout()
