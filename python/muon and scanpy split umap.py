import numpy as np
import muon as mu
from matplotlib import pyplot as plt

def split_umap(mudata, split_by, color_by=None, ncols = 4, basis = 'X_umap', **kwargs):
    """
        This can work on both mudata and anndata objects
    """
    if color_by is None:
        color_by = split_by
    color_col = mudata.obs[color_by]
    split_col = mudata.obs[split_by]
    split_cats = split_col.cat.categories if split_col.dtype.name == 'category' else sorted(split_col.unique())

    figsize = kwargs.pop('figsize', plt.rcParams.get('figure.figsize'))

    nrow = int(np.ceil(len(split_cats) / ncols))
    fig, axs = plt.subplots(nrow, ncols, figsize = figsize)
    axs = axs.flatten()
    
    temp_col_name = 'TEMP_COL_NAME_123'
    for ax, cat in zip(axs, split_cats):
        mudata.obs[temp_col_name] = color_col
        mudata.obs.loc[split_col!=cat, temp_col_name] = np.nan
        mu.pl.embedding(mudata, color = temp_col_name, na_color = '#F3F3F3', ax=ax, show=False, title=cat, na_in_legend=False, basis = basis, **kwargs)
    del mudata.obs[temp_col_name]
    
    for ax in axs[len(split_cats):]:
        ax.axis('off')
    fig.tight_layout()
