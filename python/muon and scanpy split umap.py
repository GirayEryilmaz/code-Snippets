import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt

def muon_split_umap(mudata, color_by, split_by, ncol = 4, **kwargs):
    """
        This can work on both mudata and anndata objects
    """
    color_col = mudata.obs[color_by]
    split_col = mudata.obs[split_by]
    split_cats = split_col.cat.categories if split_col.dtype.name == 'category' else sorted(split_col.unique())

    figsize = kwargs.pop('figsize', plt.rcParams.get('figure.figsize'))

    nrow = int(np.ceil(len(split_cats) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize = figsize)
    axs = axs.flatten()
    
    temp_col = 'TEMP_COL_NAME_123'
    for ax, cat in zip(axs, split_cats):
        mudata.obs[temp_col] = color_col
        mudata.obs.loc[split_col!=cat, cat] = np.nan
        mu.pl.embedding(mudata, color = cat, na_color = '#F3F3F3', ax=ax, show=False, title=cat, na_in_legend=False, **kwargs)
    del mudata.obs[temp_col]
    
    for ax in axs[len(split_cats):]:
        ax.axis('off')
    fig.tight_layout()
