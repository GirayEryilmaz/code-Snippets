import numpy as np
# import muon as mu
import scanpy as sc
from matplotlib import pyplot as plt

from pandas.api.types import is_numeric_dtype
import warnings

def split_umap(adata, split_by, color_by=None, ncols = 4, basis = 'X_umap', remove_unused_color = False, title = '', title_color = 'black', title_fontsize = 18, remove_legends_except_last = False, na_color = '#F3F3F3', use_raw =True, suptitle_kwargs = {}, **kwargs):
    """
        This can work on both adata and anndata objects
    """
    if color_by is None:
        color_by = split_by
    color_col = sc.get.obs_df(adata, color_by, use_raw=use_raw)

    if is_numeric_dtype(color_col) and remove_unused_color:
        warnings.warn("remove_unused_color is not compatible with numeric color_col. Will ignore remove_unused_color")
        remove_unused_color=False
        
    split_col = adata.obs[split_by]
    split_cats = split_col.cat.categories if split_col.dtype.name == 'category' else sorted(split_col.unique())
    
    nrows = int(np.ceil(len(split_cats) / ncols))
    figsize = kwargs.pop('figsize', (ncols*5, nrows*5))
    

    fig, axs = plt.subplots(nrows, ncols, figsize = figsize)
    axs = axs.flatten()
    
    temp_col_name = 'TEMP_COL_NAME_123'
    for ax, cat in zip(axs, split_cats):
        adata.obs[temp_col_name] = color_col
        adata.obs.loc[split_col!=cat, temp_col_name] = np.nan
        if remove_unused_color: 
            adata.obs[temp_col_name] = adata.obs[temp_col_name].cat.remove_unused_categories()
        sc.pl.embedding(adata, color = temp_col_name, na_color = na_color, ax=ax, show=False, title=cat, na_in_legend=False, basis = basis, **kwargs)
    del adata.obs[temp_col_name]
    
    for ax in axs[len(split_cats):]:
        ax.axis('off')

    
    if remove_legends_except_last:
        for ax in axs[:len(split_cats)-1]:
            ax.get_legend().remove()

    fig.suptitle(title, color = title_color, fontsize = title_fontsize, **suptitle_kwargs)
    fig.tight_layout()
    return fig, axs
