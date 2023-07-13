# See my comment in the issue: https://github.com/scverse/scanpy/issues/1956#issuecomment-1444361551 (Link date : 12 July 2023)
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt

# get data
adata = sc.datasets.pbmc68k_reduced()

# prepare df
metadata1 = 'phase'
metadata2 = 'bulk_labels'
gene = 'CD8A'

def mean_and_fract_cells_by_gene(df, gene):
    fraction_of_cells = (df[gene]>0).sum()/len(df)
    mean = df.mean(numeric_only=True)[0]
    return pd.Series({'fraction_of_cells' : fraction_of_cells*100, 'Mean': mean})

df = sc.get.obs_df(adata, keys=[gene, metadata1, metadata2])

df = df.groupby([metadata1,metadata2]).apply(lambda df: mean_and_fract_cells_by_gene(df, gene=gene)).reset_index()

# plot
fig, ax = plt.subplots(figsize = (5, 5))
df.plot.scatter(x=metadata1, y=metadata2, s='fraction_of_cells', c='Mean', cmap='Reds', edgecolors='grey', ax=ax);
ax.tick_params(axis='x', rotation=90)
ax.set_title('Dotplot')
# can remove the grid with: ax.grid(False)
fig.set_tight_layout(True)
fig.savefig('dotplot.png', facecolor='white')
