import scanpy as sc
import pandas as pd
import anndata

# ... load the data here ... #

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)


adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, batch_key='batch')

sc.tl.pca(adata, svd_solver='arpack')

# Skip if you are not going to use harmony for batch correction
sc.external.pp.harmony_integrate(adata, 'batch')

sc.pp.neighbors(adata, n_pcs=25)

sc.tl.umap(adata)

# Continue with clustering...
sc.tl.leiden(adata, resolution=0.3, key_added='l030')

# ... Annotate Cells here using marker genes ... #

# Save progress!
adata.write_h5ad('adata.h5ad')
