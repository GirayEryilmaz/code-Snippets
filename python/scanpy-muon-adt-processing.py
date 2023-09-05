import muon as mu
import scanpy as sc

adt_T = sc.read('raw_T_proteins.h5ad')
adt_T = adt_T[:, adt_T.varm['is_in_common'] & (~adt_T.varm['is_control'])]
mu.prot.pp.clr(adt_T)
sc.pp.scale(adt_T)
sc.tl.pca(adt_T, svd_solver='arpack')
batch_key='pool'
sc.external.pp.bbknn(adt_T, batch_key=batch_key, n_pcs=20)
sc.tl.umap(adt_T)
sc.tl.leiden(adt_T, key_added=f'bbknn_{batch_key}_leiden_020', resolution=0.2)
sc.tl.leiden(adt_T, key_added=f'bbknn_{batch_key}_leiden_050', resolution=0.5)
sc.tl.leiden(adt_T, key_added=f'bbknn_{batch_key}_leiden_080', resolution=0.8)
sc.tl.leiden(adt_T, key_added=f'bbknn_{batch_key}_leiden_100', resolution=1.0)

sc.pl.umap(adt_T, color = [f'bbknn_pool_leiden_{num}' for num in ['020', '050', '080', '100']])

adt_T.write('T_proteins.h5ad')
