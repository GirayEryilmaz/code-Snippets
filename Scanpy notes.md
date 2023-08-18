Scanpy notes:

You can not run

```
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

```

two times on the same counts. This will corrupt the data and you will get an error when you run 
`sc.pp.highly_variable_genes(adata)` afterwards.

If you want to find variable genes on a subset that has already been processed with these steps, you can skip th above two and directly do `sc.pp.highly_variable_genes`. However keep in mind that (asuming you have removed non-higly variable genes from the object as a part of the standat pipeline) you will be recalculating the `highly_variable_genes` within the subset of genes you are left with. If you want to find the `highly_variable_genes` for that subset of cells among all the cells, you need to bring back the `raw` as `subset.raw.to_adata()`. Note that the standart pipeline saves the `raw` after normalization and log1p so you won't need to (and in fact wont be able to do the again as that would corrupt your data.)

You can try saving raw count to layers too. However keep in midn that `raw` keeps genes when you remove them from the adata object, layers don't.


A note about reading anndata objects, you will need to do `adata.uns['log1p']["base"] = None` everytime after you read data beacuse `None`s are not properly serialized while saving the data. This is a bug! There is a fix pending. See (this scanpy issue: KeyError: 'base'.)[https://github.com/scverse/scanpy/issues/2239]