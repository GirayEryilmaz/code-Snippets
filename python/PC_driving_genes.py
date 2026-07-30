def get_pc_driving_genes(adata, n, m):
    """
    Return the top driving genes for the first m principal components.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with PCA loadings stored in ``adata.varm["PCs"]``.
    n : int
        Number of genes to return for each principal component.
    m : int
        Number of principal components to include.

    Returns
    -------
    dict[str, list[str]]
        Dictionary mapping ``"PC1"``, ``"PC2"``, ..., ``"PCm"`` to the
        corresponding top n genes, ranked by absolute PCA loading.
    """
    loadings = adata.varm["PCs"]
    genes = adata.var_names.to_numpy()

    return {
        f"PC{i + 1}": genes[
            np.argsort(np.abs(loadings[:, i]))[-n:][::-1]
        ].tolist()
        for i in range(m)
    }

def get_pc_driving_gene(adata, n, m):
    """
    Return the first m driving genes for the nth principal component.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object with PCA loadings stored in ``adata.varm["PCs"]``.
    n : int
        Principal component number, using one-based indexing; ``n=1`` selects PC1.
    m : int
        Number of top driving genes to return.

    Returns
    -------
    list[str]
        Gene names ranked by absolute PCA loading for the selected component.
    """
    loadings = adata.varm["PCs"][:, n - 1]
    order = np.argsort(np.abs(loadings))[-m:][::-1]

    return adata.var_names[order].tolist()


for pc_num, genes in get_pc_driving_genes(adata, 2, 100).items():
    sc.pl.umap(adata, frameon=False, color = genes, show=False, use_raw=False)
    plt.suptitle(pc_num)
    plt.show()
    plt.close()
