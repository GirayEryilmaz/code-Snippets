def scanpy_violin_plots(adata, feat_list, groupby, ncols=4, **kwargs):
    nrows = np.ceil(len(feat_list)/ncols).astype(int)
    save = kwargs.pop('save', False)
    figsize = (ncols*5, nrows*5)
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize = figsize)
    for i, (feat, ax) in enumerate(zip(feat_list, axs.flatten())):
        sc.pl.violin(adata, feat, groupby=groupby, show=False, ax=ax, **kwargs)
    plt.tight_layout()
    if save:
        path = 'scanpy_violin.pdf' if type(save) == bool else save
        plt.savefig(path)

# Example run
scanpy_violin_plots(adata,
                    list(chain(*markers_dict.values())),
                    groupby='leiden', 
                    jitter=False, 
                    save='violin.png')
