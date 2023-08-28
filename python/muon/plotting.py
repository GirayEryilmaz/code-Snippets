#Plot a umap from a modality
mu.pl.embedding(pbmc, basis='adt:X_umap')


# Dot plot
def mu_dotplot(mudata, var_names, groupby, prefer_raw = True, title='Dotplot', rotation=90, grid=False, figsize=(10, 3), swap_axes=False):
    def mu_get_features_df(mudata, keys, metadata_col, prefer_raw):
        dfs = [mudata.obs[[metadata_col]]]
        for mod_key, mod in mudata.mod.items():
            var_names = mod.raw.var_names if mod.raw else mod.var_names
            var_names_of_mod = var_names[var_names.isin(keys)]
            use_raw = prefer_raw and bool(mod.raw)
            df = sc.get.obs_df(mod, keys = var_names_of_mod.values.tolist(), use_raw=use_raw)
            dfs.append(df)
        return pd.concat(dfs, axis=1)
    def mean_and_fract_cells(df):
        fraction_of_cells = (df>0).sum()/len(df)
        mean = df.mean(numeric_only=True)
        return pd.DataFrame({'fraction_of_cells' : fraction_of_cells*100, 'Mean': mean})
    df = mu_get_features_df(mudata, keys = var_names, metadata_col = groupby, prefer_raw = prefer_raw)
    df = df.groupby([groupby]).apply(lambda df: mean_and_fract_cells(df)).reset_index()
    fig, ax = plt.subplots(figsize = figsize)
    x, y = groupby, 'level_1'
    if swap_axes:
        x, y = y, x
    df.reset_index().plot.scatter(x=x, y=y, s='fraction_of_cells', c='Mean', cmap='Reds', edgecolors='grey', ax=ax, fontsize=12);
    ax.tick_params(axis='x', rotation=rotation)
    ax.set_title(title)
    ax.set_ylabel('Features')
    ax.grid(grid)
    fig.set_tight_layout(True)
