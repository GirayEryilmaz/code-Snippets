#Plot a umap from a modality
mu.pl.embedding(pbmc, basis='adt:X_umap')


# Dot plot
def mu_get_features_df(mudata, keys, metadata_col, prefer_raw = True):
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

groupby = 'rna:T_subset_bbknn_pool_leiden_100'
features = ['Gene1', 'Gene2', 'Protein1']
df = mu_get_features_df(Monocyte, keys = features, metadata_col = groupby)
df = df.groupby([groupby]).apply(lambda df: mean_and_fract_cells(df)).reset_index()

fig, ax = plt.subplots(figsize = (10, 3))
df.reset_index().plot.scatter(x=groupby, y='level_1', s='fraction_of_cells', c='Mean', cmap='Reds', edgecolors='grey', ax=ax);
ax.tick_params(axis='x', rotation=90)
ax.set_title('Dotplot')
ax.grid(False)
fig.set_tight_layout(True)
