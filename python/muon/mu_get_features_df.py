def mu_get_features_df(mudata, keys, prefer_raw = True):
    """
      Similar to scanpy.get.obs_df but less sophisticated for now
    """
    dfs = []
    for mod_key, mod in mudata.mod.items():
        var_names = mod.raw.var_names if mod.raw else mod.var_names
        var_names_of_mod = var_names[var_names.isin(keys)]
        use_raw = prefer_raw and bool(mod.raw)
        df = sc.get.obs_df(mod, keys = var_names_of_mod.values.tolist(), use_raw=use_raw)
        dfs.append(df)
    return pd.concat(dfs, axis=1)
