def cel_props_plot(adata, x_axis, to_stack):
  """
  Expects anndata (or mudata) and two column names from the obs dataframe
  """
    group = adata.obs[[x_axis, to_stack]].groupby([x_axis, to_stack])
    df = group.value_counts()
    percentages = (df/df.groupby(x_axis).sum())
    percentages = percentages.unstack()

    percentages.plot.bar(stacked=True, figsize=(10,10)).legend(bbox_to_anchor=(1.0, 0.5));
