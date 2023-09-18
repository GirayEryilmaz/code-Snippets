from matplotlib import pyplot as plt

def cel_props_plot(adata, x_axis, to_stack, figsize = (5, 5), ax=None, add_total=False):
    group = adata.obs[[x_axis, to_stack]].groupby([x_axis, to_stack])
    df = group.value_counts()
    percentages = (df/df.groupby(x_axis).sum())
    percentages = percentages.unstack()
    print(percentages.head())
    if add_total:
        totals = adata.obs[to_stack].value_counts()
        totals = totals/totals.sum()
        percentages.loc['Total', :] = totals

    if ax:
        percentages.plot.bar(stacked=True, ax=ax).legend(bbox_to_anchor=(1.0, 0.5));
        ax.grid(axis='x')
    else:
        percentages.plot.bar(stacked=True, figsize=figsize).legend(bbox_to_anchor=(1.0, 0.5));
        plt.grid(axis='x')
