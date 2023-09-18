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


# example:
fig, ax = plt.subplots(1, 1, figsize = (6, 5))
cel_props_plot(adata, x_axis='clusters', to_stack='Batch', ax = ax, add_total=True)
ax.set_xlabel(None);
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));

# Example 2
fig, axs = plt.subplots(1, 2, figsize = (12, 5))
# axs = axs.flatten()
cel_props_plot(adata, x_axis='Cluster', to_stack='Batch', ax = axs[0], add_total=True)
cel_props_plot(adata, x_axis='Cluster', to_stack='Person', ax = axs[1], add_total=True)
axs[0].set_title('Batch')
axs[0].set_title('Person')
for ax in axs:
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));
    ax.set_xlabel(None);
fig.subplots_adjust(wspace=0.5)
