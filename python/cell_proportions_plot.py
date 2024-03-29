from matplotlib import pyplot as plt

def cel_props_plot(adata, x_axis, to_stack, figsize = (5, 5), ax=None, add_total=False):
    group = adata.obs[[x_axis, to_stack]].groupby([x_axis, to_stack])
    df = group.value_counts()
    percentages = (df/df.groupby(x_axis).sum())
    percentages = percentages.unstack()
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
axs[1].set_title('Person')
for ax in axs:
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));
    ax.set_xlabel(None);
fig.subplots_adjust(wspace=0.5)

def longitudinal_cell_prop_plot(adata, cell_type_of_interest, facade = 'Vaccine', group_by = ['person', 'Visit'], annots_col = 'subset_annotations', wspace=0.3, figsize = (25, 3)):
    ncats = len(adata.obs[facade].cat.categories)
    fig, axs = plt.subplots(1, ncats, figsize = figsize)
    for ax, vax in  zip(axs, adata.obs[facade].cat.categories):
        adata_vax = adata[adata.obs[facade] == vax]
        df = adata_vax.obs[group_by + [annots_col]].groupby(group_by).value_counts()
        percentages = df/df.groupby(group_by).sum()
        
        percentages_of_celltype = percentages.xs(cell_type_of_interest, level=2, drop_level=True).dropna(axis = 'rows')
        
        percentages_of_celltype = percentages_of_celltype.unstack()
        
        percentages_of_celltype.T.plot.line(ax = ax).legend(bbox_to_anchor=(1.0, 0.5));
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));
        ax.set_title(f'{vax} - {cell_type_of_interest}')
        ax.grid(axis='x')
        ax.tick_params(axis=u'both', which=u'both',length=0)
    fig.subplots_adjust(wspace=wspace)

def longitudinal_cell_prop_pbmc_plot(adata, pbmc, cell_type_of_interest, facade = 'Vaccine', group_by = ['person', 'Visit'], annots_col = 'subset_annotations', wspace=0.3, figsize = (25, 3)):
    ncats = len(adata.obs[facade].cat.categories)
    fig, axs = plt.subplots(1, ncats, figsize = figsize)
    for ax, vax in  zip(axs, adata.obs[facade].cat.categories):
        adata_vax = adata[adata.obs[facade] == vax]
        pbmc_vax = pbmc[pbmc.obs[facade] == vax]
        df = adata_vax.obs[group_by + [annots_col]].groupby(group_by).value_counts()
        percentages = df/(pbmc_vax.obs[group_by].groupby(group_by).value_counts())
        
        percentages_of_celltype = percentages.xs(cell_type_of_interest, level=2, drop_level=True).dropna(axis = 'rows')
        
        percentages_of_celltype = percentages_of_celltype.unstack()
        
        percentages_of_celltype.T.plot.line(ax = ax).legend(bbox_to_anchor=(1.0, 0.5));
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));
        ax.set_title(f'{vax} - {cell_type_of_interest}')
        ax.grid(axis='x')
        ax.tick_params(axis=u'both', which=u'both',length=0)
    fig.subplots_adjust(wspace=wspace)

def longitudinal_cell_prop_box_plot(adata, pbmc, cell_type_of_interest, facade = 'Vaccine', group_by = ['Visit', 'person'], 
                                    annots_col = 'subset_annotations', include = 'all', wspace=0.3, figsize = (25, 3), rot = 90, share_y = True):
    ncats = len(adata.obs[facade].cat.categories)
    fig, axs = plt.subplots(1, ncats, figsize = figsize)

    perc_dfs = []
    
    for ax, vax in  zip(axs, adata.obs[facade].cat.categories):
        adata_vax = adata[adata.obs[facade] == vax]
        pbmc_vax = pbmc[pbmc.obs[facade] == vax]
        df = adata_vax.obs[group_by + [annots_col]].groupby(group_by).value_counts()
        percentages = df/(pbmc_vax.obs[group_by].groupby(group_by).value_counts())
        percentages_of_celltype = percentages.xs(cell_type_of_interest, level=2, drop_level=True).dropna(axis = 'rows')
        percentages_of_celltype = percentages_of_celltype.unstack()
        percentages_of_celltype = percentages_of_celltype.T
        
        if include == 'all':
            pass
        else:
            to_include = [i for i in include if i in percentages_of_celltype.columns]
            percentages_of_celltype = percentages_of_celltype[to_include]

        perc_dfs.append(percentages_of_celltype)

    for ax, vax, df in zip(axs, adata.obs[facade].cat.categories, perc_dfs):
        ax = sns.boxplot(data=df, color = 'white', ax=ax)
        ax = sns.swarmplot(data=df, color = 'gray', ax=ax)
        ax = sns.lineplot(data=df.T, ax=ax)
        ax.tick_params(axis='x', labelrotation=90)
        ax.set_xlabel(None);
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5));
        ax.set_title(f'{vax} - {cell_type_of_interest}')
        ax.grid(axis='x')
        ax.tick_params(axis=u'both', which=u'both',length=0)

    if share_y:
        y_max = max([df.values.max() for df in perc_dfs])
        y_max = y_max*1.05
        for ax in axs:
            ax.set_ylim(top=y_max, bottom = 0)
    fig.subplots_adjust(wspace=wspace)
