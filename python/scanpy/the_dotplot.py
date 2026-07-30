dp = sc.pl.dotplot(adata, 
                   groupby='subsets', 
                   var_group_rotation = 0, 
                   cmap = 'Greens',
                   var_names=['CD4', 'LYZ'],
                   figsize=(2, 2),
                   use_raw=False,
                   return_fig=True
)
dp.add_totals()

# Make gene names (x tick labels) italic
ax = dp.get_axes()["mainplot_ax"]
for label in ax.get_xticklabels():
    label.set_fontstyle("italic")

fig = dp.fig
fig.savefig(".svg", dpi=300, format='svg', bbox_inches="tight")
