
# Get labels and axis
plotdf = pbmc.obs[['subset_annotations']].copy()
plotdf['Umap_0'] = pbmc.obsm['X_umap'][:, 0].copy()
plotdf['Umap_1'] = pbmc.obsm['X_umap'][:, 1].copy()

plotdf.to_csv('pbmc_for_interactive_plotting.csv')


# plot
import plotly.express as px
import pandas as pd

df = pd.read_csv('pbmc_for_interactive_plotting.csv')
df = df.sample(frac=1).reset_index()
fig = px.scatter(df , x='Umap_0', y='Umap_1', 
                  color='subset_annotations',
                  category_orders = {'color' : sorted(df['subset_annotations'].unique())},
                  hover_name = 'subset_annotations')
fig.update_traces(marker_size=2)
fig.show()
