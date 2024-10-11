import plotly.express as px
import scanpy as sc

def umap(adata, color, dot_size, hover_data = [], facet = None, ncols=1, height=700, width=900, **kwargs):
    keys = [color, *hover_data]
    
    if facet:
        keys.append(facet)
    
    df = sc.get.obs_df(adata, keys=keys, obsm_keys=[("X_umap", 0), ("X_umap", 1)])
    fig = px.scatter(
        df,
        x='X_umap-0',
        y='X_umap-1',
        color=color,
        labels={'X_umap-0': 'UMAP 1', 'X_umap-1': 'UMAP 2'},
        facet_col=facet,
        facet_col_wrap = ncols,
        hover_data = hover_data,
        height = height,
        width=width,
        **kwargs
    )
    fig.update_traces(marker_size=dot_size)
    return fig
