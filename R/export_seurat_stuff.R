
library(Seurat)
library(data.table)

# You have to take transpose of this if you want to import this to scanpy
# It is also possible to export this as a sparse matrix which makes more sense ususally.
data_to_write_out <- as.data.frame(as.matrix(seuratObject@assays$RNA@data))
fwrite(x = data_to_write_out, file = "normalized_data.csv")

# Metadata
seuratObject[['celltypes']] = seuratObject@active.ident
write.csv(seuratObject[[c('celltypes', 'orig.ident')]], 'metadata.csv')

# Feature names (gene, protein etc)
write.csv(rownames(seuratObject), 'feature_names.csv', row.names = F)

# The umap
write.csv(seuratObject[['umap']]@cell.embeddings, 'umap.csv')
