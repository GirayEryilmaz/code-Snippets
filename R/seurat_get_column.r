pbmc[['column']] -> returns the column from the metadata as a DATAFRAME
pbmc[['column', drop=TRUE]] or pbmc$column -> returns the column from the metadata as a VECTOR

# multiple columns:
pbmc[[c('col1', 'col2')]]
