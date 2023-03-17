pbmc = readRDS('pmbc.rds')

n = 3000
pbmc_subsampled = pbmc[, sample(colnames(pbmc), size = n, replace=F)]

saveRDS(pbmc_subsampled, 'pbmc_subsampled.rds')
