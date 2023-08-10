save_plot = function(plot, path, w = 10, h = 10){
    pdf(path, height = h, width = w)
    print(plot)
    dev.off()
}

default_opts = options()
sadjust <- function(w = NULL, h = NULL, res=100, ...) {
    if(!is.null(h)){
        options(repr.plot.height = h)
    }
    
    if(!is.null(w)){
        options(repr.plot.width = w)
    }
    options(repr.plot.res = res)
    
}
sreset = function() {
    options(default_opts)
}


mFeaturePlot = function(title="", ...){
    suppressMessages(FeaturePlot(...)) +
    plot_annotation(title = title) &  theme(plot.title = element_text(hjust = 0.5, size= 32))
}

mDimPlot = function(title ='', ...){
    do_raster = TRUE
    suppressMessages(DimPlot(raster = do_raster, ...)) +
    plot_annotation(title = title) &  theme(plot.title = element_text(hjust = 0.5, size= 32))
}

mVlnPlot = function(title ='', ...){
    suppressMessages(VlnPlot(...)) +
    plot_annotation(title = title) &  theme(plot.title = element_text(hjust = 0.5, size= 16))
}


make_feature_plots = function(pbmc, markers, path, dim_plot_reduction = 'umap', temp_folder_loc = '/tmp'){
    temp_dir_name = paste0(replicate(1, sample(c(LETTERS,letters), 15)), collapse = '')
    temp_folder_path = file.path(temp_folder_loc, temp_dir_name)
    dir.create(temp_folder_path)
    plots_per_row = 5
    p = mDimPlot(object = pbmc, reduction = dim_plot_reduction)
    i = 0
    for(key in names(markers)){
        i = i + 1
        gene_names = markers[[key]]
        num_plots = length (gene_names) + 1
        num_rows = ceiling(num_plots / plots_per_row)
        
        do_raster = TRUE
        pdf(file.path(temp_folder_path, sprintf('%02d.pdf',i)), height = num_rows*5, width = 25)
        print(mFeaturePlot(pbmc, features = c(markers[[key]]), ncol=plots_per_row, title=key, raster=do_raster, order=TRUE) + p & NoAxes())
        dev.off()
    }

    system(paste0('pdfunite ', temp_folder_path, '/*.pdf ', path))
    system(paste0('rm -r ', temp_folder_path))
}


make_vln_plots = function(pbmc, markers, path, temp_folder_loc = '/tmp'){
    temp_dir_name = paste0(replicate(1, sample(c(LETTERS,letters), 15)), collapse = '')
    temp_folder_path = file.path(temp_folder_loc, temp_dir_name)
    dir.create(temp_folder_path)
    num_clus = length(unique(Idents(pbmc)))
    i = 0
    for(key in names(markers)){
        i = i + 1
        gene_names = markers[[key]]
        num_genes = length(gene_names)
        pdf(file.path(temp_folder_path, sprintf('%02d.pdf',i)), height = num_clus*1, width = num_genes*1 + 3)
        print(mVlnPlot(pbmc, features = c(markers[[key]]), stack = TRUE, flip = FALSE, title = key) & NoLegend())
        dev.off()
    }

    system(paste0('pdfunite ', temp_folder_path, '/*.pdf ', path))
    system(paste0('rm -r ', temp_folder_path))
}

tp_by_sc = function(seurat_obj, title, clusters_column){
    temp = seurat_obj@meta.data
    temp = group_by(temp, timepoints, !!sym((clusters_column)))
    temp = summarise(temp, n = n())
    temp = mutate(temp, percentage = (n/sum(n)) * 100)
    temp = as.data.frame(ungroup(temp))

    suppressWarnings(ggplot(temp, aes(x = timepoints, y = percentage, fill = !!sym(clusters_column))) +
        geom_col() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, size=20),
             axis.text.y = element_text(size=15)) & plot_annotation(title = title))

}
