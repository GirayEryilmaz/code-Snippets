# use the file as: source('/path/plotting_with_seurat.R')
the_colors = c('NK_C1' = '#fee000',
         'NK_XCL1' = '#f2e4a0',
         'NK_Cycling' = '#ccb72d',
         
         'CD8_Naive' = '#f37421',
         'CD8_TEMRA' = '#d28529',
         'CD8_GZMK' = '#fba919',
         'CD8_MAIT' = '#fbb36a',
         'Naive_Cycling' = '#ffd592',
         'Cytotoxic_Cycling' = '#7f7f7f',
         'CD8_GD' = '#80622f',
         
         'CD4_ISG_hi' = '#697d35',
         'CD4_CTL' = '#90aa3c',
         'CD4_Naive' = '#193a1c',
         'CD4_Memory_0' = '#1c572b',
         'CD4_Memory_1' = '#1c7b3d',
         'CD4_Memory_2_JUN' = '#3cb54a',
         'CD4_Memory_4' = '#74c168',
         'CD4_Memory_5' = '#b3daab',
         
         'Tr_B' = '#307ec1',
         'Naive_B' = '#41b8ea',
         'MBC' = '#283779',
         'MBC_ABC' = '#4c459c',
         'MBC_IgM' = '#96daf7',
         'ISG_hi_B' = '#005b96',
         
         'CD4_TREG_1_S100A4' = '#137d82',
         'CD4_IKZF1' = '#56bbbf',
         'CD4_TREG_0' = '#a8dde3',
         'CD4_TREG_4_JUN' = '#0b5657',
         'CD4_TREG_2_CTLA4' = '#c19952',
         'CD4_TEMRA' = '#50C878',
         
         'HSPC' = '#b0479a',
         
         'pDC' = '#a5a4a4',
         'Plasma' = '#232323',
         'MgK_pDC' = '#932169',
         
         'cDC2' = '#d84598',
         'cDC1' = '#771215',
         'AXL_DC' = '#a41e21',
         'Mono_cDC' = '#ed2024',
         
         'CD14_ISG++' = '#f15d64',
         'CD14' = '#f6a2a7',
         'CD16' = '#f9d3d7'
)

colors_25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

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


make_feature_plots = function(pbmc, markers, path, reduction = 'umap', temp_folder_loc = '/tmp'){
    temp_dir_name = paste0(replicate(1, sample(c(LETTERS,letters), 15)), collapse = '')
    temp_folder_path = file.path(temp_folder_loc, temp_dir_name)
    dir.create(temp_folder_path)
    plots_per_row = 5
    p = mDimPlot(object = pbmc, reduction = reduction)
    i = 0
    for(key in names(markers)){
        i = i + 1
        gene_names = markers[[key]]
        num_plots = length (gene_names) + 1
        num_rows = ceiling(num_plots / plots_per_row)
        
        do_raster = TRUE
        pdf(file.path(temp_folder_path, sprintf('%02d.pdf',i)), height = num_rows*5, width = 25)
        print(mFeaturePlot(pbmc, features = c(markers[[key]]), ncol=plots_per_row, title=key, raster=do_raster, order=TRUE, reduction=reduction) + p & NoAxes())
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
