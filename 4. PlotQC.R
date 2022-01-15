# PlotQC ----
## Create .pdf and .jpeg plots of metrics used for filtering before and after

PlotQC <- function(seurat_merged, 
                   seurat_filtered, 
                   out_path = getwd()) {
  
  # before filtering
  
  pdf(paste0(out_path,"qc_before_filtering.pdf"), width = 20, height = 5, bg = "white")
  print(VlnPlot(seurat_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "orig.ident", ncol = 4)&theme(aspect.ratio = 1))
  dev.off()
  
  jpeg(filename = paste0(data_path,"qc_before_filtering.jpeg"), width = 2480, height = 620, pointsize = 12, quality = 600, bg = "white", res = NA)#A5 = 2480 X 1748, A6 = 1748 X 1240, A4 = 3508 X 2480
  print(VlnPlot(samples_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "orig.ident", ncol = 4)&theme(aspect.ratio = 1))
  dev.off()
  
  # after filtering
  
  pdf(paste0(data_path,"qc_after_filtering.pdf"), width = 20, height = 5, bg = "white")
  print(VlnPlot(seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "orig.ident", ncol = 4)&theme(aspect.ratio = 1))
  dev.off()
  
  jpeg(filename = paste0(data_path,"qc_after_filtering.jpeg"), width = 2480, height = 620, pointsize = 12, quality = 600, bg = "white", res = NA)#A5 = 2480 X 1748, A6 = 1748 X 1240, A4 = 3508 X 2480
  print(VlnPlot(samples_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), group.by = "orig.ident", ncol = 4)&theme(aspect.ratio = 1))
  dev.off()
  
}

##eg.how to use function

PlotQC(seurat_merged = samples_merged, seurat_filtered = samples_filtered, data_path)

