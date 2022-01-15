# PrintPlots ---- 
## Plotting of some basic plots for the analysis

PrintPlots <- function(de_analysis_results,
                       out_path) {
  
  require(Seurat)
  require(scran)
  require(ggplot2)
  require(ggrepel)
  require(ggnewscale)
  require(ggbeeswarm)
  require(ggthemes)
  require(scales)
  require(RColorBrewer)
  require(Polychrome)
  require(viridis)
  require(ggpubr)
  require(reshape)
  require(cowplot)
  require(patchwork)
  require(purrr)
  require(tidyverse)
  require(rJava)
  require(xlsx)
  
  for (seurat in 1:length(de_analysis_results)) {
    
    tryCatch({ 
      
      # create all the plots first
      
      message("Creating UMAP...")
      umap <- DimPlot(de_analysis_results[[seurat]]$seurat_object, reduction = "umap", label = T, label.size = 4)&theme(aspect.ratio=1)
      message("Creating TSNE...")
      tsne <- DimPlot(de_analysis_results[[seurat]]$seurat_object, reduction = "tsne", label = T, label.size = 4)&theme(aspect.ratio=1)
      message("Creating PCA...")
      pca <- DimPlot(de_analysis_results[[seurat]]$seurat_object, reduction = "pca", label = T, label.size = 4)&theme(aspect.ratio=1)
      message("Creating DOTPLOT...")
      dotplot <- DotPlot(de_analysis_results[[seurat]]$seurat_object,
                         assay = "RNA",
                         features = de_analysis_results[[seurat]][["clustering_markers_top10"]],
                         cluster.idents = F,
                         cols = c("white", "red"),
                         scale.by = "radius") + RotatedAxis() + theme(aspect.ratio = 0.15, 
                                                                      legend.text=element_text(size=8),
                                                                      legend.title=element_text(size=8))
      message("Creating HEATMAP...")
      heatmap <- DoHeatmap(de_analysis_results[[seurat]]$seurat_object, 
                           features = de_analysis_results[[seurat]][["clustering_markers_top10"]])+ theme(aspect.ratio = 2, 
                                                                                                          legend.text=element_text(size=8),
                                                                                                          legend.title=element_text(size=8))
      message("Creating FEATUREPLOTS...")
      featureplots <- FeaturePlot(de_analysis_results[[seurat]]$seurat_object, reduction = "umap", features = de_analysis_results[[seurat]][["clustering_markers_top10"]], label = T, label.size = 2, ncol = 10)&theme(aspect.ratio=1)
      
      message("Printing plots....")
      
      # print plots in PDF
      pdf(paste0(out_path, "/", de_analysis_results[[seurat]]$seurat_object$orig.ident[[1]],"_results_plots.pdf"), width = 30, height = 20)
      
      # umap
      print(umap)
      
      # tsne
      print(tsne)
      
      # pca
      print(pca)
      
      # dotplot of top 10 markers in each cluster
      print(dotplot)
      
      # heatmap
      print(heatmap)
      
      # featureplots of top markers
      print(featureplots)
      
      dev.off()
      dev.off()
      dev.off()
      
      write.xlsx(de_analysis_results[[seurat]][["clustering_markers"]], file = paste0(out_path, "/", de_analysis_results[[seurat]]$seurat_object$orig.ident[1], "_clustering_markers.xlsx"))
      
      write.xlsx(de_analysis_results[[seurat]][["clustering_markers_top10"]], file = paste0(out_path, "/", de_analysis_results[[seurat]]$seurat_object$orig.ident[1], "_clustering_markers_top10.xlsx"))
      
      remove(umap, tsne, pca, dotplot, heatmap, featureplots)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
}

##eg.how to use function

PrintPlots(samples_degenes, out_path = "C:/Users/Fotini/Desktop")

