# DifferentiallyExpressedGenes ---- 
## Standard 'Seurat' FindAllMarkers wrapper function for the identification of 
## differentially expressed genes between clusters using the clustering result of preference (the appropriate Ident)
## This can be applied in a single Seurat object or a list of Seurat objects.
## The output of this function is a list of lists with the results for each Seurat object.

DifferentiallyExpressedGenes <- function(seurat_analyzed, #list of Seurat objects or single Seurat object
                                         clustering = "RNA_snn_res.0.4"){ #the clustering ident to use for the de analysis
  
  
  require(Seurat)
  require(scran)
  
  if (!is.list(seurat_analyzed)) {
    
    print("Converting into a list...")
    
    temp <- list()
    
    temp[[print(as.character(substitute(seurat_analyzed)))]] <- seurat_analyzed
    
    seurat_analyzed <- temp
    
  } 
  
  de_analysis_results <- list()
  
  for (seurat in 1:length(seurat_analyzed)) {
    
    tryCatch({
      
      de_analysis_results[[print(names(seurat_analyzed[seurat]))]][["seurat_object"]] <- seurat_analyzed[[seurat]]
      
      # count cells in each cluster of cells
      
      Idents(object = de_analysis_results[[seurat]]$seurat_object) <- clustering
      
      de_analysis_results[[seurat]][["clustering_cell_counts"]] <- FetchData(de_analysis_results[[seurat]]$seurat_object, vars = c("ident", "orig.ident")) %>%
        dplyr::count(ident, orig.ident) %>%
        tidyr::spread(ident, n)
      
      # identify marker genes in each cluster
      
      DefaultAssay(de_analysis_results[[seurat]]$seurat_object) <- "RNA"
      
      de_analysis_results[[seurat]][["clustering_markers"]] <- FindAllMarkers(de_analysis_results[[seurat]]$seurat_object,
                                                                              assay = assay,
                                                                              only.pos = T,
                                                                              min.pct = 0.2,
                                                                              logfc.threshold = 0.5,
                                                                              verbose = F)
      # top 10 markers of each cluster
      
      de_analysis_results[[seurat]][["clustering_markers_top10"]] <- de_analysis_results[[seurat]][["clustering_markers"]] %>%
        group_by(cluster) %>%
        top_n(n = 10, wt = avg_log2FC)#choosing top highest expressed markers with the highest log2FC in each cluster 
      
      de_analysis_results[[seurat]][["clustering_markers_top10"]] <- unique(de_analysis_results[[seurat]][["clustering_markers_top10"]]$gene)
      
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  return(de_analysis_results)
  
}


##eg.how to use function

samples_degenes <- DifferentiallyExpressedGenes(samples_integrated, "integrated_snn_res.0.2")

