# SeparateSeuratAnalysis ----
## Standard 'Seurat' analysis and clustering of cells in each Seurat object or merged Seurat objects

SeparateSeuratAnalysis <- function(seurat_filtered, 
                                   n_variable_feats, #n_variable_feats between 2000-3000 recommended
                                   npcs) {
								   
  require(Seurat)
  require(scran)
  require(dplyr)
  
  if (!is.list(seurat_filtered)) {
    
    seurat_filtered <- Seurat::SplitObject(seurat_filtered, split.by = "orig.ident")
    
  } 
  
  select <- dplyr::select
  
  #load("C:/Users/Fotini/Documents/Meduoa_research/scrna_seq/data/cycle.rda")#load cell cycle genes must be downloaded and loaded
  
  for (seurat in 1:length(seurat_filtered)) {
    
    tryCatch({
      
      set.seed(1993)
      
      print("Analysis of sample: ")
      
      print(seurat_filtered[[seurat]]$orig.ident[[1]])
      
      print("Normalizing...")
      
      seurat_filtered[[seurat]] <- Seurat::NormalizeData(seurat_filtered[[seurat]])#, assay = NULL, normalization.method = "LogNormalize", scale.factor = 10000, margin = 1, verbose = F)
      
      #print("Cell Cycle Scoring...")
      
      #seurat_filtered[[seurat]] <- Seurat::CellCycleScoring(seurat_filtered[[seurat]], g2m.features=g2m_genes, s.features=s_genes)
      
      print("Finding Variable Genes...")
      
      seurat_filtered[[seurat]] <- Seurat::FindVariableFeatures(seurat_filtered[[seurat]], selection.method = "vst", nfeatures = n_variable_feats, verbose = F)
      
      print("Scaling Data Matrix...")
      
      seurat_filtered[[seurat]] <- Seurat::ScaleData(seurat_filtered[[seurat]], features = rownames(seurat_filtered[[seurat]]), verbose = F)#added vars.to.regress, vars.to.regress = c("S.Score", "G2M.Score")
      
      print("Performing PCA...")
      
      seurat_filtered[[seurat]] <- Seurat::RunPCA(object = seurat_filtered[[seurat]], features = VariableFeatures(object = seurat_filtered[[seurat]]), do.print = F, npcs = npcs, verbose = F, seed.use = 1993)
      
      print("Executing UMAP Algorithm...")
      
      seurat_filtered[[seurat]] <- Seurat::RunUMAP(seurat_filtered[[seurat]], reduction = "pca", dims = 1:npcs, verbose = F, seed.use = 1993)
      
      print("Executing TSNE Algorithm...")
      
      seurat_filtered[[seurat]] <- Seurat::RunTSNE(seurat_filtered[[seurat]], reduction = "pca", dims = 1:npcs, verbose = F, seed.use = 1993)
      
      print("Finding Neighbors...")
      
      seurat_filtered[[seurat]] <- Seurat::FindNeighbors(seurat_filtered[[seurat]], dims = 1:npcs, verbose = F)
      
      print("Clustering...")
      
      seurat_filtered[[seurat]] <- Seurat::FindClusters(seurat_filtered[[seurat]], verbose = F, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
      
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  return(seurat_filtered)
}

##eg.how to use function

samples_analyzed <- SeparateSeuratAnalysis(samples_filtered, 2000, 30)
