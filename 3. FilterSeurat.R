# FilterSeurat ----
## Filter all cells in a merged Seurat object with universal filtering parameters with this function

FilterSeurat <- function(seurat_merged, 
                         min_genes = 200, 
                         max_genes = 2500, 
                         max_mito_genes = 15, 
                         min_non_zero_cells_per_gene = 10) {
  
  # filter cells
  
  print(dim(seurat_merged))# cell numbers before filtering
  
  seurat_filtered <- subset(x = seurat_merged, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mito_genes)
  
  print(dim(seurat_filtered))# cell numbers after filtering
  
  # filter genes
  
  counts <- GetAssayData(object = seurat_filtered, slot = "counts")# Output a logical vector for every gene on whether the more than zero counts per cell
  
  nonzero <- (counts > 0)
  
  keep_genes <- Matrix::rowSums(nonzero) >= min_non_zero_cells_per_gene # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  
  filtered_counts <- counts[keep_genes, ] # Only keeping those genes expressed in more than 10 cells
  
  seurat_filtered <- CreateSeuratObject(filtered_counts, meta.data = seurat_filtered@meta.data)# Reassign to filtered seurat object
  
  print(dim(seurat_filtered))# cell and gene numbers after filtering
  
  return(seurat_filtered)
}

##eg.how to use function

samples_filtered <- FilterSeurat(samples_merged, min_genes =  200, max_genes =  2500, max_mito_genes =  15, min_non_zero_cells_per_gene = 10)
