
# MergeSeurat ----
## Merge all of the seurat objects in a list

MergeSeurat <- function(seurat_list) {
  
  print(seurat_list)
  
  seurat_merged <- merge(x = seurat_list[[1]],
                         y = seurat_list[c(2:length(seurat_list))],
                         add.cell.id = names(seurat_list),
                         project = "Waldenstrom",
                         merge.data = F)#merge.data=TRUE will keep normalized data matrices and raw
  
  seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-") # calculate mitochondrial gene percentage per cell
  
  seurat_merged[["percent.rb"]] <- PercentageFeatureSet(seurat_merged, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^RACK1") # calculate ribosomal gene percentage per cell
  
  return(seurat_merged)
  
}

##eg.how to use function

samples_merged <- MergeSeurat(samples)
