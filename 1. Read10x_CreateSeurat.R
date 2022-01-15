# Read10x_CreateSeurat ----

## Specify a path of the directory of the Cell Ranger filtered_feature_bc_matrix
## output and it will return a list all of the respective Seurat objects
## Note: Name the filtered_feature_bc_matrix directories accordingly to have nicely named Seurat objects that represent the biology of each sample

Read10x_CreateSeurat <- function(dir_path) { 
  
  set.seed(1993)
  
  require(Seurat)
  
  objects <- as.list(list.dirs(dir_path, full.names = FALSE, recursive=F))
  
  seurat_list <- list()
  
  for (object in objects) {
    
    print(object)
    
    seurat_data <- Read10X(data.dir = paste0(dir_path, object))
    
    seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                     min.features = 0,
                                     min.cells = 0,
                                     project = object)
    
    seurat_list[[object]] <- seurat_obj
    
    remove(seurat_data, seurat_obj)
  }
  
  return(seurat_list)
  
}

##eg.how to use function

data_path <- "C:/Users/Fotini/Documents/Meduoa_research/scrna_seq/data/test/"

samples <- Read10x_CreateSeurat(data_path)
