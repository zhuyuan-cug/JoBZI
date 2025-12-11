prepare_data <- function(count_matrix) {
  counts <- as.matrix(count_matrix)
  rownames(counts) <- iconv(rownames(counts), to = "ASCII//TRANSLIT", sub = " ")
  rownames(counts) <- gsub("|", "-", rownames(counts), fixed = TRUE)
  rownames(counts) <- gsub(" +", "-", rownames(counts))
  rownames(counts)[rownames(counts) == ""] <- "UnknownGene" 
  rownames(counts) <- make.unique(rownames(counts))
  counts <- counts[!is.na(rownames(counts)), ]
  sce <- CreateSeuratObject(counts = counts)
  sce <- NormalizeData(sce)   
  sce <- ScaleData(sce)       
  sce <- FindVariableFeatures(object = sce, selection.method = "vst", nfeatures = 2000)
  sce <- RunPCA(sce, npcs=20)  
  return(sce)
}
detect_technical_zeros <- function(sce, k=20, quantile_thresh) {
  norm_data <- as.matrix(GetAssayData(sce, layer = "data"))
  counts <- as.matrix(GetAssayData(sce, layer = "counts"))
  pca_emb <- as.matrix(Embeddings(sce, "pca"))
  
  knn <- get.knn(pca_emb, k=k)
  knn_index <- knn$nn.index
  
  technical_zeros <- matrix(FALSE, nrow(counts), ncol(counts))
  rownames(technical_zeros) <- rownames(counts)
  colnames(technical_zeros) <- colnames(counts)
  
  gene_activity <- rowMeans(norm_data > 0)
  
  for(gene_idx in 1:nrow(counts)) {
    zero_cells <- which(counts[gene_idx, ] == 0)
    if(length(zero_cells) == 0) next
    gene_active <- gene_activity[gene_idx] > quantile(gene_activity, quantile_thresh)
    neighbor_expr <- sapply(zero_cells, function(cell) {
      neighbors <- knn_index[cell, ]
      mean(norm_data[gene_idx, neighbors] > 0)
    })
    tech_flags <- gene_active & (neighbor_expr > quantile_thresh)
    technical_zeros[gene_idx, zero_cells] <- tech_flags
  }
  
  return(technical_zeros)
}


get_zero_index <- function(input_matrix, k=20, quantile_thresh=0.2){
  count_matrix <- input_matrix
  count_matrix <- as.matrix(count_matrix)
  sce <- prepare_data(count_matrix)
  tech_zeros <- detect_technical_zeros(sce, k=k, quantile_thresh=quantile_thresh)
  PP1 <- matrix(as.numeric(tech_zeros), nrow = nrow(tech_zeros))
  M <- 1-PP1
  return(M)
}


