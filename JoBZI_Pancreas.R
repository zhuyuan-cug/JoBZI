library(Seurat)
library(ggplot2)
library(S4Vectors)
library(BiocParallel)
library(BiocNeighbors)
library(aricode)
library(mclust)
library(FNN)
library(Matrix)
source('c:/hsy/BatchEffect/utils.R')
source('c:/hsy/BatchEffect/PS.R')
source('c:/hsy/BatchEffect/zeroIdentify.R')


pancreas.data <- readRDS(file = "c:/hsy/BatchEffect/Data/Pancreas/pancreas_expression_matrix.rds") # 34363*6321
metadata <- readRDS(file = "c:/hsy/BatchEffect/Data/Pancreas/pancreas_metadata.rds")


pancreas <- CreateSeuratObject(counts = pancreas.data,    # 21215 * 6321
                               meta.data = metadata, min.cells = 30,
                               min.features = 200)
pancreas$percent.mito <- PercentageFeatureSet(object = pancreas, pattern = "^MT-") # 21215 * 6321
pancreas <- subset(pancreas, subset = percent.mito < 5) # 21215 * 5842

# 归一化
pancreas <- NormalizeData(pancreas)

pancreas <- FindVariableFeatures(pancreas, nfeatures = 2000)
markers <- c("REG1A", "SST", "GHRL", "VWF",
             "INS", "GCG", "NKX6-1", "PDX1", "IAPP",
             "PCSK2", "AMY2A",
             "MAFA", "MAFB", "IRX1", "IRX2", "ARX",
             "KRT19", "CPA1", "PDGFRB", "HHEX", "ESR1")
hvg <- scan("c:/hsy/BatchEffect/Data/Pancreas/HVGs_1000.txt", what = "c")


pancreas.list <- SplitObject(object = pancreas, split.by = "tech")
# celseq 1004, celseq2 2285, fluidigmc1 638, smartseq2 2394
# perform log-normalization individually for each batch
celseq2 <- NormalizeData(object = pancreas.list[["celseq2"]])
smartseq2 <- NormalizeData(object = pancreas.list[["smartseq2"]])
celseq <- NormalizeData(object = pancreas.list[["celseq"]])
fluidigmc1 <- NormalizeData(object = pancreas.list[["fluidigmc1"]])
B1 <- as.matrix(smartseq2[["RNA"]]@data)
B2 <- as.matrix(celseq2[["RNA"]]@data)
B3 <- as.matrix(celseq[["RNA"]]@data)
B4 <- as.matrix(fluidigmc1[["RNA"]]@data)

batches <- list(B1, B2, B3, B4)
psm.data <- ps(batches = batches, markers = markers, hvg = hvg, merge.order = 1:4)

psm.data <- t(psm.data)
split_sizes <- c(2394, 2285, 1004, 159)  # 每个批次的列数
split_indices <- split(seq_len(ncol(psm.data)), rep(1:4, split_sizes))

# 按列拆分 matrix_data
new_B1 <- psm.data[, split_indices[[1]], drop=FALSE]
new_B2 <- psm.data[, split_indices[[2]], drop=FALSE]
new_B3 <- psm.data[, split_indices[[3]], drop=FALSE]
new_B4 <- psm.data[, split_indices[[4]], drop=FALSE]

index_B2 <- get_zero_index(B2, k = 60, quantile_thresh = 0.8)
index_B3 <- get_zero_index(B3, k = 60, quantile_thresh = 0.8)
index_B4 <- get_zero_index(B4, k = 60, quantile_thresh = 0.8)

index_B2 <- as.matrix(index_B2)
index_B3 <- as.matrix(index_B3)
index_B4 <- as.matrix(index_B4)

new_B20 <- new_B2 * index_B2
new_B30 <- new_B3 * index_B3
new_B40 <- new_B4 * index_B4

all.equal(new_B20, new_B2)
batch_list <- list(new_B30, new_B20, new_B40, new_B1)

new_data <- do.call(cbind, batch_list)


psm <- pancreas
dim(new_data)
dim(psm[['RNA']]@data)

psm[["RNA"]]@data <- new_data
psm <- ScaleData(psm, verbose = T)
psm <- RunPCA(psm, verbose = T)
psm <- RunUMAP(psm, reduction = "pca", dims = 1:30, check_duplicates = FALSE)

# visualization
DimPlot(psm, reduction = "umap", group.by = "celltype")+
  theme(axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  xlab(NULL)+ylab(NULL)+
  labs(title = "Celltype of Our")

DimPlot(psm, reduction = "umap", group.by = "tech")+
  theme(axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_blank())+
  xlab(NULL)+ylab(NULL)+
  labs(title = "Batch of Our")

psm <- FindNeighbors(psm, dims = 1:30)
psm <- FindClusters(psm, resolution = 0.05)

psm@meta.data$true_labels <- as.numeric(as.factor(psm@meta.data$celltype))
psm@meta.data$seurat_clusters <- as.numeric(as.character(psm@meta.data$seurat_clusters))

true_labels <- psm@meta.data$true_labels
predicted_clusters <- psm@meta.data$seurat_clusters

# Calculate ARI using the adjustedRandIndex function
ari <- adjustedRandIndex(true_labels, predicted_clusters)
print(ari)
nmi_value <- NMI(true_labels, predicted_clusters)
print(nmi_value)


calculate_sbatch <- function(distance_matrix, batch_labels) {
  sbatch_values <- c()
  
  for (i in 1:nrow(distance_matrix)) {
    same_batch <- which(batch_labels == batch_labels[i])
    other_batch <- which(batch_labels != batch_labels[i])
    
    a_i <- mean(distance_matrix[i, same_batch])
    b_i <- mean(distance_matrix[i, other_batch]) 
    
    sbatch_i <- (b_i - a_i) / max(a_i, b_i)
    sbatch_values <- c(sbatch_values, sbatch_i)
  }
  
  return(sbatch_values)
}

calculate_aswbatch <- function(distance_matrix, batch_labels) {
  sbatch_values <- calculate_sbatch(distance_matrix, batch_labels)
  asw_batch <- mean(sbatch_values)
  return(1 - asw_batch)
}

distance_matrix <- as.matrix(dist(psm@reductions$pca@cell.embeddings))
batch_labels <- psm@meta.data$tech 

asw_batch_value <- calculate_aswbatch(distance_matrix, batch_labels)
print(asw_batch_value)



