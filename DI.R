library(clusterCrit)

raw <- psm
# 获取 UMAP 坐标和标签
coords <- Embeddings(raw, reduction = "umap")
labels <- raw$seurat_clusters
celltypes <- raw$CellType

# 计算每个细胞类型上的 Dunn Index
dunn_table <- data.frame()

for (ctype in unique(celltypes)) {
  idx <- which(celltypes == ctype)
  
  # 至少需要2个聚类类别才能计算 DI
  if (length(unique(labels[idx])) > 1) {
    di <- intCriteria(as.matrix(coords[idx, ]), as.integer(labels[idx]), "Dunn")$dunn
    neg_log2_di <- -log2(di)
    dunn_table <- rbind(dunn_table, data.frame(Celltype = ctype, DI = di, NegLog2DI = neg_log2_di))
  } else {
    dunn_table <- rbind(dunn_table, data.frame(Celltype = ctype, DI = NA, NegLog2DI = NA))
  }
}

print(dunn_table)
