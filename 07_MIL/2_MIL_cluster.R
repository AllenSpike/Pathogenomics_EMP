library(jsonlite)
library(Seurat)
library(ggplot2)
library(stringr)
library(harmony)
library(dplyr)
library(patchwork)


# 数据读入 -----------------------------------------------------------------
pred <- fromJSON('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/Clam_mb_prediction.json')
att <- fromJSON('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/Clam_mb_attention.json')
feat <- fromJSON('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/Clam_mb_feature.json')

att <- lapply(att, function(x){
  names(x) <- pred[[1]][["slide"]]
  return(x)
})
feat <- lapply(feat, function(x){
  names(x) <- pred[[1]][["slide"]]
  return(x)
})

# # 汇总5个模型
# att_ls <- lapply(1:175, function(i) {
#   dfs_to_merge <- lapply(att, function(dim) dim[[i]])
#   merged_df <- do.call(cbind, dfs_to_merge)
#   merged_df <- as.data.frame(merged_df)
#   merged_df$id <- str_split(rownames(merged_df), '\\.', simplify = T)[,1]
#   merged_df <- aggregate(.~id,mean,data=merged_df)
#   merged_df <- merged_df[,-1]
#   return(merged_df)
# })
# feat_ls <- lapply(1:175, function(i) {
#   dfs_to_merge <- lapply(feat, function(dim) dim[[i]])
#   do.call(cbind, dfs_to_merge)
# })
# names(att_ls) <- names(feat_ls) <- pred[[1]][["slide"]]

# 选择最佳模型
att_ls <- att[[4]]
feat_ls <- feat[[4]]


# 无监督聚类 ------------------------------------------------------------------
seurat.list <- lapply(names(feat_ls), function(x){
  data <- t(feat_ls[[x]])
  CreateSeuratObject(counts = data, project = x, min.cells = 0, min.features = 0,
                     meta.data = att_ls[[x]])
})
names(seurat.list) <- names(feat_ls)

scobj <- merge(x=seurat.list[[1]], y = seurat.list[-1], project = "clam_mb")
scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj), reduction.name = "pca")

scobj <- RunHarmony(scobj, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")
pcs_idx <- scobj@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:12]))

# ElbowPlot(scobj, reduction = "pca", ndims = 50)
scobj <- FindNeighbors(scobj, reduction = "harmony", dims = pcs_idx)
scobj <- FindClusters(scobj, algorithm = 1, resolution = seq(0.1, 1, by = 0.1))
scobj <- RunUMAP(scobj, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model=TRUE)

scobj@meta.data[["seurat_clusters"]] <- scobj@meta.data[["RNA_snn_res.0.2"]]
table(scobj$orig.ident, scobj$seurat_clusters)

DimPlot(scobj, reduction = "umap", group.by = "orig.ident", raster = FALSE) + theme(legend.position = "none")
p1 <- DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE, group.by = "seurat_clusters", raster=FALSE, 
              # cols = c('0' = '#f9a411', '1' = '#6a75c2', '2' = '#f14166')
              cols = c('dodgerblue', 'tomato', 'mediumseagreen', 'darkorange', 'purple', 'gold')
              ) +
  labs(title = "CACNA1C classification patch atlas")

p2 <- ggplot(scobj@meta.data, aes(x = seurat_clusters, y = X0, fill = seurat_clusters)) +
  geom_boxplot() +
  theme_bw() + 
  labs(y = 'Attention_1', x = 'Cluster', fill = 'Cluster') + 
  scale_fill_manual(values = c('dodgerblue', 'tomato', 'mediumseagreen', 'darkorange', 'purple', 'gold')) +
  ggplot(scobj@meta.data, aes(x = seurat_clusters, y = X1, fill = seurat_clusters)) +
  geom_boxplot() +
  theme_bw() + 
  labs(y = 'Attention_2', x = 'Cluster', fill = 'Cluster') +
  scale_fill_manual(values = c('dodgerblue', 'tomato', 'mediumseagreen', 'darkorange', 'purple', 'gold')) +
  ggplot(scobj@meta.data, aes(x = seurat_clusters, y = X2, fill = seurat_clusters)) +
  geom_boxplot() +
  theme_bw() + 
  labs(y = 'Attention_3', x = 'Cluster', fill = 'Cluster') + 
  scale_fill_manual(values = c('dodgerblue', 'tomato', 'mediumseagreen', 'darkorange', 'purple', 'gold')) +
  plot_layout(ncol = 3, guides = 'collect')

p1/p2 + plot_layout(heights = c(3,1))


save(scobj, file = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/mil/Clam_mb_seurat.RData')


# 空转patch所属cluster-----------------------------------------------------------------
#### 构造seruat对象 ####
feat_gdph <- fromJSON('/home/ljc/0_project/0_ESCC/0_slideflow/gdph/Clam_mb_feature.json')
feat_gdph <- lapply(feat_gdph, function(x){
  names(x) <- c('1909933-6', '1930377-13')
  return(x)
})
feat_gdph_ls <- lapply(1:2, function(i) {
  dfs_to_merge <- lapply(feat_gdph, function(dim) dim[[i]])
  do.call(cbind, dfs_to_merge)
})
names(feat_gdph_ls) <- c('1909933-6', '1930377-13')

seurat.list <- lapply(names(feat_gdph_ls), function(x){
  data <- t(feat_gdph_ls[[x]])
  CreateSeuratObject(counts = data, project = x, min.cells = 0, min.features = 0)
})
names(seurat.list) <- names(feat_gdph_ls)

scobj_gdph <- merge(x=seurat.list[[1]], y = seurat.list[-1], project = "clam_mb")
scobj_gdph <- NormalizeData(scobj_gdph)
scobj_gdph <- FindVariableFeatures(scobj_gdph, selection.method = "vst", nfeatures = 2000)
scobj_gdph <- ScaleData(scobj_gdph, features = rownames(scobj_gdph))
scobj_gdph <- RunPCA(scobj_gdph, features = VariableFeatures(object = scobj_gdph), reduction.name = "pca")

scobj_gdph <- RunHarmony(scobj_gdph, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")
pcs_idx <- scobj_gdph@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj_gdph@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:10]))

scobj_gdph <- FindNeighbors(scobj_gdph, reduction = "harmony", dims = pcs_idx)
scobj_gdph <- RunUMAP(scobj_gdph, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model=TRUE)
DimPlot(scobj_gdph)

#### 映射 ####
anchors <- FindTransferAnchors(
  reference = scobj,
  query = scobj_gdph,
  normalization.method = "LogNormalize",
  reference.reduction = "harmony"
)

scobj_gdph <- MapQuery(
  anchorset = anchors,
  query = scobj_gdph,
  reference = scobj,
  reference.reduction = "harmony",
  reduction.model = "umap",
  refdata = list(
    a = 'RNA_snn_res.0.1',
    b = 'RNA_snn_res.0.2',
    c = 'RNA_snn_res.0.3',
    d = 'RNA_snn_res.0.4',
    e = 'RNA_snn_res.0.5',
    f = 'RNA_snn_res.0.6',
    g = 'RNA_snn_res.0.7',
    h = 'RNA_snn_res.0.8',
    i = 'RNA_snn_res.0.9',
    j = 'RNA_snn_res.1'
  )
)

DimPlot(scobj_gdph, reduction = 'umap',group.by = 'predicted.b')
DimPlot(scobj_gdph, reduction = 'ref.umap',group.by = 'predicted.b')

#### 输出 ####
patch <- filter(scobj_gdph@meta.data, orig.ident == '1909933-6')
write.csv(patch, file = '/home/ljc/0_project/0_ESCC/0_slideflow/gdph/Clam_mb_1909933-6.csv', row.names = F)
patch <- filter(scobj_gdph@meta.data, orig.ident == '1930377-13')
write.csv(patch, file = '/home/ljc/0_project/0_ESCC/0_slideflow/gdph/Clam_mb_1930377-13.csv', row.names = F)

save(scobj_gdph, file = '/home/ljc/0_project/0_ESCC/0_slideflow/gdph/Clam_mb_seurat.RData')