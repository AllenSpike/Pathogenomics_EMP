library(Seurat)
library(readxl)
library(pbapply)
library(harmony)
library(AUCell)
library(data.table)
library(openxlsx)
library(stringr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggVennDiagram)
library(ggpmisc)
library(scRNAtoolVis)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggsankey)
library(scCustomize)
library(monocle)
library(future)
library(gghalves)
library(ggpointdensity)
library(viridis)
library(mclust)
library(GSVA)
library(CytoTRACE2)
library(clusterProfiler)
library(GseaVis)
library(ggSCvis)

options(future.globals.maxSize = 8000 * 1024^2)

#### EMT进展构建 ####
scobj_epi_mes <- subset(scobj_epi, Cell2 == 'Mes')

# 重处理
scobj_epi_mes <- NormalizeData(scobj_epi_mes)
scobj_epi_mes <- FindVariableFeatures(scobj_epi_mes, selection.method = "vst", nfeatures = 2000)
scobj_epi_mes <- ScaleData(scobj_epi_mes,
                           features = VariableFeatures(scobj_epi_mes))
scobj_epi_mes <- RunPCA(scobj_epi_mes, features = VariableFeatures(object = scobj_epi_mes), reduction.name = "pca")
scobj_epi_mes <- RunHarmony(scobj_epi_mes, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")

pcs_idx <- scobj_epi_mes@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj_epi_mes@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:10]))

scobj_epi_mes <- RunUMAP(scobj_epi_mes, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model = T)
scobj_epi_mes <- FindNeighbors(scobj_epi_mes, reduction = "harmony", dims = pcs_idx)
scobj_epi_mes <- FindClusters(scobj_epi_mes, resolution = seq(0.1, 1, 0.1), algorithm = 1)
scobj_epi_mes$seurat_clusters <- as.character(scobj_epi_mes$RNA_snn_res.0.1)

DimPlot(scobj_epi_mes, group.by = 'seurat_clusters', label = T, label.box = T)


##### 上皮分值和间充质分值 ##### 
# Genomic and microenvironmental heterogeneity shaping epithelial-to-mesenchymal trajectories in cancer
epi <- c('KRT14', 'KRT17', 'KRT6A', 'KRT5', 'KRT19', 'KRT8', 'KRT16', 'KRT18', 'KRT6B', 'KRT15', 'KRT6C', 'KRTCAP3', 'SFN', 'EPCAM')
mes <- c('VIM', 'CDH2', 'FOXC2', 'SNAI1', 'SNAI2', 'TWIST1', 'FN1', 'ITGB6', 'MMP2', 'MMP3', 'MMP9', 'SOX10', 'GSC', 'ZEB1', 'ZEB2', 'TWIST2')

scobj_epi_mes <- AddModuleScore(scobj_epi_mes, features = list(epi), name = 'EPI')
scobj_epi_mes <- AddModuleScore(scobj_epi_mes, features = list(mes), name = 'MES')
scobj_epi_mes$EMT1 <- scobj_epi_mes$MES1 - scobj_epi_mes$EPI1 # EMT分值计算为间充质分值减去上皮分值

DimPlot(scobj_epi_mes, group.by = 'seurat_clusters', label = T, label.box = F,
        cols = colorRampPalette(brewer.pal(7, "Set2"))(11)) + ggtitle("Mes\n(n = 5,986)")

p1 <- FeaturePlot(scobj_epi_mes, 'EPI1') + ggtitle("Epithelia score")
p2 <- FeaturePlot(scobj_epi_mes, 'MES1') + ggtitle("Mesenchymal score")
p3 <- FeaturePlot(scobj_epi_mes, 'EMT1') + ggtitle("EMT score")
p1 + p2 + p3 + plot_layout(ncol = 3)

p1 <- ggplot(scobj_epi_mes@meta.data, aes(x = seurat_clusters, y = EPI1, fill = seurat_clusters)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "", y = "Epithelia score", x = "Cell type", fill = 'Cell type') +
  scale_fill_manual(values = colorRampPalette(brewer.pal(7, "Set2"))(11))
p2 <- ggplot(scobj_epi_mes@meta.data, aes(x = seurat_clusters, y = MES1, fill = seurat_clusters)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "", y = "Mesenchymal score", x = "Cell type", fill = 'Cell type') +
  scale_fill_manual(values = colorRampPalette(brewer.pal(7, "Set2"))(11))
p3 <- ggplot(scobj_epi_mes@meta.data, aes(x = seurat_clusters, y = EMT1, fill = seurat_clusters)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "", y = "EMT score", x = "Cell type", fill = 'Cell type') +
  scale_fill_manual(values = colorRampPalette(brewer.pal(7, "Set2"))(11))
p1 + p2 + p3 + plot_layout(ncol = 3, guides = 'collect')


##### 拟时序 #####
# 创建cds对象
cell_ann <- scobj_epi_mes@meta.data
cell_ann$celltype <- as.factor(as.character(scobj_epi_mes$seurat_clusters)) 
cell_ann$orig.ident <- as.factor(cell_ann$orig.ident) # 细胞注释
gene_ann <- data.frame(gene_short_name = rownames(scobj_epi_mes@assays$RNA),
                       row.names = rownames(scobj_epi_mes@assays$RNA)) # 基因注释
ct <- GetAssayData(scobj_epi_mes, slot = 'counts') # 表达矩阵
sccds <- newCellDataSet(cellData = as(ct, 'sparseMatrix'),
                        phenoData = new("AnnotatedDataFrame", data = cell_ann),
                        featureData = new("AnnotatedDataFrame", data = gene_ann),
                        expressionFamily = negbinomial.size())
# 预处理
sccds <- estimateSizeFactors(sccds)
sccds <- estimateDispersions(sccds)
gc()

# 排序基因
scobj_epi_mes <- SetIdent(scobj_epi_mes, value = 'seurat_clusters')
marker_mes <- FindAllMarkers(scobj_epi_mes)
marker_mes_f <- group_by(marker_mes, cluster) %>%
  top_n(50, wt = avg_log2FC)
ordering_genes <- unique(marker_mes_f$gene)

# 降维排序
sccds <- setOrderingFilter(sccds, ordering_genes)
sccds <- reduceDimension(sccds,
                         max_components = 2,
                         verbose = T,
                         norm_method = 'vstExprs',
                         reduction_method = 'DDRTree')
gc()

source('./script/OderCells_2.R')
sccds <- orderCells2(sccds)

# 可视化分化轨迹
table(sccds$State, sccds$orig.ident) # 查看State与样本的关系
table(sccds$State, sccds$seurat_clusters) # 查看State与细胞类型的关系

sccds <- orderCells2(sccds, root_state = 2) # 指定起始节点再次排序

plot_cell_trajectory(sccds, cell_size = 1.5, color_by = "seurat_clusters") +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set2"))(11)) +
  labs(color = 'Mes subclusters')
plot_cell_trajectory(sccds, cell_size = 1.5, color_by = "Pseudotime")

plot_cell_trajectory(sccds, cell_size = 1.5, color_by = "State") + 
  plot_cell_trajectory(sccds, cell_size = 1.5, color_by = "seurat_clusters")　+ 
  plot_cell_trajectory(sccds, cell_size = 1.5, color_by = "Pseudotime") 


##### EMT阶段 #####
# 检查拟时序分布
ggplot(pData(sccds), aes(x = Pseudotime)) +
  geom_density() 
ggplot(pData(sccds), aes(x = Pseudotime, y = State, fill = State)) +
  ggridges::geom_density_ridges_gradient() 

# 划分EMT阶段
mclust <- Mclust(pData(sccds)$Pseudotime, G = 3, modelNames = 'E')
summary(mclust)

sccds$EMT_state <- paste0('EMT_', as.character(mclust[["classification"]]))
table(pData(sccds)$State, pData(sccds)$EMT_state)

ggplot(pData(sccds), aes(x = Pseudotime, y = EMT_state, fill = EMT_state)) +
  ggridges::geom_density_ridges_gradient() 

# EMT阶段之间的cutoff
find_crossing = function(dens1, dens2) {
  
  # 确保使用相同的x值范围
  x <- seq(max(min(dens1$x), min(dens2$x)), 
           min(max(dens1$x), max(dens2$x)), 
           length.out = 1000)
  
  # 插值密度值
  y1 <- approx(dens1$x, dens1$y, x)$y
  y2 <- approx(dens2$x, dens2$y, x)$y
  
  # 寻找差值符号变化的位置
  sign_change <- which(diff(sign(y1 - y2)) != 0)
  
  if(length(sign_change) == 0) return(NA)
  
  # 返回第一个交点
  x[sign_change[1]]
}

cutoff_emt1_2 <- find_crossing(density(filter(pData(sccds), EMT_state == 'EMT_1')$Pseudotime), density(filter(pData(sccds), EMT_state == 'EMT_2')$Pseudotime))
cutoff_emt2_3 <- find_crossing(density(filter(pData(sccds), EMT_state == 'EMT_2')$Pseudotime), density(filter(pData(sccds), EMT_state == 'EMT_3')$Pseudotime))

# 重命名
sccds$EMT_state <- ifelse(pData(sccds)$EMT_state == 'EMT_1', 'EMT_early', 
                          ifelse(pData(sccds)$EMT_state == 'EMT_2', 'EMT_stable', 'EMT_diff'))

ggplot(pData(sccds), aes(x = Pseudotime, fill = EMT_state)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(cutoff_emt1_2, cutoff_emt2_3), 
             linetype = "dashed", color = "red") +
  annotate("text", x = cutoff_emt1_2, y = 1, label = round(cutoff_emt1_2, 2), vjust = -1) +
  annotate("text", x = cutoff_emt2_3, y = 1, label = round(cutoff_emt2_3, 2), vjust = -1) +
  labs(y = 'Probability') + 
  theme_bw() +
  scale_fill_manual(values = c('EMT_early' = "#b6dd49", 'EMT_stable' = "#f9e762", 'EMT_diff' = "#74cc64"))

# 桑基图
df <- data.frame(
  Seurat_cluster = as.character(sccds$seurat_clusters),
  time = sccds$Pseudotime,
  Pseudo_state = sccds$State,
  EMT = sccds$EMT1,
  EMT_state = sccds$EMT_state
)

df2 <- make_long(df, Seurat_cluster, Pseudo_state, EMT_state)
ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = .6,
              node.color = "gray30") +
  geom_sankey_label(size = 3, color = "white", fill = "gray40") +
  scale_fill_viridis_d(drop = FALSE) +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5))

# EMT进展示意图
ggplot(df, mapping = aes(x = time, y = EMT, color = Seurat_cluster)) +
  annotate("rect", xmin = 0, xmax = cutoff_emt1_2, ymin = -Inf, ymax = -3, alpha = 1, fill = "#b6dd49") +
  annotate("rect", xmin = cutoff_emt1_2, xmax = cutoff_emt2_3, ymin = -Inf, ymax = -3, alpha = 1, fill = "#f9e762") +
  annotate("rect", xmin = cutoff_emt2_3, xmax = 11.5, ymin = -Inf, ymax = -3, alpha = 1, fill = "#74cc64") +
  geom_point(size = 1) +
  geom_smooth(data = df, aes(x = time, y = EMT), method = "loess", linetype = "dashed", color = 'black', se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set2"))(11)) +
  theme_bw() +
  labs(x = "Pseudo-time", y = "EMT score", color = 'Cell')

ggplot(pData(sccds), aes(x = Pseudotime, y = celltype, fill = celltype)) +
  annotate("rect", xmin = 0, xmax = cutoff_emt1_2, ymin = -Inf, ymax = Inf, alpha = 0.7, fill = "#b6dd49") +
  annotate("rect", xmin = cutoff_emt1_2, xmax = cutoff_emt2_3, ymin = -Inf, ymax = Inf, alpha = 0.7, fill = "#f9e762") +
  annotate("rect", xmin = cutoff_emt2_3, xmax = 11.5, ymin = -Inf, ymax = Inf, alpha = 0.7, fill = "#74cc64") +
  ggridges::geom_density_ridges_gradient() +
  scale_fill_manual(values = colorRampPalette(brewer.pal(7, "Set2"))(11)) +
  labs(fill = 'Cell', y = 'Cell') +
  theme_bw()

