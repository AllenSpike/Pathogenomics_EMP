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


#### 单细胞数据处理 ####
# 元信息
ref_meta <- fread('/data/chenrz/scseq_gse160269/GSE160269_series_matrix.txt.gz', fill=TRUE)
ref_meta <- t(ref_meta[62:101,])
colnames(ref_meta) <- ref_meta[1,]
ref_meta <- as.data.frame(ref_meta[-1,])

# 基因注释
sheets <- excel_sheets('/data/chenrz/scseq_gse160269/41467_2021_25539_MOESM9_ESM.xlsx')
ref_anno <- lapply(sheets, function(x) read.xlsx('/data/chenrz/scseq_gse160269/41467_2021_25539_MOESM9_ESM.xlsx', sheet = x, startRow = 3))
names(ref_anno) <- sheets
ref_anno[['Epithelia']] <- read.xlsx('/data/chenrz/scseq_gse160269/41467_2021_25539_MOESM8_ESM.xlsx', sheet = 1, startRow = 3)

# 表达谱
cells <- c('Bcell', 'Tcell', 'Endothelial', 'Epithelia', 'Fibroblast', 'FRC', 'Myeloid', 'Pericytes')
exp_ls <- pblapply(cells, function(x){
  ref <- fread(paste0('/data/chenrz/scseq_gse160269/GSE160269_UMI_matrix_', x, '.txt.gz'))
  exp <- as.matrix(ref[,-1])
  rownames(exp) <- ref$V1
  return(exp)
})
names(exp_ls) <- cells

# 创建seurat对象并重新处理
scobj.ls <- pblapply(cells, function(x){
  CreateSeuratObject(counts = as(exp_ls[[x]], 'dgCMatrix'), project = x, min.cells = 3, min.features = 200)
})
names(scobj.ls) <- cells

scobj <- merge(x=scobj.ls[[1]], y=scobj.ls[-1], project = "10x")
rm(scobj.ls)
gc()

scobj$Cell <- scobj$orig.ident
scobj$orig.ident <- str_extract_all(colnames(scobj), '^P\\d+[A-Z]-\\D', simplify = T)[,1]
scobj$percent.mt <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj$percent.ribo <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")

scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj), reduction.name = "pca")
scobj <- RunHarmony(scobj, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:10, reduction.name = "umap", return.model = T)
DimPlot(scobj, group.by = 'Cell', raster = F, label = T, label.box = T)

## 
scobj[["RNA"]] <- as(object = scobj[["RNA"]], Class = "Assay")
save(scobj, file = '/data/chenrz/scseq_gse160269/scobj.RData')

mycols = c("Fibroblast" = "#b64ba8",
           "Myeloid" = "#dcb08d",
           "Mast cell" = "#bdbdbd",
           "Endothelial" = "#97cb2b",
           "FRC" = "#f6cacc",
           "Pericytes" = "#d33b5f",
           "Bcell" = "#0065a1",
           "Tcell" = "#ff9f00",
           "Epithelia" = "#e42f44")

DimPlot(scobj, group.by = 'Cell', raster = T, label = T, cols = mycols, pt.size = 1) + 
  labs(title = "ESCC atlas of Zhang_2021\n(n = 208,659)")

## 
epi <- c('EPCAM', 'SFN', 'KRT5', 'KRT14') # 上皮细胞
fibro <- c('FN1', 'DCN', 'COL1A1', 'COL1A2', 'COL3A1', 'COL6A1') # 成纤维细胞
frc <- c('CCL21', 'IL7', 'PDPN') # 成纤维网状细胞
endo <- c('VWF', 'PECAM1', 'ENG', 'CDH5') # 内皮细胞
peri <- c('NOTCH3', 'MCAM', 'CD146', 'RGS5') # 周细胞
tcell <- c('CD2', 'CD3D', 'CD3E', 'CD3G') # T细胞
bcell <- c('CD19', 'CD79A', 'MS4A1', 'JCHAIN', 'MZB1') # B细胞
Mye <- c('LYZ','CD68','TYROBP') # 髓系细胞

markers <- unique(c(bcell, endo, epi, fibro, frc, Mye, peri, tcell))

DotPlot(scobj, assay = 'RNA', 
        features = markers, 
        group.by = "Cell") + 
  RotatedAxis() +
  scale_x_discrete("") +
  scale_y_discrete("") +
  theme(
    panel.border = element_rect(color="black"),
    panel.spacing = unit(1, "mm"),
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = 'outlet', 
    axis.line = element_blank(),
    axis.text.x = element_text(size=10),
    text = element_text(size=12),
  ) +
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))


#### 单细胞亚群注释 ####
# Dissecting esophageal squamous-cell carcinoma ecosystem by single-cell transcriptomic analysis

##### FRC #####
scobj_frc <- subset(scobj, Cell == 'FRC')
scobj_frc$Cell2 <- 'FRC'

##### Epi #####
## 重处理
scobj_epi <- subset(scobj, Cell == 'Epithelia')
scobj_epi <- NormalizeData(scobj_epi)
scobj_epi <- FindVariableFeatures(scobj_epi, selection.method = "vst", nfeatures = 2000)
scobj_epi <- ScaleData(scobj_epi,
                       features = VariableFeatures(scobj_epi))
scobj_epi <- RunPCA(scobj_epi, features = VariableFeatures(object = scobj_epi), reduction.name = "pca")
scobj_epi <- RunHarmony(scobj_epi, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")

pcs_idx <- scobj_epi@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj_epi@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:10]))

scobj_epi <- RunUMAP(scobj_epi, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model = T)
scobj_epi <- FindNeighbors(scobj_epi, reduction = "harmony", dims = pcs_idx)
scobj_epi <- FindClusters(scobj_epi, resolution = 1, algorithm = 1)
scobj_epi$seurat_clusters <- as.character(scobj_epi$RNA_snn_res.1)
DimPlot(scobj_epi, group.by = 'seurat_clusters', label = T, label.box = T)

gc()

# 标记基因分值
marker <- as.list(ref_anno$Epithelia)
exprMatrix <- as(GetAssayData(scobj_epi, assay = "RNA", slot = "scale.data"), "dgCMatrix")
cells_rankings <- AUCell_buildRankings(exprMatrix)
cells_AUC <- AUCell_calcAUC(marker, cells_rankings)
auc_matrix <- getAUC(cells_AUC)

scobj_epi$Cell2 <- apply(auc_matrix, 2, function(x){
  rownames(auc_matrix)[which.max(x)]}
)

table(scobj_epi$Cell2)
DimPlot(scobj_epi, group.by = 'Cell2', label = T, label.box = F, raster = F)

save(scobj_epi, file = '/data/chenrz/scseq_gse160269/scobj_epi.RData')

DimPlot(scobj_epi, group.by = 'Cell2', raster = T, label = T, repel = T,
        cols = colorRampPalette(brewer.pal(8, "Set2"))(8), pt.size = 1.2) + 
  labs(title = "Epithelia cells\n(n = 44,730)")


##### Fibro #####
## 重处理
scobj_fibro <- subset(scobj, Cell %in% c('Fibroblast', 'Pericytes'))
scobj_fibro <- NormalizeData(scobj_fibro)
scobj_fibro <- FindVariableFeatures(scobj_fibro, selection.method = "vst", nfeatures = 2000)
scobj_fibro <- ScaleData(scobj_fibro,
                         features = VariableFeatures(scobj_fibro))
scobj_fibro <- RunPCA(scobj_fibro, features = VariableFeatures(object = scobj_fibro), reduction.name = "pca")
scobj_fibro <- RunHarmony(scobj_fibro, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")

pcs_idx <- scobj_fibro@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj_fibro@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:10]))

scobj_fibro <- RunUMAP(scobj_fibro, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model = T)
scobj_fibro <- FindNeighbors(scobj_fibro, reduction = "harmony", dims = pcs_idx)
scobj_fibro <- FindClusters(scobj_fibro, resolution = 1, algorithm = 1)
scobj_fibro$seurat_clusters <- as.character(scobj_fibro$RNA_snn_res.1)
DimPlot(scobj_fibro, group.by = 'seurat_clusters', label = T, label.box = T)

# 标记基因分值
marker <- split(ref_anno$S6d_Fibroblasts$gene, ref_anno$S6d_Fibroblasts$cluster)
exprMatrix <- as(GetAssayData(scobj_fibro, assay = "RNA", slot = "scale.data"), "dgCMatrix")
cells_rankings <- AUCell_buildRankings(exprMatrix)
cells_AUC <- AUCell_calcAUC(marker, cells_rankings)
auc_matrix <- getAUC(cells_AUC)

scobj_fibro$Cell2 <- apply(auc_matrix, 2, function(x){
  rownames(auc_matrix)[which.max(x)]}
)

table(scobj_fibro$Cell2)
DimPlot(scobj_fibro, group.by = 'Cell2', label = T, label.box = T, raster = F)

save(scobj_fibro, file = '/data/chenrz/scseq_gse160269/scobj_fibro.RData')

DimPlot(scobj_fibro, group.by = 'Cell2', raster = T, label = T, repel = T,
        cols = colorRampPalette(brewer.pal(8, "Set3"))(9), pt.size = 1.3) + 
  labs(title = "Fibroblast cells\n(n = 40,315)")


##### Endo #####
## 重处理
scobj_endo <- subset(scobj, Cell == 'Endothelial')
scobj_endo <- NormalizeData(scobj_endo)
scobj_endo <- FindVariableFeatures(scobj_endo, selection.method = "vst", nfeatures = 2000)
scobj_endo <- ScaleData(scobj_endo,
                        features = VariableFeatures(scobj_endo))
scobj_endo <- RunPCA(scobj_endo, features = VariableFeatures(object = scobj_endo), reduction.name = "pca")
scobj_endo <- RunHarmony(scobj_endo, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")

pcs_idx <- scobj_endo@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj_endo@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:10]))

scobj_endo <- RunUMAP(scobj_endo, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model = T)
scobj_endo <- FindNeighbors(scobj_endo, reduction = "harmony", dims = pcs_idx)
scobj_endo <- FindClusters(scobj_endo, resolution = 1, algorithm = 1)
scobj_endo$seurat_clusters <- as.character(scobj_endo$RNA_snn_res.1)
DimPlot(scobj_endo, group.by = 'seurat_clusters', label = T, label.box = T)

# 标记基因分值
marker <- split(ref_anno$S6e_Endothelial$gene, ref_anno$S6e_Endothelial$cluster)
exprMatrix <- as(GetAssayData(scobj_endo, assay = "RNA", slot = "scale.data"), "dgCMatrix")
cells_rankings <- AUCell_buildRankings(exprMatrix)
cells_AUC <- AUCell_calcAUC(marker, cells_rankings)
auc_matrix <- getAUC(cells_AUC)

scobj_endo$Cell2 <- apply(auc_matrix, 2, function(x){
  rownames(auc_matrix)[which.max(x)]}
)

table(scobj_endo$Cell2)
DimPlot(scobj_endo, group.by = 'Cell2', label = T, label.box = T, raster = F)

save(scobj_endo, file = '/data/chenrz/scseq_gse160269/scobj_endo.RData')

DimPlot(scobj_endo, group.by = 'Cell2', raster = T, label = T, repel = T,
        cols = colorRampPalette(brewer.pal(8, "Set2"))(10), pt.size = 1.5) + 
  labs(title = "Endothelia cells\n(n = 11,267)")


##### Myeloid #####
## 重处理
scobj_mye <- subset(scobj, Cell == 'Myeloid')
scobj_mye <- NormalizeData(scobj_mye)
scobj_mye <- FindVariableFeatures(scobj_mye, selection.method = "vst", nfeatures = 2000)
scobj_mye <- ScaleData(scobj_mye,
                       features = VariableFeatures(scobj_mye))
scobj_mye <- RunPCA(scobj_mye, features = VariableFeatures(object = scobj_mye), reduction.name = "pca")
scobj_mye <- RunHarmony(scobj_mye, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")

pcs_idx <- scobj_mye@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj_mye@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:10]))

scobj_mye <- RunUMAP(scobj_mye, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model = T)
scobj_mye <- FindNeighbors(scobj_mye, reduction = "harmony", dims = pcs_idx)
scobj_mye <- FindClusters(scobj_mye, resolution = 1, algorithm = 1)
scobj_mye$seurat_clusters <- as.character(scobj_mye$RNA_snn_res.1)
DimPlot(scobj_mye, group.by = 'seurat_clusters', label = T, label.box = T)

# 标记基因分值
marker <- split(ref_anno$S6c_Myeloid$gene, ref_anno$S6c_Myeloid$cluster)
exprMatrix <- as(GetAssayData(scobj_mye, assay = "RNA", slot = "scale.data"), "dgCMatrix")
cells_rankings <- AUCell_buildRankings(exprMatrix)
cells_AUC <- AUCell_calcAUC(marker, cells_rankings)
auc_matrix <- getAUC(cells_AUC)

scobj_mye$Cell2 <- apply(auc_matrix, 2, function(x){
  rownames(auc_matrix)[which.max(x)]}
)

table(scobj_mye$Cell2)
DimPlot(scobj_mye, group.by = 'Cell2', label = T, label.box = T, raster = F)

save(scobj_mye, file = '/data/chenrz/scseq_gse160269/scobj_mye.RData')

DimPlot(scobj_mye, group.by = 'Cell2', raster = T, label = T, repel = T,
        cols = colorRampPalette(brewer.pal(8, "Set1"))(11), pt.size = 1.5) + 
  labs(title = "Myeloid cells\n(n = 19,273)")


##### Bcells #####
## 重处理
scobj_bcell <- subset(scobj, Cell == 'Bcell')
scobj_bcell <- NormalizeData(scobj_bcell)
scobj_bcell <- FindVariableFeatures(scobj_bcell, selection.method = "vst", nfeatures = 2000)
scobj_bcell <- ScaleData(scobj_bcell,
                         features = VariableFeatures(scobj_bcell))
scobj_bcell <- RunPCA(scobj_bcell, features = VariableFeatures(object = scobj_bcell), reduction.name = "pca")
scobj_bcell <- RunHarmony(scobj_bcell, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")

pcs_idx <- scobj_bcell@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj_bcell@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:10]))

scobj_bcell <- RunUMAP(scobj_bcell, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model = T)
scobj_bcell <- FindNeighbors(scobj_bcell, reduction = "harmony", dims = pcs_idx)
scobj_bcell <- FindClusters(scobj_bcell, resolution = 1, algorithm = 1)
scobj_bcell$seurat_clusters <- as.character(scobj_bcell$RNA_snn_res.1)
DimPlot(scobj_bcell, group.by = 'seurat_clusters', label = T, label.box = T)

# 标记基因分值
marker <- split(ref_anno$S6b_Bcells$gene, ref_anno$S6b_Bcells$cluster)
exprMatrix <- as(GetAssayData(scobj_bcell, assay = "RNA", slot = "scale.data"), "dgCMatrix")
cells_rankings <- AUCell_buildRankings(exprMatrix)
cells_AUC <- AUCell_calcAUC(marker, cells_rankings)
auc_matrix <- getAUC(cells_AUC)

scobj_bcell$Cell2 <- apply(auc_matrix, 2, function(x){
  rownames(auc_matrix)[which.max(x)]}
)

table(scobj_bcell$Cell2)
DimPlot(scobj_bcell, group.by = 'Cell2', label = T, label.box = T, raster = F)

save(scobj_bcell, file = '/data/chenrz/scseq_gse160269/scobj_bcell.RData')

DimPlot(scobj_bcell, group.by = 'Cell2', raster = T, label = T, repel = T,
        cols = colorRampPalette(brewer.pal(8, "Set1"))(11), pt.size = 1.5) + 
  labs(title = "B cells\n(n = 22,477)")


##### Tcells #####
## 重处理
scobj_tcell <- subset(scobj, Cell == 'Tcell')
scobj_tcell <- NormalizeData(scobj_tcell)
scobj_tcell <- FindVariableFeatures(scobj_tcell, selection.method = "vst", nfeatures = 2000)
scobj_tcell <- ScaleData(scobj_tcell,
                         features = VariableFeatures(scobj_tcell))
scobj_tcell <- RunPCA(scobj_tcell, features = VariableFeatures(object = scobj_tcell), reduction.name = "pca")
scobj_tcell <- RunHarmony(scobj_tcell, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")

pcs_idx <- scobj_tcell@reductions$harmony@stdev
names(pcs_idx) <- colnames(scobj_tcell@reductions$harmony@cell.embeddings)
pcs_idx_df <- as.data.frame(pcs_idx)
pcs_idx_df$rank <- rank(-pcs_idx_df$pcs_idx)
pcs_idx_df$name <- rownames(pcs_idx_df)
pcs_idx_df <- arrange(pcs_idx_df, rank)
ggplot(pcs_idx_df, aes(x = rank, y = pcs_idx)) +
  geom_point()
pcs_idx <- as.numeric(gsub('harmony_', '', pcs_idx_df$name[1:10]))

scobj_tcell <- RunUMAP(scobj_tcell, reduction = "harmony", dims = pcs_idx, reduction.name = "umap", return.model = T)
scobj_tcell <- FindNeighbors(scobj_tcell, reduction = "harmony", dims = pcs_idx)
scobj_tcell <- FindClusters(scobj_tcell, resolution = 1, algorithm = 1)
scobj_tcell$seurat_clusters <- as.character(scobj_tcell$RNA_snn_res.1)
DimPlot(scobj_tcell, group.by = 'seurat_clusters', label = T, label.box = T)

# 标记基因分值
marker <- split(ref_anno$S6a_TCells$gene, ref_anno$S6a_TCells$cluster)
exprMatrix <- as(GetAssayData(scobj_tcell, assay = "RNA", slot = "scale.data"), "dgCMatrix")
cells_rankings <- AUCell_buildRankings(exprMatrix)
cells_AUC <- AUCell_calcAUC(marker, cells_rankings)
auc_matrix <- getAUC(cells_AUC)

scobj_tcell$Cell2 <- apply(auc_matrix, 2, function(x){
  rownames(auc_matrix)[which.max(x)]}
)

table(scobj_tcell$Cell2)
DimPlot(scobj_tcell, group.by = 'Cell2', label = T, label.box = T, raster = F)

save(scobj_tcell, file = '/data/chenrz/scseq_gse160269/scobj_tcell.RData')

DimPlot(scobj_tcell, group.by = 'Cell2', raster = T, label = T, repel = T,
        cols = colorRampPalette(brewer.pal(8, "Set1"))(11), pt.size = 1.5) + 
  labs(title = "T cells\n(n = 69,278)")