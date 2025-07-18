library(Seurat)
library(readxl)
library(pbapply)
library(data.table)
library(openxlsx)
library(stringr)
library(dplyr)
library(ComplexHeatmap)
library(patchwork)
library(cowplot)
library(ggplot2)
library(ggExtra)
library(ggspavis)
library(ggalluvial)
library(clusterProfiler)
library(VisiumIO)
library(SpotSweeper)
library(scuttle)
library(spacexr)
library(MERINGUE)

load('/data/chenrz/LRinter.dataframe.RData')
load('/data/chenrz/scseq_gse160269/scobj.RData')
load('/data/chenrz/scseq_gse160269/emt_scobj_sccds.RData')

anno_reactome <- read.gmt('/data/chenrz/resource/reactome.v2024.1.Hs.symbols.gmt')
anno_hallmark <- read.gmt('/data/chenrz/resource/hallmark.v2023.2.Hs.symbols.gmt')
anno_go <- read.gmt('/data/chenrz/resource/go.v2024.1.Hs.symbols.gmt')

sheets <- excel_sheets('/data/chenrz/scseq_gse160269/41467_2021_25539_MOESM9_ESM.xlsx')
ref_anno <- lapply(sheets, function(x) read.xlsx('/data/chenrz/scseq_gse160269/41467_2021_25539_MOESM9_ESM.xlsx', sheet = x, startRow = 3))
names(ref_anno) <- sheets
ref_anno[['Epithelia']] <- read.xlsx('/data/chenrz/scseq_gse160269/41467_2021_25539_MOESM8_ESM.xlsx', sheet = 1, startRow = 3)

ncbi <- read.delim('/data/chenrz/resource/Homo_sapiens.gene_info', check.names = F) # 基因文件
ncbi$Ensembl <- unlist(pbapply(ncbi, 1, function(x){
  index <- gregexpr('ENSG\\d+', x[6])
  a = regmatches(x[6], index)[[1]]
  if (index[[1]][1] == -1) {
    a = '-'
  } else if (index[[1]][1] > 1) {
    a = paste0(a, collapse = '|')
  }
  return(a)
}))
ncbi_f <- filter(ncbi, ncbi$Ensembl != '-')

options(future.globals.maxSize = 8000 * 1024^2)


# marker -----------------------------------------------------------------------
marker <- rbind(arrange(ref_anno[["S6a_TCells"]], desc(avg_logFC)), 
                arrange(ref_anno[["S6b_Bcells"]], desc(avg_logFC)),
                arrange(ref_anno[["S6c_Myeloid"]], desc(avg_logFC)), 
                arrange(ref_anno[["S6e_Endothelial"]], desc(avg_logFC)), 
                arrange(ref_anno[["S6d_Fibroblasts"]], desc(avg_logFC)))
marker_1 <- reshape2::melt(as.matrix(ref_anno[["Epithelia"]]), varnames = c("Row", "Col"), value.name = "Value")
marker_2 <- matrix(nrow = nrow(marker_1), ncol = ncol(marker))
marker_2[,1] <- as.character(marker_1[,2])
marker_2[,2] <- marker_1[,3]
colnames(marker_2) <- colnames(marker)
marker <- rbind(marker, marker_2)
marker <- group_by(marker, cluster) %>% slice_head(n = 5)

marker_emt <- FindAllMarkers(scobj_epi_mes, group.by = 'EMT_state')


# ViusumHD ---------------------------------------------------------------------
# ViusumHD对应病理切片240211911.mrxs
# S047801_2 sct会卡死
S047801_16 <- Load10X_Spatial(data.dir = '/data/chenrz/stseq_data/S047801/', bin.size = 16)
S047801_8 <- Load10X_Spatial(data.dir = '/data/chenrz/stseq_data/S047801/', bin.size = 8)

S047801_16 <- S047801_16[,which(colSums(GetAssayData(S047801_16))!=0)]
S047801_8 <- S047801_8[,which(colSums(GetAssayData(S047801_8))!=0)]

##### 标准化 #####
S047801_16 <- SCTransform(S047801_16, assay = "Spatial.016um", verbose = TRUE)
S047801_16 <- ScaleData(S047801_16, features = rownames(S047801_16))

S047801_8 <- SCTransform(S047801_8, assay = "Spatial.008um", verbose = TRUE)

save(S047801_16, file = '/data/chenrz/stseq_data/S047801/S047801_016um.RData')
save(S047801_8, file = '/data/chenrz/stseq_data/S047801/S047801_008um.RData')

##### QC #####
S047801_16_spe <- TENxVisiumHD(
  spacerangerOut='/data/chenrz/stseq_data/S047801/', 
  processing="filtered",
  format="h5", 
  images="lowres", 
  bin_size="016") |>
  import()

gs <- rowData(S047801_16_spe)$Symbol
rownames(S047801_16_spe) <- make.unique(gs)
colnames(S047801_16_spe) <- S047801_16_spe$barcode

mt <- grepl("^MT-", rownames(S047801_16_spe))
S047801_16_spe <- addPerCellQCMetrics(S047801_16_spe, subsets=list(mt=mt)) # 线粒体基因比例

# determine outliers based on 
# - low log-library size
# - few uniquely detected features
# - high mitochondrial count fraction
S047801_16_spe <- localOutliers(S047801_16_spe, metric="sum", direction="lower", log=TRUE)
S047801_16_spe <- localOutliers(S047801_16_spe, metric="detected", direction="lower", log=TRUE)
S047801_16_spe <- localOutliers(S047801_16_spe, metric="subsets_mt_percent", direction="higher", log=TRUE)

S047801_16_spe$discard <- 
  S047801_16_spe$sum_outliers | 
  S047801_16_spe$detected_outliers | 
  S047801_16_spe$subsets_mt_percent_outliers

plotCoords(S047801_16_spe, annotate="sum_outliers") + ggtitle("low_lib_size") +
  plotCoords(S047801_16_spe, annotate="detected_outliers") + ggtitle("low_n_features") +
  plotCoords(S047801_16_spe, annotate="discard") + ggtitle("discard") +
  plot_layout(nrow=1, guides="collect") & theme(
    plot.title=element_text(hjust=0.5),
    legend.key.size=unit(0, "lines")) &
  guides(col=guide_legend(override.aes=list(size=3))) & 
  scale_color_manual("discard", values=c("lavender", "purple"))

table(S047801_16_spe$discard) # 3118

# 去掉低质量spot
S047801_16 <- S047801_16[,setdiff(colnames(S047801_16), colnames(S047801_16_spe)[S047801_16_spe$discard])]
save(S047801_16_spe, file = '/data/chenrz/S047801_spe_016um.RData')


##### 细胞类型反卷积 #####
# 参考数据集
ref <- subset(scobj, downsample = 2000) # 下采样
ref$Cell2 <- sub("/","_",ref$Cell2)
sort(table(ref@meta.data[["Cell2"]]))

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$Cell2)
cluster <- droplevels(cluster)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA 
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# 反卷积数据集
counts <- GetAssayData(S047801_16, layer = 'counts')
coords <- GetTissueCoordinates(S047801_16)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

# 反卷积
RCTD <- create.RCTD(query, reference, max_cores = 10)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

save(RCTD, file = '/data/chenrz/stseq_data/S047801/RCTD_016um.RData')

# 可视化关注细胞类型
S047801_16 <- AddMetaData(S047801_16, metadata = RCTD@results$results_df)
SpatialDimPlot(S047801_16, group.by = "first_type", pt.size.factor = 3.5, image.alpha = 1) + labs(fill = 'Cell type')

# 上皮细胞
S047801_16_epi <- subset(S047801_16, first_type %in% c('AP', 'Cycling', 'Mes', 'Mucosal', 'Oxd', 'QP', 'Stress', 'TD'))
S047801_16_epi$first_type <- as.character(S047801_16_epi$first_type)

mycols <- unique(colorRampPalette(brewer.pal(8, "Set2"))(8))
names(mycols) <- unique(S047801_16_epi$first_type)
res_abundance <- as.data.frame(table(S047801_16_epi$first_type)/length(S047801_16_epi$first_type))
res_abundance$Var2 <- 'A'
res_abundance <- arrange(res_abundance, desc(Freq))
res_abundance$Var1 <- factor(res_abundance$Var1, levels = res_abundance$Var1)
S047801_16_epi$first_type <- factor(S047801_16_epi$first_type, levels = res_abundance$Var1)

p1 <- SpatialDimPlot(S047801_16_epi, group.by = "first_type", pt.size.factor = 4, image.alpha = .5, cols = mycols) + 
  labs(fill = 'Epithelial cell')
p1_2 <- ggplot(res_abundance, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_stratum(width = 0.7, color='white', size=0.6) + 
  geom_alluvium(alpha = 0.3, width = 0.7, curve_type = "linear") +
  theme_void() +
  scale_fill_manual(values = mycols) + 
  labs(x = NA, y = NA, fill = "Cell type") +
  theme(, legend.position = "none") +
  scale_y_reverse() +
  coord_flip()

p1 / p1_2 + plot_layout(ncol = 1, heights = c(10, 1))  


# 间质细胞
S047801_16_fib <- subset(S047801_16, first_type %in% c('CAF1', 'CAF2', 'CAF3', 'CAF4', 'NAF1', 'NAF2', 'NMF', 'Pericyte', 'VSMC', 'FRC'))
S047801_16_fib$first_type <- as.character(S047801_16_fib$first_type)

mycols <- unique(colorRampPalette(brewer.pal(8, "Set2"))(10))
names(mycols) <- unique(S047801_16_fib$first_type)
res_abundance <- as.data.frame(table(S047801_16_fib$first_type)/length(S047801_16_fib$first_type))
res_abundance$Var2 <- 'A'
res_abundance <- arrange(res_abundance, desc(Freq))
res_abundance$Var1 <- factor(res_abundance$Var1, levels = res_abundance$Var1)
S047801_16_fib$first_type <- factor(S047801_16_fib$first_type, levels = res_abundance$Var1)

p2 <- SpatialDimPlot(S047801_16_fib, group.by = "first_type", pt.size.factor = 4, image.alpha = .5, cols = mycols) + 
  labs(fill = 'Fibroblast cell')
p2_2 <- ggplot(res_abundance, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_stratum(width = 0.7, color='white', size=0.6) + 
  geom_alluvium(alpha = 0.3, width = 0.7, curve_type = "linear") +
  theme_void() +
  scale_fill_manual(values = mycols) + 
  labs(x = NA, y = NA, fill = "Cell type") +
  theme(, legend.position = "none") +
  scale_y_reverse() +
  coord_flip()

p2 / p2_2 + plot_layout(ncol = 1, heights = c(10, 1)) 


# 内皮细胞
S047801_16_endo <- subset(S047801_16, first_type %in% c('NEC1', 'NEC2', 'TEC1', 'TEC2', 'TEC3', 'TEC4', 'TEC5'))
S047801_16_endo$first_type <- as.character(S047801_16_endo$first_type)

mycols <- unique(colorRampPalette(brewer.pal(8, "Set2"))(7))
names(mycols) <- unique(S047801_16_endo$first_type)
res_abundance <- as.data.frame(table(S047801_16_endo$first_type)/length(S047801_16_endo$first_type))
res_abundance$Var2 <- 'A'
res_abundance <- arrange(res_abundance, desc(Freq))
res_abundance$Var1 <- factor(res_abundance$Var1, levels = res_abundance$Var1)
S047801_16_endo$first_type <- factor(S047801_16_endo$first_type, levels = res_abundance$Var1)

p3 <- SpatialDimPlot(S047801_16_endo, group.by = "first_type", pt.size.factor = 8, image.alpha = .5, cols = mycols) + 
  labs(fill = 'Endothelial cell')
p3_2 <- ggplot(res_abundance, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_stratum(width = 0.7, color='white', size=0.6) + 
  geom_alluvium(alpha = 0.3, width = 0.7, curve_type = "linear") +
  theme_void() +
  scale_fill_manual(values = mycols) + 
  labs(x = NA, y = NA, fill = "Cell type") +
  theme(, legend.position = "none") +
  scale_y_reverse() +
  coord_flip()

p3 / p3_2 + plot_layout(ncol = 1, heights = c(10, 1)) 

# 淋巴细胞
S047801_16_lym <- subset(S047801_16, first_type %in% c(unique(ref_anno[["S6a_TCells"]][["cluster"]]), 
                                                       unique(ref_anno[["S6b_Bcells"]][["cluster"]]),
                                                       unique(ref_anno[["S6c_Myeloid"]][["cluster"]])))
S047801_16_lym$first_type <- as.character(S047801_16_lym$first_type)

mycols <- unique(colorRampPalette(brewer.pal(8, "Set1"))(25))
names(mycols) <- unique(S047801_16_lym$first_type)
res_abundance <- as.data.frame(table(S047801_16_lym$first_type)/length(S047801_16_lym$first_type))
res_abundance$Var2 <- 'A'
res_abundance <- arrange(res_abundance, desc(Freq))
res_abundance$Var1 <- factor(res_abundance$Var1, levels = res_abundance$Var1)
S047801_16_lym$first_type <- factor(S047801_16_lym$first_type, levels = res_abundance$Var1)

p4 <- SpatialDimPlot(S047801_16_lym, group.by = "first_type", pt.size.factor = 8, image.alpha = .5, cols = mycols) + 
  labs(fill = 'Immune cell')
p4_2 <- ggplot(res_abundance, aes(x = Var2, y = Freq, fill = Var1, stratum = Var1, alluvium = Var1)) +
  geom_stratum(width = 0.7, color='white', size=0.6) + 
  geom_alluvium(alpha = 0.3, width = 0.7, curve_type = "linear") +
  theme_void() +
  scale_fill_manual(values = mycols) + 
  labs(x = NA, y = NA, fill = "Cell type") +
  theme(, legend.position = "none") +
  scale_y_reverse() +
  coord_flip()

p4 / p4_2 + plot_layout(ncol = 1, heights = c(10, 1)) 


##### EMT类型反卷积 #####
# 参考数据集
ref <- subset(scobj_epi_mes, downsample = 2000) # 下采样
ref$EMT_state <- sub("/","_",ref$EMT_state)
sort(table(ref$EMT_state))
table(ref$EMT_state, ref$seurat_clusters)

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$EMT_state)
cluster <- droplevels(cluster)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA 
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# 反卷积数据集
S047801_16_mes <- subset(S047801_16, first_type == 'Mes' | second_type == 'Mes')

counts <- GetAssayData(S047801_16_mes, layer = 'counts')
coords <- GetTissueCoordinates(S047801_16_mes)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

# 反卷积
RCTD_mes <- create.RCTD(query, reference, max_cores = 3)
RCTD_mes <- run.RCTD(RCTD_mes, doublet_mode = "doublet")

save(RCTD_mes, file = '/data/chenrz/stseq_data/S047801/RCTD_mes_016.RData')

table(RCTD_mes@results$results_df$first_type)
table(RCTD_mes@results$results_df$second_type)

table(RCTD_mes@results$results_df$spot_class, RCTD_mes@results$results_df$first_type)
table(RCTD_mes@results$results_df$spot_class, RCTD_mes@results$results_df$second_type)

# 可视化关注细胞类型
S047801_16_mes <- AddMetaData(S047801_16_mes, metadata = RCTD_mes@results$results_df)
S047801_16_mes$type <- as.character(S047801_16_mes$first_type)
S047801_16_mes$type <- ifelse(S047801_16_mes$second_type == 'EMT_stable', 'EMT_stable', S047801_16_mes$type) # 手动调整结果

SpatialDimPlot(S047801_16_mes, group.by = "first_type", pt.size.factor = 4, image.alpha = 0)
SpatialDimPlot(S047801_16_mes, group.by = "second_type", pt.size.factor = 4, image.alpha = 0)
SpatialDimPlot(S047801_16_mes, group.by = "type", pt.size.factor = 4, image.alpha = 0) # 符合预期

table(S047801_16_mes$type)

# 上皮间充质分值分布
epi <- c('KRT14', 'KRT17', 'KRT6A', 'KRT5', 'KRT19', 'KRT8', 'KRT16', 'KRT18', 'KRT6B', 'KRT15', 'KRT6C', 'KRTCAP3', 'SFN', 'EPCAM')
mes <- c('VIM', 'CDH2', 'FOXC2', 'SNAI1', 'SNAI2', 'TWIST1', 'FN1', 'ITGB6', 'MMP2', 'MMP3', 'MMP9', 'SOX10', 'GSC', 'ZEB1', 'ZEB2', 'TWIST2')
S047801_16_mes <- AddModuleScore(S047801_16_mes, features = list(epi), name = 'Epi')
S047801_16_mes <- AddModuleScore(S047801_16_mes, features = list(mes), name = 'Mes')
S047801_16_mes$EMT <- S047801_16_mes$Mes1 - S047801_16_mes$Epi1

SpatialFeaturePlot(S047801_16_mes, features = 'Epi1', pt.size.factor = 4, image.alpha = 0)
SpatialFeaturePlot(S047801_16_mes, features = 'Mes1', pt.size.factor = 4, image.alpha = 0)


##### 组织区域划分 #####
label <- fread('/data/chenrz/STARCH/s047801/labels_s047801.csv')
position <- fread('/data/chenrz/STARCH/s047801/tissue_positions_list.csv')
position$V1 <- paste0(sprintf("%.1f", position$array_row), 'x', sprintf("%.1f", position$array_col))
position <- left_join(position, label, by = 'V1')

# 定义肿瘤间质交界区
position_f <- filter(position, V2 == 1) %>% group_by(array_col) %>% arrange(array_row) %>% slice_head(n = 3)
position$V2 <- ifelse(position$barcode %in% position_f$barcode, 3, position$V2)
S047801_16$CNV <- position$V2[match(colnames(S047801_16), position$barcode)]
SpatialDimPlot(S047801_16, group.by = 'CNV', pt.size.factor = 4)

# 0 - 侵袭区域
# 1 - 肿瘤区域
# 2 - 间质和正常上皮区域
# 3 - 肿瘤间质交界区域

S047801_16_mes$CNV <- position$V2[match(colnames(S047801_16_mes), position$barcode)]
S047801_16_mes$CNV <- ifelse(S047801_16_mes$CNV == 0, 'Normal_Tumor',
                             ifelse(S047801_16_mes$CNV == 1, 'Tumor',
                                    ifelse(S047801_16_mes$CNV == 2, 'Stroma', 'Tumor_Stroma')))
S047801_16_mes$CNV <- factor(S047801_16_mes$CNV, levels = c('Normal_Tumor', 'Tumor', 'Tumor_Stroma', 'Stroma'))

SpatialDimPlot(S047801_16_mes, group.by = 'CNV', pt.size.factor = 6, image.alpha = 0) +
  labs(fill = 'Region') +
  theme_bw() +
  labs(x = '', y = '')

SpatialDimPlot(S047801_16_mes, group.by = "type", pt.size.factor = 6, image.alpha = 0) + 
  labs(fill = 'Mes cell') +
  theme_bw() +
  labs(x = '', y = '')

# 定量
ht <- as.matrix(table(S047801_16_mes$first_type, S047801_16_mes$CNV))
ht <- as.data.frame(apply(ht, 2, function(x){x/sum(x)}))
mark_mat <- matrix(FALSE, nrow = nrow(ht), ncol = ncol(ht))
for (i in 1:ncol(ht)) {
  mark_mat[which.max(ht[,i]),i] <- TRUE
}

Heatmap(ht, show_row_dend = F, show_column_dend = F, 
        cluster_rows = F, cluster_columns = F,
        col = circlize::colorRamp2(c(0, 0.5), c('white' ,'red')),
        border = TRUE, 
        heatmap_legend_param = list(title = "Proportion"),
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (mark_mat[i, j]) {
            grid.points(
              x = x, y = y, 
              pch = 8,  # 星号形状
              size = unit(4, "mm"),
              gp = gpar(col = "black")
            )
          }
        }
)


# Steroseq ----------------------------------------------------------------
read_steroseq = function(bin_path){
  count <- Read10X(paste0(bin_path, '/filtered_feature_bc_matrix'))
  stseq <- CreateSeuratObject(counts = count, assay = 'Spatial')
  image.list <- mapply(
    Read10X_Image,
    paste0(bin_path, '/spatial'),
    assay = "Spatial",
    slice = "slice1",
    MoreArgs = list(filter.matrix = TRUE)
  )
  names(image.list) <- c("slice1")
  stseq@images <- image.list
  stseq@images[["slice1"]]@boundaries[["centroids"]]@cells <- as.character(format(as.numeric(stseq@images[["slice1"]]@boundaries[["centroids"]]@cells), scientific = FALSE))
  
  genemap <- ncbi_f$Symbol
  names(genemap) <- ncbi_f$Ensembl
  gene <- intersect(rownames(stseq), ncbi_f$Ensembl)
  stseq <- stseq[gene,]
  rownames(stseq) <- recode(rownames(stseq), !!!genemap)
  
  stseq <- stseq[,unname(which(colSums(GetAssayData(stseq))!=0))]
  
  return(stseq)
  
}

# E1(刘振鑫 - _13_1_103028.svs)
E1_100 <- read_steroseq('/data/chenrz/stseq_data/E1/E1/bin100/')
E1_50 <- read_steroseq('/data/chenrz/stseq_data/E1/E1/bin50/')
E1_20 <- read_steroseq('/data/chenrz/stseq_data/E1/E1/bin20/')

# E2(刘三妹 - E012101.svs)
E2_100 <- read_steroseq('/data/chenrz/stseq_data/E2/E2/bin100/')
E2_50 <- read_steroseq('/data/chenrz/stseq_data/E2/E2/bin50/')
E2_20 <- read_steroseq('/data/chenrz/stseq_data/E2/E2/bin20/')


##### 标准化 #####
# E1
E1_100 <- SCTransform(E1_100, assay = "Spatial", verbose = TRUE)
E1_50 <- SCTransform(E1_50, assay = "Spatial", verbose = TRUE)
E1_20 <- SCTransform(E1_20, assay = "Spatial", verbose = TRUE)

save(E1_100, file = '/data/chenrz/stseq_data/E1/E1/E1_100bin.RData')
save(E1_50, file = '/data/chenrz/stseq_data/E1/E1/E1_50bin.RData')
save(E1_20, file = '/data/chenrz/stseq_data/E1/E1/E1_20bin.RData')

# E2
E2_100 <- SCTransform(E2_100, assay = "Spatial", verbose = TRUE)
E2_50 <- SCTransform(E2_50, assay = "Spatial", verbose = TRUE)
E2_20 <- SCTransform(E2_20, assay = "Spatial", verbose = TRUE)

save(E2_100, file = '/data/chenrz/stseq_data/E2/E2/E2_100bin.RData')
save(E2_50, file = '/data/chenrz/stseq_data/E2/E2/E2_50bin.RData')
save(E2_20, file = '/data/chenrz/stseq_data/E2/E2/E2_20bin.RData')


##### QC #####
qc_steroseq = function(bin_path){
  stero_spe <- TENxVisium(
    resources = paste0(bin_path, '/filtered_feature_bc_matrix'),
    spatialResource = paste0(bin_path, '/spatial'),
    sample_id = 'sample',
    images = 'lowres') |>
    import()
  gs <- rowData(stero_spe)$Symbol
  rownames(stero_spe) <- make.unique(gs)
  colnames(stero_spe) <- stero_spe$barcode
  
  mt <- grepl("^MT-", rownames(stero_spe))
  stero_spe <- addPerCellQCMetrics(stero_spe, subsets=list(mt=mt))
  stero_spe <- localOutliers(stero_spe, metric="sum", direction="lower", log=TRUE)
  stero_spe <- localOutliers(stero_spe, metric="detected", direction="lower", log=TRUE)
  stero_spe <- localOutliers(stero_spe, metric="subsets_mt_percent", direction="higher", log=TRUE)
  
  stero_spe$discard <- 
    stero_spe$sum_outliers | 
    stero_spe$detected_outliers | 
    stero_spe$subsets_mt_percent_outliers
  
  return(stero_spe)
  
}

# qc
E1_100_spe_qc <- qc_steroseq('/data/chenrz/stseq_data/E1/E1/bin100/')
E1_50_spe_qc <- qc_steroseq('/data/chenrz/stseq_data/E1/E1/bin50/')
E1_20_spe_qc <- qc_steroseq('/data/chenrz/stseq_data/E1/E1/bin20/')

save(E1_100_spe_qc, file = '/data/chenrz/stseq_data/E1/E1/E1_spe_100bin.RData')
save(E1_50_spe_qc, file = '/data/chenrz/stseq_data/E1/E1/E1_spe_50bin.RData')
save(E1_20_spe_qc, file = '/data/chenrz/stseq_data/E1/E1/E1_spe_20bin.RData')

E2_100_spe_qc <- qc_steroseq('/data/chenrz/stseq_data/E2/E2/bin100/')
E2_50_spe_qc <- qc_steroseq('/data/chenrz/stseq_data/E2/E2/bin50/')
E2_20_spe_qc <- qc_steroseq('/data/chenrz/stseq_data/E2/E2/bin20/')

save(E2_100_spe_qc, file = '/data/chenrz/stseq_data/E2/E2/E2_spe_100bin.RData')
save(E2_50_spe_qc, file = '/data/chenrz/stseq_data/E2/E2/E2_spe_50bin.RData')
save(E2_20_spe_qc, file = '/data/chenrz/stseq_data/E2/E2/E2_spe_20bin.RData')

# 可视化
qc <- E2_20_spe_qc
plotCoords(qc, annotate="sum_outliers") + ggtitle("low_lib_size") +
  plotCoords(qc, annotate="detected_outliers") + ggtitle("low_n_features") +
  plotCoords(qc, annotate="discard") + ggtitle("discard") +
  plot_layout(nrow=1, guides="collect") & theme(
    plot.title=element_text(hjust=0.5),
    legend.key.size=unit(0, "lines")) &
  guides(col=guide_legend(override.aes=list(size=3))) & 
  scale_color_manual("discard", values=c("lavender", "purple"))

# 去掉低质量spot
E1_100 <- E1_100[,setdiff(colnames(E1_100), colnames(E1_100_spe_qc)[E1_100_spe_qc$discard])]
E1_50 <- E1_50[,setdiff(colnames(E1_50), colnames(E1_50_spe_qc)[E1_50_spe_qc$discard])]
E1_20 <- E1_20[,setdiff(colnames(E1_20), colnames(E1_20_spe_qc)[E1_20_spe_qc$discard])]

E2_100 <- E2_100[,setdiff(colnames(E2_100), colnames(E2_100_spe_qc)[E2_100_spe_qc$discard])]
E2_50 <- E2_50[,setdiff(colnames(E2_50), colnames(E2_50_spe_qc)[E2_50_spe_qc$discard])]
E2_20 <- E2_20[,setdiff(colnames(E2_20), colnames(E2_20_spe_qc)[E2_20_spe_qc$discard])]


##### 细胞类型反卷积 #####
deconv_main = function(refdata, spatial) {
  
  # 参考数据集
  ref <- subset(refdata, downsample = 2000) # 下采样
  ref$Cell2 <- sub("/","_",ref$Cell2)
  counts <- ref[["RNA"]]$counts
  cluster <- as.factor(ref$Cell2)
  cluster <- droplevels(cluster)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA 
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)
  
  # 反卷积数据集
  counts <- GetAssayData(spatial, layer = 'counts')
  coords <- GetTissueCoordinates(spatial)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  # 反卷积
  RCTD <- create.RCTD(query, reference, max_cores = 10)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  return(RCTD)
  
}

# steroseq的20bin对应visiumhd的16um
RCTD_E1_50 <- deconv_main(scobj, E1_50)
RCTD_E2_50 <- deconv_main(scobj, E2_50)

RCTD_E1_20 <- deconv_main(scobj, E1_20)
RCTD_E2_20 <- deconv_main(scobj, E2_20)

save(RCTD_E1_50, file = '/data/chenrz/stseq_data/E1/E1/RCTD_E1_50bin.RData')
save(RCTD_E2_50, file = '/data/chenrz/stseq_data/E2/E2/RCTD_E2_50bin.RData')

save(RCTD_E1_20, file = '/data/chenrz/stseq_data/E1/E1/RCTD_E1_20bin.RData')
save(RCTD_E2_20, file = '/data/chenrz/stseq_data/E2/E2/RCTD_E2_20bin.RData')

gc()

# 可视化
E1_50 <- AddMetaData(E1_50, metadata = RCTD_E1_50@results$results_df)
E2_50 <- AddMetaData(E2_50, metadata = RCTD_E2_50@results$results_df)

SpatialDimPlot(E1_50, group.by = "first_type", pt.size.factor = 2.5, image.alpha = 1) + labs(fill = 'Cell type')
SpatialDimPlot(E2_50, group.by = "first_type", pt.size.factor = 2.5, image.alpha = 1) + labs(fill = 'Cell type')

E1_20 <- AddMetaData(E1_20, metadata = RCTD_E1_20@results$results_df)
E2_20 <- AddMetaData(E2_20, metadata = RCTD_E2_20@results$results_df)

SpatialDimPlot(E1_20, group.by = "first_type", pt.size.factor = 2.5, image.alpha = 1) + labs(fill = 'Cell type')
SpatialDimPlot(E2_20, group.by = "first_type", pt.size.factor = 2.5, image.alpha = 1) + labs(fill = 'Cell type')


##### EMT类型反卷积 #####
deconv_emt = function(refdata, spatial) {
  
  # 参考数据集
  ref <- subset(refdata, downsample = 2000) # 下采样
  ref$EMT_state <- sub("/","_",ref$EMT_state)
  counts <- ref[["RNA"]]$counts
  cluster <- as.factor(ref$EMT_state)
  cluster <- droplevels(cluster)
  names(cluster) <- colnames(ref)
  nUMI <- ref$nCount_RNA 
  names(nUMI) <- colnames(ref)
  reference <- Reference(counts, cluster, nUMI)
  
  # 反卷积数据集
  counts <- GetAssayData(spatial, layer = 'counts')
  coords <- GetTissueCoordinates(spatial)
  colnames(coords) <- c("x", "y")
  coords[is.na(colnames(coords))] <- NULL
  query <- SpatialRNA(coords, counts, colSums(counts))
  
  # 反卷积
  RCTD <- create.RCTD(query, reference, max_cores = 3)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
  
  return(RCTD)
  
}

E1_50_mes <- subset(E1_50, first_type == 'Mes' | second_type == 'Mes')
E1_50_mes <- SCTransform(E1_50_mes, assay = "Spatial", verbose = TRUE)
RCTD_E1_50_mes <- deconv_emt(scobj_epi_mes, E1_50_mes)

E2_50_mes <- subset(E2_50, first_type == 'Mes' | second_type == 'Mes')
E2_50_mes <- SCTransform(E2_50_mes, assay = "Spatial", verbose = TRUE)
RCTD_E2_50_mes <- deconv_emt(scobj_epi_mes, E2_50_mes)

E1_20_mes <- subset(E1_20, first_type == 'Mes' | second_type == 'Mes')
E1_20_mes <- SCTransform(E1_20_mes, assay = "Spatial", verbose = TRUE)
RCTD_E1_20_mes <- deconv_emt(scobj_epi_mes, E1_20_mes)

E2_20_mes <- subset(E2_20, first_type == 'Mes' | second_type == 'Mes')
E2_20_mes <- SCTransform(E2_20_mes, assay = "Spatial", verbose = TRUE)
RCTD_E2_20_mes <- deconv_emt(scobj_epi_mes, E2_20_mes)

save(RCTD_E1_20_mes, file = '/data/chenrz/stseq_data/E1/E1/RCTD_E1_mes_20bin.RData')
save(RCTD_E2_20_mes, file = '/data/chenrz/stseq_data/E2/E2/RCTD_E2_mes_20bin.RData')


# 可视化
E1_50_mes <- AddMetaData(E1_50_mes, metadata = RCTD_E1_50_mes@results$results_df)
table(E1_50_mes$second_type)
SpatialDimPlot(E1_50_mes, group.by = "second_type", pt.size.factor = 3, image.alpha = 1)

E2_50_mes <- AddMetaData(E2_50_mes, metadata = RCTD_E2_50_mes@results$results_df)
table(E2_50_mes$second_type)
SpatialDimPlot(E2_50_mes, group.by = "second_type", pt.size.factor = 3, image.alpha = 1)

####
E1_20_mes <- AddMetaData(E1_20_mes, metadata = RCTD_E1_20_mes@results$results_df)
table(E1_20_mes$second_type)
SpatialDimPlot(E1_20_mes, group.by = "second_type", pt.size.factor = 4, image.alpha = 1)

E2_20_mes <- AddMetaData(E2_20_mes, metadata = RCTD_E2_20_mes@results$results_df)
table(E2_20_mes$second_type)
SpatialDimPlot(E2_20_mes, group.by = "second_type", pt.size.factor = 4, image.alpha = 1)

####
## steroseq反卷积出来很多区域没有细胞
## E1存活，E2死亡
## E1和E2的EMT_diff含量有很大差别，能否作为一个论点？
####


##### 组织区域划分 #####
# E1
label <- fread('/data/chenrz/STARCH/e1/labels_e1.csv')
position <- fread('/data/chenrz/STARCH/e1/tissue_positions_list.csv')
position$V7 <- paste0(sprintf("%.1f", position$V3), 'x', sprintf("%.1f", position$V4))
position <- left_join(position, label, by = c('V7' = 'V1'))
E1_20$CNV <- position$V2.y[match(colnames(E1_20), position$V1)]

SpatialDimPlot(E1_20, group.by = 'CNV', pt.size.factor = 3, image.alpha = 0) +
  labs(fill = 'Region') +
  theme_bw() +
  labs(x = '', y = '')

# E2
label <- fread('/data/chenrz/STARCH/E2/labels_e2.csv')
position <- fread('/data/chenrz/STARCH/e2/tissue_positions_list.csv')
position$V7 <- paste0(sprintf("%.1f", position$V3), 'x', sprintf("%.1f", position$V4))
position <- left_join(position, label, by = c('V7' = 'V1'))
E2_20$CNV <- position$V2.y[match(colnames(E2_20), position$V1)]

SpatialDimPlot(E2_20, group.by = 'CNV', pt.size.factor = 3, image.alpha = 0) +
  labs(fill = 'Region') +
  theme_bw() +
  labs(x = '', y = '')
