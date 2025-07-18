library(Seurat)
library(readxl)
library(pbapply)
library(data.table)
library(openxlsx)
library(stringr)
library(dplyr)
library(ComplexHeatmap)
library(patchwork)
library(harmony)
library(cowplot)
library(ggplot2)
library(ggExtra)
library(ggspavis)
library(ggalluvial)
library(clusterProfiler)
library(mclust)
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

norm = function(x){log2(x - min(x) + 1)}


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

epi <- c('KRT14', 'KRT17', 'KRT6A', 'KRT5', 'KRT19', 'KRT8', 'KRT16', 'KRT18', 'KRT6B', 'KRT15', 'KRT6C', 'KRTCAP3', 'SFN', 'EPCAM')
mes <- c('VIM', 'CDH2', 'FOXC2', 'SNAI1', 'SNAI2', 'TWIST1', 'FN1', 'ITGB6', 'MMP2', 'MMP3', 'MMP9', 'SOX10', 'GSC', 'ZEB1', 'ZEB2', 'TWIST2')

marker_emt <- FindAllMarkers(scobj_epi_mes, group.by = 'EMT_state')
marker_emt <- arrange(marker_emt, desc(avg_log2FC))


# ViusumHD ---------------------------------------------------------------------
load('/data/chenrz/stseq_data/S047801/S047801_016um.RData')
load('/data/chenrz/stseq_data/S047801/S047801_spe_016um.RData')
load('/data/chenrz/stseq_data/S047801/RCTD_016um.RData')
load('/data/chenrz/stseq_data/S047801/RCTD_mes_016.RData')

# QC
S047801_16 <- S047801_16[,setdiff(colnames(S047801_16), colnames(S047801_16_spe)[S047801_16_spe$discard])]

# Deconv
S047801_16 <- AddMetaData(S047801_16, metadata = RCTD@results$results_df)
S047801_16 <- ScaleData(S047801_16, features = rownames(S047801_16))

S047801_16_mes <- subset(S047801_16, first_type == 'Mes' | second_type == 'Mes')
S047801_16_mes <- ScaleData(S047801_16_mes, features = rownames(S047801_16_mes))
S047801_16_mes <- AddMetaData(S047801_16_mes, metadata = RCTD_mes@results$results_df)
S047801_16_mes$type <- as.character(S047801_16_mes$first_type)
S047801_16_mes$type <- ifelse(S047801_16_mes$second_type == 'EMT_stable', 'EMT_stable', S047801_16_mes$type) # 手动调整结果

S047801_16_mes <- AddModuleScore(S047801_16_mes, features = list(epi), name = 'Epi')
S047801_16_mes <- AddModuleScore(S047801_16_mes, features = list(mes), name = 'Mes')
S047801_16_mes$EMT <- S047801_16_mes$Mes1 - S047801_16_mes$Epi1

ggplot(S047801_16_mes@meta.data, aes(x = type, y = EMT, fill = type)) + 
  geom_boxplot() + 
  scale_x_discrete(limits = c("EMT_early", "EMT_stable", "EMT_diff")) +
  theme_bw() +
  coord_flip() +
  labs(x = 'Mes cell', y = 'EMT score', fill = 'Mes cell')

# CNV
label <- fread('/data/chenrz/STARCH/s047801/labels_s047801.csv')
position <- fread('/data/chenrz/STARCH/s047801/tissue_positions_list.csv')
position$V1 <- paste0(sprintf("%.1f", position$array_row), 'x', sprintf("%.1f", position$array_col))
position <- left_join(position, label, by = 'V1')
position_f <- filter(position, V2 == 1) %>% group_by(array_col) %>% arrange(array_row) %>% slice_head(n = 3)
position$V2 <- ifelse(position$barcode %in% position_f$barcode, 3, position$V2)

S047801_16$CNV <- position$V2[match(colnames(S047801_16), position$barcode)]
S047801_16_mes$CNV <- position$V2[match(colnames(S047801_16_mes), position$barcode)]


##### 受体表达分布 #####
p1 <- SpatialFeaturePlot(S047801_16, features = 'IGF1R', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('IGF1R') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p2 <- SpatialFeaturePlot(S047801_16, features = 'PLXNA1', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('PLXNA1') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p3 <- SpatialFeaturePlot(S047801_16, features = 'ITGA3', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGA3') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p4 <- SpatialFeaturePlot(S047801_16, features = 'ITGB4', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGB4') +
  theme(legend.position = "right") +
  labs(fill = "SCT")

p1 + p2 + p3 + p4 + plot_layout(ncol = 4, guides = 'collect')

# ITGA3和ITGB4在Mes spot中的表达分布
p5 <- SpatialFeaturePlot(S047801_16_mes, features = 'ITGA3', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1)) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGA3') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p6 <- SpatialFeaturePlot(S047801_16_mes, features = 'ITGB4', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1)) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGB4') +
  theme(legend.position = "right") +
  labs(fill = "SCT")

p5 + p6 + plot_layout(ncol = 2, guides = 'collect')

tmp <- as.data.frame(S047801_16_mes@images$slice1$centroids@coords)
rownames(tmp) <- S047801_16_mes@images$slice1$centroids@cells
tmp <- left_join(data.frame(tmp, ID = rownames(tmp)), 
                 data.frame(S047801_16_mes@meta.data, ID = rownames(S047801_16_mes@meta.data)), 
                 by = 'ID')
rownames(tmp) <- tmp$ID
tmp$ITGA3 <- S047801_16_mes@assays$SCT@scale.data['ITGA3',]
tmp$ITGB4 <- S047801_16_mes@assays$SCT@scale.data['ITGB4',]
tmp$ITGA3_high <- ifelse(tmp$ITGA3 > quantile(tmp$ITGA3, 0.95), TRUE, FALSE)
tmp$ITGB4_high <- ifelse(tmp$ITGB4 > quantile(tmp$ITGB4, 0.95), TRUE, FALSE)
tmp$high <- ifelse(tmp$ITGA3_high == TRUE, 'ITGA3',
                   ifelse(tmp$ITGB4_high == TRUE, 'ITGB4', 'Low'))

p7 <- ggplot(filter(tmp, high != 'Low'), aes(x = y/100, fill = high)) + 
  geom_density(alpha = 0.5) +
  theme_bw() + 
  labs(fill = 'Hotspot', x = 'X axis', y = 'Probability')
p8 <- ggplot(filter(tmp, high != 'Low'), aes(x = x/100, fill = high)) + 
  geom_density(alpha = 0.5) +
  theme_bw() + 
  labs(fill = 'Hotspot', x = 'Y axis', y = 'Probability')

p7 / p8 + plot_layout(guides = "collect")


##### 划分EMT_early类型 #####
tmp_f <- filter(tmp, type == 'EMT_early')

mclust1 <- Mclust(tmp_f$ITGA3, G = 2)
mclust2 <- Mclust(tmp_f$ITGB4, G = 2)
tmp_f$type1 <- as.character(mclust1[["classification"]])
tmp_f$type2 <- as.character(mclust2[["classification"]])
tmp_f$type3 <- ifelse(tmp_f$type1 == '1' & tmp_f$type2 == '1', '1',
                      ifelse(tmp_f$type1 == '1' & tmp_f$type2 == '2', '2',
                             ifelse(tmp_f$type1 == '2' & tmp_f$type2 == '1', '3',
                                    ifelse(tmp_f$type1 == '2' & tmp_f$type2 == '2', '4', NA))))

tmp_f1 <- stack(tmp_f[,c('ITGA3', 'ITGB4')])
tmp_f1$type <- rep(tmp_f$type3,)
group_by(tmp_f1, type, ind) %>% summarise(mean = mean(values))
tmp_f1$type <- ifelse(tmp_f1$type == '1', 'ITGA3-ITGB4-',
                      ifelse(tmp_f1$type == '2', 'ITGA3-ITGB4+', 
                             ifelse(tmp_f1$type == '3', 'ITGA3+ITGB4-',
                                    ifelse(tmp_f1$type == '4', 'ITGA3+ITGB4+', NA))))
p8 <- ggplot(tmp_f1, aes(x = type, y = values, fill = ind, colour = ind)) + 
  geom_violin() +
  theme_bw() +
  labs(fill = 'Gene', colour = 'Gene', x = 'Type', y = 'Expression level')

S047801_16_mes_early <- subset(S047801_16_mes, type == 'EMT_early' & CNV %in% c(0,1))
S047801_16_mes_early$type1 <- tmp_f$type3[match(colnames(S047801_16_mes_early),tmp_f$ID)]
S047801_16_mes_early$type2 <- ifelse(S047801_16_mes_early$type1 == '1', 'ITGA3-ITGB4-',
                                     ifelse(S047801_16_mes_early$type1 == '2', 'ITGA3-ITGB4+', 
                                            ifelse(S047801_16_mes_early$type1 == '3', 'ITGA3+ITGB4-',
                                                   ifelse(S047801_16_mes_early$type1 == '4', 'ITGA3+ITGB4+', NA))))
table(S047801_16_mes_early$type2)

p9 <- ggplot(S047801_16_mes_early@meta.data, aes(x = type2, y = EMT, fill = type2)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "Type", y = 'EMT score', fill = 'Type') +
  scale_fill_manual(values = c('ITGA3-ITGB4-' = 'gray',
                               'ITGA3-ITGB4+' = '#66C2A5',
                               'ITGA3+ITGB4-' = '#FC8D62',
                               'ITGA3+ITGB4+' = '#E78AC3'))

p10 <- SpatialDimPlot(S047801_16_mes_early, group.by = 'type2', pt.size.factor = 5, image.alpha = 0, 
               cols = c('ITGA3-ITGB4-' = 'gray',
                        'ITGA3-ITGB4+' = '#66C2A5',
                        'ITGA3+ITGB4-' = '#FC8D62',
                        'ITGA3+ITGB4+' = '#E78AC3')) +
  theme_bw() +
  labs(fill = 'Type', x = '', y = '')

p9 + p10 + plot_layout(ncol = 2, width = c(2,2))
p8 + p9

####
## EMT_early中存在ITGA3/ITGB4的表达分布互斥倾向
####


##### 受配体细胞共定位测试 #####
###### ITGB4+EMT_early自分泌促进EMT ######
# 区分spot类型
coords <- as.data.frame(S047801_16@images$slice1.016um@boundaries$centroids@coords)
coords$ID <- S047801_16@images[["slice1.016um"]]@boundaries[["centroids"]]@cells
nnmatrix <- RANN::nn2(coords[,1:2], k = 120)$nn.idx
coords <- left_join(coords,
                    data.frame(S047801_16@meta.data, ID = rownames(S047801_16@meta.data)),
                    by = 'ID')

hotspot <- colnames(S047801_16_mes_early)[S047801_16_mes_early$type2 == 'ITGA3-ITGB4+'] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

S047801_16$spot <- ifelse(colnames(S047801_16) %in% hotspot, 'hotspot',
                          ifelse(colnames(S047801_16) %in% neighbor1, 'neighbor1', 
                                 ifelse(colnames(S047801_16) %in% neighbor2, 'neighbor2', 
                                        ifelse(colnames(S047801_16) %in% neighbor3, 'neighbor3',
                                               ifelse(colnames(S047801_16) %in% neighbor4, 'neighbor4',
                                                      ifelse(colnames(S047801_16) %in% neighbor5, 'neighbor5', 'else'))))))

# spot类型之间配体表达量对比
tmp_f_subset <- S047801_16@meta.data
tmp_f_subset$LAMA3 <- norm(S047801_16@assays$SCT@scale.data['LAMA3',])
tmp_f_subset$LAMB3 <- norm(S047801_16@assays$SCT@scale.data['LAMB3',])
tmp_f_subset$LAMC2 <- norm(S047801_16@assays$SCT@scale.data['LAMC2',])

p11 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMA3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p12 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMB3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p13 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMC2)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)

p11 + p12 + p13 + plot_layout(nrow = 1)

# spot类型之间EMT分值对比
S047801_16_mes_early$EMT <- S047801_16_mes_early$Mes1 - S047801_16_mes_early$Epi1
S047801_16_mes_early$spot <- ifelse(colnames(S047801_16_mes_early) %in% hotspot, 'hotspot',
                                    ifelse(colnames(S047801_16_mes_early) %in% neighbor1, 'neighbor1', 
                                           ifelse(colnames(S047801_16_mes_early) %in% neighbor2, 'neighbor2', 
                                                  ifelse(colnames(S047801_16_mes_early) %in% neighbor3, 'neighbor3',
                                                         ifelse(colnames(S047801_16_mes_early) %in% neighbor4, 'neighbor4',
                                                                ifelse(colnames(S047801_16_mes_early) %in% neighbor5, 'neighbor5', 'else'))))))
p14 <- SpatialDimPlot(S047801_16, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                    'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) + ggtitle('ITGB4+EMT_early')
p15 <- SpatialDimPlot(S047801_16_mes_early, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                              'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
p16 <- ggplot(S047801_16_mes_early@meta.data, aes(x = spot, y = EMT, fill = spot)) + 
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                               'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) +
  theme(axis.text.x = element_blank()) +
  labs(x = 'Spot', y = 'EMT score', fill = 'Spot')

p14 + p16
p15
p17 <- p16 + ggtitle('ITGA3-ITGB4+EMT_early')
p17

mean(filter(S047801_16_mes_early@meta.data, spot == 'hotspot')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor1')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor2')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor3')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor4')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor5')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'else')$EMT)


###### ITGA3+EMT_early自分泌促进MET ######
# 区分spot类型
hotspot <- colnames(S047801_16_mes_early)[S047801_16_mes_early$type2 == 'ITGA3+ITGB4-'] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

# spot类型分布可视化
S047801_16$spot <- ifelse(colnames(S047801_16) %in% hotspot, 'hotspot',
                          ifelse(colnames(S047801_16) %in% neighbor1, 'neighbor1', 
                                 ifelse(colnames(S047801_16) %in% neighbor2, 'neighbor2', 
                                        ifelse(colnames(S047801_16) %in% neighbor3, 'neighbor3',
                                               ifelse(colnames(S047801_16) %in% neighbor4, 'neighbor4',
                                                      ifelse(colnames(S047801_16) %in% neighbor5, 'neighbor5', 'else'))))))
SpatialDimPlot(S047801_16, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                              'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
# spot类型之间基因表达对比
tmp_f_subset <- S047801_16@meta.data
tmp_f_subset$LAMA3 <- norm(S047801_16@assays$SCT@scale.data['LAMA3',])
tmp_f_subset$LAMB3 <- norm(S047801_16@assays$SCT@scale.data['LAMB3',])
tmp_f_subset$LAMC2 <- norm(S047801_16@assays$SCT@scale.data['LAMC2',])

p18 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMA3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p19 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMB3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p20 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMC2)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)

p18 + p19 + p20 + plot_layout(nrow = 1)

# EMT_early spot之间EMT分值对比
S047801_16_mes_early$EMT <- S047801_16_mes_early$Mes1 - S047801_16_mes_early$Epi1
S047801_16_mes_early$spot <- ifelse(colnames(S047801_16_mes_early) %in% hotspot, 'hotspot',
                                    ifelse(colnames(S047801_16_mes_early) %in% neighbor1, 'neighbor1', 
                                           ifelse(colnames(S047801_16_mes_early) %in% neighbor2, 'neighbor2', 
                                                  ifelse(colnames(S047801_16_mes_early) %in% neighbor3, 'neighbor3',
                                                         ifelse(colnames(S047801_16_mes_early) %in% neighbor4, 'neighbor4',
                                                                ifelse(colnames(S047801_16_mes_early) %in% neighbor5, 'neighbor5', 'else'))))))
p21 <- SpatialDimPlot(S047801_16, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                    'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) + ggtitle('ITGA3+EMT_early')
p22 <- SpatialDimPlot(S047801_16_mes_early, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                              'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
p23 <- ggplot(S047801_16_mes_early@meta.data, aes(x = spot, y = EMT, fill = spot)) + 
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                               'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) +
  theme(axis.text.x = element_blank()) +
  labs(x = 'Spot', y = 'EMT score', fill = 'Spot')

p21 + p23
p22
p24 <- p23 + ggtitle('ITGA3+ITGB4-EMT_early')
p24

mean(filter(S047801_16_mes_early@meta.data, spot == 'hotspot')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor1')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor2')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor3')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor4')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'neighbor5')$EMT)
mean(filter(S047801_16_mes_early@meta.data, spot == 'else')$EMT)

# 合并箱线图
p8 + p9 + p17 + p24 + plot_layout(guides = 'collect')


###### ITGA3+EMT_early与CAF互作 ######
# 区分spot类型
hotspot <- colnames(S047801_16_mes_early)[S047801_16_mes_early$type2 == 'ITGA3+ITGB4-'] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

# spot类型分布可视化
S047801_16$spot <- ifelse(colnames(S047801_16) %in% hotspot, 'hotspot',
                          ifelse(colnames(S047801_16) %in% neighbor1, 'neighbor1', 
                                 ifelse(colnames(S047801_16) %in% neighbor2, 'neighbor2', 
                                        ifelse(colnames(S047801_16) %in% neighbor3, 'neighbor3',
                                               ifelse(colnames(S047801_16) %in% neighbor4, 'neighbor4',
                                                      ifelse(colnames(S047801_16) %in% neighbor5, 'neighbor5', 'else'))))))

# spot类型之间基因表达对比
tmp_f_subset <- S047801_16@meta.data
tmp_f_subset$FN1 <- norm(S047801_16@assays$SCT@scale.data['FN1',])
tmp_f_subset$PLAU <- norm(S047801_16@assays$SCT@scale.data['PLAU',])

ggplot(tmp_f_subset, aes(x = spot, y = FN1)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
ggplot(tmp_f_subset, aes(x = spot, y = PLAU)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)


# Steroseq ----------------------------------------------------------------
#### E1 ####
load('/data/chenrz/stseq_data/E1/E1/E1_20bin.RData')
load('/data/chenrz/stseq_data/E1/E1/E1_spe_20bin.RData')
load('/data/chenrz/stseq_data/E1/E1/RCTD_E1_20bin.RData')
load('/data/chenrz/stseq_data/E1/E1/RCTD_E1_mes_20bin.RData')

# QC
E1_20 <- E1_20[,setdiff(colnames(E1_20), colnames(E1_20_spe_qc)[E1_20_spe_qc$discard])]

# Deconv
E1_20 <- AddMetaData(E1_20, metadata = RCTD_E1_20@results$results_df)
E1_20 <- ScaleData(E1_20, features = rownames(E1_20))

E1_20_mes <- subset(E1_20, first_type == 'Mes' | second_type == 'Mes')
E1_20_mes <- ScaleData(E1_20_mes, features = rownames(E1_20_mes))
E1_20_mes <- AddMetaData(E1_20_mes, metadata = RCTD_E1_20_mes@results$results_df)
E1_20_mes$type <- as.character(E1_20_mes$first_type)
E1_20_mes$type <- ifelse(E1_20_mes$second_type == 'EMT_stable', 'EMT_stable', E1_20_mes$type) # 手动调整结果

E1_20_mes <- AddModuleScore(E1_20_mes, features = list(epi), name = 'Epi')
E1_20_mes <- AddModuleScore(E1_20_mes, features = list(mes), name = 'Mes')
E1_20_mes$EMT <- E1_20_mes$Mes1 - E1_20_mes$Epi1

SpatialDimPlot(E1_20_mes, group.by = 'type', pt.size.factor = 5, image.alpha = 0) + theme_bw()
ggplot(E1_20_mes@meta.data, aes(x = type, y = EMT, fill = type)) + 
  geom_boxplot() + 
  scale_x_discrete(limits = c("EMT_early", "EMT_stable", "EMT_diff")) +
  theme_bw() +
  coord_flip() +
  labs(x = 'Mes cell', y = 'EMT score', fill = 'Mes cell')

# CNV
label <- fread('/data/chenrz/STARCH/e1/labels_e1.csv')
position <- fread('/data/chenrz/STARCH/e1/tissue_positions_list.csv')
position$V7 <- paste0(sprintf("%.1f", position$V3), 'x', sprintf("%.1f", position$V4))
position <- left_join(position, label, by = c('V7' = 'V1'))

E1_20$CNV <- position$V2.y[match(colnames(E1_20), position$V1)]
E1_20_mes$CNV <- position$V2.y[match(colnames(E1_20_mes), position$V1)]


##### 受体表达分布 #####
p1 <- SpatialFeaturePlot(E1_20, features = 'IGF1R', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('IGF1R') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p2 <- SpatialFeaturePlot(E1_20, features = 'PLXNA1', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('PLXNA1') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p3 <- SpatialFeaturePlot(E1_20, features = 'ITGA3', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGA3') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p4 <- SpatialFeaturePlot(E1_20, features = 'ITGB4', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGB4') +
  theme(legend.position = "right") +
  labs(fill = "SCT")

p1 + p2 + p3 + p4 + plot_layout(ncol = 4, guides = 'collect')


# ITGA3和ITGB4在Mes spot中的表达分布
p5 <- SpatialFeaturePlot(E1_20_mes, features = 'ITGA3', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1)) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGA3') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p6 <- SpatialFeaturePlot(E1_20_mes, features = 'ITGB4', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1)) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGB4') +
  theme(legend.position = "right") +
  labs(fill = "SCT")

p5 + p6 + plot_layout(ncol = 2, guides = 'collect')

tmp <- as.data.frame(E1_20_mes@images$slice1$centroids@coords)
rownames(tmp) <- E1_20_mes@images$slice1$centroids@cells
tmp <- left_join(data.frame(tmp, ID = rownames(tmp)), 
                 data.frame(E1_20_mes@meta.data, ID = rownames(E1_20_mes@meta.data)), 
                 by = 'ID')
rownames(tmp) <- tmp$ID
tmp$ITGA3 <- E1_20_mes@assays$SCT@scale.data['ITGA3',]
tmp$ITGB4 <- E1_20_mes@assays$SCT@scale.data['ITGB4',]
tmp$ITGA3_high <- ifelse(tmp$ITGA3 > quantile(tmp$ITGA3, 0.95), TRUE, FALSE)
tmp$ITGB4_high <- ifelse(tmp$ITGB4 > quantile(tmp$ITGB4, 0.95), TRUE, FALSE)
tmp$high <- ifelse(tmp$ITGA3_high == TRUE, 'ITGA3',
                   ifelse(tmp$ITGB4_high == TRUE, 'ITGB4', 'Low'))

p7 <- ggplot(filter(tmp, high != 'Low'), aes(x = y/100, fill = high)) + 
  geom_density(alpha = 0.5) +
  theme_bw() + 
  labs(fill = 'Hotspot', x = 'X axis', y = 'Probability')
p8 <- ggplot(filter(tmp, high != 'Low'), aes(x = x/100, fill = high)) + 
  geom_density(alpha = 0.5) +
  theme_bw() + 
  labs(fill = 'Hotspot', x = 'Y axis', y = 'Probability')

p7 / p8 + plot_layout(guides = "collect")


##### 划分EMT_early类型 #####
tmp_f <- filter(tmp, type == 'EMT_early')

mclust1 <- Mclust(tmp_f$ITGA3, G = 2)
mclust2 <- Mclust(tmp_f$ITGB4, G = 2)
tmp_f$type1 <- as.character(mclust1[["classification"]])
tmp_f$type2 <- as.character(mclust2[["classification"]])
tmp_f$type3 <- ifelse(tmp_f$type1 == '1' & tmp_f$type2 == '1', '1',
                      ifelse(tmp_f$type1 == '1' & tmp_f$type2 == '2', '2',
                             ifelse(tmp_f$type1 == '2' & tmp_f$type2 == '1', '3',
                                    ifelse(tmp_f$type1 == '2' & tmp_f$type2 == '2', '4', NA))))

tmp_f1 <- stack(tmp_f[,c('ITGA3', 'ITGB4')])
tmp_f1$type <- rep(tmp_f$type3,)
group_by(tmp_f1, type, ind) %>% summarise(mean = mean(values))
tmp_f1$type <- ifelse(tmp_f1$type == '1', 'ITGA3-ITGB4-',
                      ifelse(tmp_f1$type == '2', 'ITGA3-ITGB4+', 
                             ifelse(tmp_f1$type == '3', 'ITGA3+ITGB4-',
                                    ifelse(tmp_f1$type == '4', 'ITGA3+ITGB4+', NA))))
p8 <- ggplot(tmp_f1, aes(x = type, y = values, fill = ind, colour = ind)) + 
  geom_violin() +
  theme_bw() +
  labs(fill = 'Gene', colour = 'Gene', x = 'Type', y = 'Expression level')

E1_20_mes_early <- subset(E1_20_mes, type == 'EMT_early' & CNV %in% c(0,1))
E1_20_mes_early$type1 <- tmp_f$type3[match(colnames(E1_20_mes_early),tmp_f$ID)]
E1_20_mes_early$type2 <- ifelse(E1_20_mes_early$type1 == '1', 'ITGA3-ITGB4-',
                                     ifelse(E1_20_mes_early$type1 == '2', 'ITGA3-ITGB4+', 
                                            ifelse(E1_20_mes_early$type1 == '3', 'ITGA3+ITGB4-',
                                                   ifelse(E1_20_mes_early$type1 == '4', 'ITGA3+ITGB4+', NA))))
table(E1_20_mes_early$type2)

p9 <- ggplot(E1_20_mes_early@meta.data, aes(x = type2, y = EMT, fill = type2)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "Type", y = 'EMT score', fill = 'Type') +
  scale_fill_manual(values = c('ITGA3-ITGB4-' = 'gray',
                               'ITGA3-ITGB4+' = '#66C2A5',
                               'ITGA3+ITGB4-' = '#FC8D62',
                               'ITGA3+ITGB4+' = '#E78AC3'))

p10 <- SpatialDimPlot(E1_20_mes_early, group.by = 'type2', pt.size.factor = 5, image.alpha = 0, 
                      cols = c('ITGA3-ITGB4-' = 'gray',
                               'ITGA3-ITGB4+' = '#66C2A5',
                               'ITGA3+ITGB4-' = '#FC8D62',
                               'ITGA3+ITGB4+' = '#E78AC3')) +
  theme_bw() +
  labs(fill = 'Type', x = '', y = '')

p9 + p10 + plot_layout(ncol = 2, width = c(2,2))
p8 + p9


##### 受配体细胞共定位测试 #####
###### ITGB4+EMT_early自分泌促进EMT ######
# 区分spot类型
coords <- as.data.frame(E1_20@images$slice1@boundaries$centroids@coords)
coords$ID <- E1_20@images[["slice1"]]@boundaries[["centroids"]]@cells
nnmatrix <- RANN::nn2(coords[,1:2], k = 120)$nn.idx
coords <- left_join(coords,
                    data.frame(E1_20@meta.data, ID = rownames(E1_20@meta.data)),
                    by = 'ID')

hotspot <- colnames(E1_20_mes_early)[E1_20_mes_early$type2 %in% c('ITGA3-ITGB4+', 'ITGA3+ITGB4+')] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

E1_20$spot <- ifelse(colnames(E1_20) %in% hotspot, 'hotspot',
                          ifelse(colnames(E1_20) %in% neighbor1, 'neighbor1', 
                                 ifelse(colnames(E1_20) %in% neighbor2, 'neighbor2', 
                                        ifelse(colnames(E1_20) %in% neighbor3, 'neighbor3',
                                               ifelse(colnames(E1_20) %in% neighbor4, 'neighbor4',
                                                      ifelse(colnames(E1_20) %in% neighbor5, 'neighbor5', 'else'))))))

# spot类型之间基因表达对比
tmp_f_subset <- E1_20@meta.data
tmp_f_subset$LAMA3 <- norm(E1_20@assays$SCT@scale.data['LAMA3',])
tmp_f_subset$LAMB3 <- norm(E1_20@assays$SCT@scale.data['LAMB3',])
tmp_f_subset$LAMC2 <- norm(E1_20@assays$SCT@scale.data['LAMC2',])

p11 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMA3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p12 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMB3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p13 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMC2)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)

p11 + p12 + p13 + plot_layout(nrow = 1)

# EMT_early spot之间EMT分值对比
E1_20_mes_early$EMT <- E1_20_mes_early$Mes1 - E1_20_mes_early$Epi1
E1_20_mes_early$spot <- ifelse(colnames(E1_20_mes_early) %in% hotspot, 'hotspot',
                                    ifelse(colnames(E1_20_mes_early) %in% neighbor1, 'neighbor1', 
                                           ifelse(colnames(E1_20_mes_early) %in% neighbor2, 'neighbor2', 
                                                  ifelse(colnames(E1_20_mes_early) %in% neighbor3, 'neighbor3',
                                                         ifelse(colnames(E1_20_mes_early) %in% neighbor4, 'neighbor4',
                                                                ifelse(colnames(E1_20_mes_early) %in% neighbor5, 'neighbor5', 'else'))))))
p14 <- SpatialDimPlot(E1_20, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                    'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) + ggtitle('ITGB4+EMT_early')
p15 <- SpatialDimPlot(E1_20_mes_early, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                              'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
p16 <- ggplot(E1_20_mes_early@meta.data, aes(x = spot, y = EMT, fill = spot)) + 
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                               'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) +
  theme(axis.text.x = element_blank()) +
  labs(x = 'Spot', y = 'EMT score', fill = 'Spot')

p14 + p16
p15
p17 <- p16 + ggtitle('ITGA3-ITGB4+EMT_early')
p17

mean(filter(E1_20_mes_early@meta.data, spot == 'hotspot')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor1')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor2')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor3')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor4')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor5')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'else')$EMT)


###### ITGA3+EMT_early自分泌促进MET ######
# 区分spot类型
hotspot <- colnames(E1_20_mes_early)[E1_20_mes_early$type2 %in% c('ITGA3+ITGB4-')] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

E1_20$spot <- ifelse(colnames(E1_20) %in% hotspot, 'hotspot',
                          ifelse(colnames(E1_20) %in% neighbor1, 'neighbor1', 
                                 ifelse(colnames(E1_20) %in% neighbor2, 'neighbor2', 
                                        ifelse(colnames(E1_20) %in% neighbor3, 'neighbor3',
                                               ifelse(colnames(E1_20) %in% neighbor4, 'neighbor4',
                                                      ifelse(colnames(E1_20) %in% neighbor5, 'neighbor5', 'else'))))))
# spot类型之间基因表达对比
tmp_f_subset <- E1_20@meta.data
tmp_f_subset$LAMA3 <- norm(E1_20@assays$SCT@scale.data['LAMA3',])
tmp_f_subset$LAMB3 <- norm(E1_20@assays$SCT@scale.data['LAMB3',])
tmp_f_subset$LAMC2 <- norm(E1_20@assays$SCT@scale.data['LAMC2',])

p18 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMA3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p19 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMB3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p20 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMC2)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)

p18 + p19 + p20 + plot_layout(nrow = 1)

# EMT_early spot之间EMT分值对比
E1_20_mes_early$EMT <- E1_20_mes_early$Mes1 - E1_20_mes_early$Epi1
E1_20_mes_early$spot <- ifelse(colnames(E1_20_mes_early) %in% hotspot, 'hotspot',
                                    ifelse(colnames(E1_20_mes_early) %in% neighbor1, 'neighbor1', 
                                           ifelse(colnames(E1_20_mes_early) %in% neighbor2, 'neighbor2', 
                                                  ifelse(colnames(E1_20_mes_early) %in% neighbor3, 'neighbor3',
                                                         ifelse(colnames(E1_20_mes_early) %in% neighbor4, 'neighbor4',
                                                                ifelse(colnames(E1_20_mes_early) %in% neighbor5, 'neighbor5', 'else'))))))
p21 <- SpatialDimPlot(E1_20, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                    'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) + ggtitle('ITGA3+EMT_early')
p22 <- SpatialDimPlot(E1_20_mes_early, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                              'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
p23 <- ggplot(E1_20_mes_early@meta.data, aes(x = spot, y = EMT, fill = spot)) + 
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                               'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) +
  theme(axis.text.x = element_blank()) +
  labs(x = 'Spot', y = 'EMT score', fill = 'Spot')

p21 + p23
p22
p24 <- p23 + ggtitle('ITGA3+ITGB4-EMT_early')
p24

mean(filter(E1_20_mes_early@meta.data, spot == 'hotspot')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor1')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor2')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor3')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor4')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'neighbor5')$EMT)
mean(filter(E1_20_mes_early@meta.data, spot == 'else')$EMT)

# 合并箱线图
p8 + p9 + p17 + p24 + plot_layout(guides = 'collect')


###### ITGA3+EMT_early与CAF互作 ######
# 区分spot类型
hotspot <- colnames(E1_20_mes_early)[E1_20_mes_early$type2 == 'ITGA3+'] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

# spot类型分布可视化
E1_20$spot <- ifelse(colnames(E1_20) %in% hotspot, 'hotspot',
                          ifelse(colnames(E1_20) %in% neighbor1, 'neighbor1', 
                                 ifelse(colnames(E1_20) %in% neighbor2, 'neighbor2', 
                                        ifelse(colnames(E1_20) %in% neighbor3, 'neighbor3',
                                               ifelse(colnames(E1_20) %in% neighbor4, 'neighbor4',
                                                      ifelse(colnames(E1_20) %in% neighbor5, 'neighbor5', 'else'))))))
SpatialDimPlot(E1_20, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                              'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
# spot类型之间基因表达对比
tmp_f_subset <- E1_20@meta.data
tmp_f_subset$FN1 <- norm(E1_20@assays$SCT@scale.data['FN1',])
tmp_f_subset$PLAU <- norm(E1_20@assays$SCT@scale.data['PLAU',])

ggplot(tmp_f_subset, aes(x = spot, y = FN1)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
ggplot(tmp_f_subset, aes(x = spot, y = PLAU)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)


#### E2 ####
load('/data/chenrz/stseq_data/E2/E2/E2_20bin.RData')
load('/data/chenrz/stseq_data/E2/E2/E2_spe_20bin.RData')
load('/data/chenrz/stseq_data/E2/E2/RCTD_E2_20bin.RData')
load('/data/chenrz/stseq_data/E2/E2/RCTD_E2_mes_20bin.RData')

# QC
E2_20 <- E2_20[,setdiff(colnames(E2_20), colnames(E2_20_spe_qc)[E2_20_spe_qc$discard])]

# Deconv
E2_20 <- AddMetaData(E2_20, metadata = RCTD_E2_20@results$results_df)
E2_20 <- ScaleData(E2_20, features = rownames(E2_20))

E2_20_mes <- subset(E2_20, first_type == 'Mes' | second_type == 'Mes')
E2_20_mes <- ScaleData(E2_20_mes, features = rownames(E2_20_mes))
E2_20_mes <- AddMetaData(E2_20_mes, metadata = RCTD_E2_20_mes@results$results_df)
E2_20_mes$type <- as.character(E2_20_mes$first_type)
E2_20_mes$type <- ifelse(E2_20_mes$second_type == 'EMT_stable', 'EMT_stable', E2_20_mes$type) # 手动调整结果

E2_20_mes <- AddModuleScore(E2_20_mes, features = list(epi), name = 'Epi')
E2_20_mes <- AddModuleScore(E2_20_mes, features = list(mes), name = 'Mes')
E2_20_mes$EMT1 <- E2_20_mes$Mes1 - E2_20_mes$Epi1

SpatialDimPlot(E2_20_mes, group.by = 'type', pt.size.factor = 5, image.alpha = 0) + theme_bw()
ggplot(E2_20_mes@meta.data, aes(x = type, y = EMT, fill = type)) + 
  geom_boxplot() + 
  scale_x_discrete(limits = c("EMT_early", "EMT_stable", "EMT_diff")) +
  theme_bw() +
  coord_flip() +
  labs(x = 'Mes cell', y = 'EMT score', fill = 'Mes cell')

# CNV
label <- fread('/data/chenrz/STARCH/e2/labels_e2.csv')
position <- fread('/data/chenrz/STARCH/e2/tissue_positions_list.csv')
position$V7 <- paste0(sprintf("%.1f", position$V3), 'x', sprintf("%.1f", position$V4))
position <- left_join(position, label, by = c('V7' = 'V1'))

E2_20$CNV <- position$V2.y[match(colnames(E2_20), position$V1)]
E2_20_mes$CNV <- position$V2.y[match(colnames(E2_20_mes), position$V1)]


##### 受体表达分布 #####
p1 <- SpatialFeaturePlot(E2_20, features = 'IGF1R', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('IGF1R') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p2 <- SpatialFeaturePlot(E2_20, features = 'PLXNA1', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('PLXNA1') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p3 <- SpatialFeaturePlot(E2_20, features = 'ITGA3', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGA3') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p4 <- SpatialFeaturePlot(E2_20, features = 'ITGB4', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1), min.cutoff = 0, max.cutoff = 1.5) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGB4') +
  theme(legend.position = "right") +
  labs(fill = "SCT")

p1 + p2 + p3 + p4 + plot_layout(ncol = 4, guides = 'collect')

# ITGA3和ITGB4在Mes spot中的表达分布
p5 <- SpatialFeaturePlot(E2_20_mes, features = 'ITGA3', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1)) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGA3') +
  theme(legend.position = "right") +
  labs(fill = "SCT")
p6 <- SpatialFeaturePlot(E2_20_mes, features = 'ITGB4', pt.size.factor = 5, stroke = .2, image.alpha = 0, alpha = c(0.1, 1)) +
  scale_fill_gradient(low = "white", high = "red", breaks = c(0, 1), labels = c("Low", "High")) +
  ggtitle('ITGB4') +
  theme(legend.position = "right") +
  labs(fill = "SCT")

p5 + p6 + plot_layout(ncol = 2, guides = 'collect')

tmp <- as.data.frame(E2_20_mes@images$slice1$centroids@coords)
rownames(tmp) <- E2_20_mes@images$slice1$centroids@cells
tmp <- left_join(data.frame(tmp, ID = rownames(tmp)), 
                 data.frame(E2_20_mes@meta.data, ID = rownames(E2_20_mes@meta.data)), 
                 by = 'ID')
rownames(tmp) <- tmp$ID
tmp$ITGA3 <- E2_20_mes@assays$SCT@scale.data['ITGA3',]
tmp$ITGB4 <- E2_20_mes@assays$SCT@scale.data['ITGB4',]
tmp$ITGA3_high <- ifelse(tmp$ITGA3 > quantile(tmp$ITGA3, 0.95), TRUE, FALSE)
tmp$ITGB4_high <- ifelse(tmp$ITGB4 > quantile(tmp$ITGB4, 0.95), TRUE, FALSE)
tmp$high <- ifelse(tmp$ITGA3_high == TRUE, 'ITGA3',
                   ifelse(tmp$ITGB4_high == TRUE, 'ITGB4', 'Low'))

p7 <- ggplot(filter(tmp, high != 'Low'), aes(x = y/100, fill = high)) + 
  geom_density(alpha = 0.5) +
  theme_bw() + 
  labs(fill = 'Hotspot', x = 'X axis', y = 'Probability')
p8 <- ggplot(filter(tmp, high != 'Low'), aes(x = x/100, fill = high)) + 
  geom_density(alpha = 0.5) +
  theme_bw() + 
  labs(fill = 'Hotspot', x = 'Y axis', y = 'Probability')

p7 / p8 + plot_layout(guides = "collect")


##### 划分EMT_early类型 #####
tmp_f <- filter(tmp, type == 'EMT_early')

mclust1 <- Mclust(tmp_f$ITGA3, G = 2)
mclust2 <- Mclust(tmp_f$ITGB4, G = 2)
tmp_f$type1 <- as.character(mclust1[["classification"]])
tmp_f$type2 <- as.character(mclust2[["classification"]])
tmp_f$type3 <- ifelse(tmp_f$type1 == '1' & tmp_f$type2 == '1', '1',
                      ifelse(tmp_f$type1 == '1' & tmp_f$type2 == '2', '2',
                             ifelse(tmp_f$type1 == '2' & tmp_f$type2 == '1', '3',
                                    ifelse(tmp_f$type1 == '2' & tmp_f$type2 == '2', '4', NA))))

tmp_f1 <- stack(tmp_f[,c('ITGA3', 'ITGB4')])
tmp_f1$type <- rep(tmp_f$type3,)
group_by(tmp_f1, type, ind) %>% summarise(mean = mean(values))
tmp_f1$type <- ifelse(tmp_f1$type == '1', 'ITGA3-ITGB4-',
                      ifelse(tmp_f1$type == '2', 'ITGA3-ITGB4+', 
                             ifelse(tmp_f1$type == '3', 'ITGA3+ITGB4-',
                                    ifelse(tmp_f1$type == '4', 'ITGA3+ITGB4+', NA))))
p8 <- ggplot(tmp_f1, aes(x = type, y = values, fill = ind, colour = ind)) + 
  geom_violin() +
  theme_bw() +
  labs(fill = 'Gene', colour = 'Gene', x = 'Type', y = 'Expression level')

E2_20_mes_early <- subset(E2_20_mes, type == 'EMT_early' & CNV %in% c(0,1,2))
E2_20_mes_early$type1 <- tmp_f$type3[match(colnames(E2_20_mes_early),tmp_f$ID)]
E2_20_mes_early$type2 <- ifelse(E2_20_mes_early$type1 == '1', 'ITGA3-ITGB4-',
                                ifelse(E2_20_mes_early$type1 == '2', 'ITGA3-ITGB4+', 
                                       ifelse(E2_20_mes_early$type1 == '3', 'ITGA3+ITGB4-',
                                              ifelse(E2_20_mes_early$type1 == '4', 'ITGA3+ITGB4+', NA))))
table(E2_20_mes_early$type2)

p9 <- ggplot(E2_20_mes_early@meta.data, aes(x = type2, y = EMT, fill = type2)) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "Type", y = 'EMT score', fill = 'Type') +
  scale_fill_manual(values = c('ITGA3-ITGB4-' = 'gray',
                               'ITGA3-ITGB4+' = '#66C2A5',
                               'ITGA3+ITGB4-' = '#FC8D62',
                               'ITGA3+ITGB4+' = '#E78AC3'))

p10 <- SpatialDimPlot(E2_20_mes_early, group.by = 'type2', pt.size.factor = 5, image.alpha = 0, 
                      cols = c('ITGA3-ITGB4-' = 'gray',
                               'ITGA3-ITGB4+' = '#66C2A5',
                               'ITGA3+ITGB4-' = '#FC8D62',
                               'ITGA3+ITGB4+' = '#E78AC3')) +
  theme_bw() +
  labs(fill = 'Type', x = '', y = '')

p8 + p9
p9 + p10 + plot_layout(ncol = 2, width = c(2,2))


##### 受配体细胞共定位测试 #####
###### ITGB4+EMT_early自分泌促进EMT ######
# 区分spot类型
coords <- as.data.frame(E2_20@images$slice1@boundaries$centroids@coords)
coords$ID <- E2_20@images[["slice1"]]@boundaries[["centroids"]]@cells
nnmatrix <- RANN::nn2(coords[,1:2], k = 120)$nn.idx
coords <- left_join(coords,
                    data.frame(E2_20@meta.data, ID = rownames(E2_20@meta.data)),
                    by = 'ID')

hotspot <- colnames(E2_20_mes_early)[E2_20_mes_early$type2 %in% c('ITGA3-ITGB4+', 'ITGA3+ITGB4+')] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

E2_20$spot <- ifelse(colnames(E2_20) %in% hotspot, 'hotspot',
                     ifelse(colnames(E2_20) %in% neighbor1, 'neighbor1', 
                            ifelse(colnames(E2_20) %in% neighbor2, 'neighbor2', 
                                   ifelse(colnames(E2_20) %in% neighbor3, 'neighbor3',
                                          ifelse(colnames(E2_20) %in% neighbor4, 'neighbor4',
                                                 ifelse(colnames(E2_20) %in% neighbor5, 'neighbor5', 'else'))))))

# spot类型之间基因表达对比
tmp_f_subset <- E2_20@meta.data
tmp_f_subset$LAMA3 <- norm(E2_20@assays$SCT@scale.data['LAMA3',])
tmp_f_subset$LAMB3 <- norm(E2_20@assays$SCT@scale.data['LAMB3',])
tmp_f_subset$LAMC2 <- norm(E2_20@assays$SCT@scale.data['LAMC2',])

p11 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMA3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p12 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMB3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p13 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMC2)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)

p11 + p12 + p13 + plot_layout(nrow = 1)

# EMT_early spot之间EMT分值对比
E2_20_mes_early$EMT <- E2_20_mes_early$Mes1 - E2_20_mes_early$Epi1
E2_20_mes_early$spot <- ifelse(colnames(E2_20_mes_early) %in% hotspot, 'hotspot',
                               ifelse(colnames(E2_20_mes_early) %in% neighbor1, 'neighbor1', 
                                      ifelse(colnames(E2_20_mes_early) %in% neighbor2, 'neighbor2', 
                                             ifelse(colnames(E2_20_mes_early) %in% neighbor3, 'neighbor3',
                                                    ifelse(colnames(E2_20_mes_early) %in% neighbor4, 'neighbor4',
                                                           ifelse(colnames(E2_20_mes_early) %in% neighbor5, 'neighbor5', 'else'))))))
p14 <- SpatialDimPlot(E2_20, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) + ggtitle('ITGB4+EMT_early')
p15 <- SpatialDimPlot(E2_20_mes_early, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                          'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
p16 <- ggplot(E2_20_mes_early@meta.data, aes(x = spot, y = EMT, fill = spot)) + 
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                               'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) +
  theme(axis.text.x = element_blank()) +
  labs(x = 'Spot', y = 'EMT score', fill = 'Spot')

p14 + p16
p15
p17 <- p16 + ggtitle('ITGA3-ITGB4+EMT_early')
p17

mean(filter(E2_20_mes_early@meta.data, spot == 'hotspot')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor1')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor2')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor3')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor4')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor5')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'else')$EMT)


###### ITGA3+EMT_early自分泌促进MET ######
# 区分spot类型
hotspot <- colnames(E2_20_mes_early)[E2_20_mes_early$type2 %in% c('ITGA3+ITGB4-')] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

E2_20$spot <- ifelse(colnames(E2_20) %in% hotspot, 'hotspot',
                     ifelse(colnames(E2_20) %in% neighbor1, 'neighbor1', 
                            ifelse(colnames(E2_20) %in% neighbor2, 'neighbor2', 
                                   ifelse(colnames(E2_20) %in% neighbor3, 'neighbor3',
                                          ifelse(colnames(E2_20) %in% neighbor4, 'neighbor4',
                                                 ifelse(colnames(E2_20) %in% neighbor5, 'neighbor5', 'else'))))))
# spot类型之间基因表达对比
tmp_f_subset <- E2_20@meta.data
tmp_f_subset$LAMA3 <- norm(E2_20@assays$SCT@scale.data['LAMA3',])
tmp_f_subset$LAMB3 <- norm(E2_20@assays$SCT@scale.data['LAMB3',])
tmp_f_subset$LAMC2 <- norm(E2_20@assays$SCT@scale.data['LAMC2',])

p18 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMA3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p19 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMB3)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
p20 <- ggplot(tmp_f_subset, aes(x = spot, y = LAMC2)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)

p18 + p19 + p20 + plot_layout(nrow = 1)

# EMT_early spot之间EMT分值对比
E2_20_mes_early$EMT <- E2_20_mes_early$Mes1 - E2_20_mes_early$Epi1
E2_20_mes_early$spot <- ifelse(colnames(E2_20_mes_early) %in% hotspot, 'hotspot',
                               ifelse(colnames(E2_20_mes_early) %in% neighbor1, 'neighbor1', 
                                      ifelse(colnames(E2_20_mes_early) %in% neighbor2, 'neighbor2', 
                                             ifelse(colnames(E2_20_mes_early) %in% neighbor3, 'neighbor3',
                                                    ifelse(colnames(E2_20_mes_early) %in% neighbor4, 'neighbor4',
                                                           ifelse(colnames(E2_20_mes_early) %in% neighbor5, 'neighbor5', 'else'))))))
p21 <- SpatialDimPlot(E2_20, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) + ggtitle('ITGA3+EMT_early')
p22 <- SpatialDimPlot(E2_20_mes_early, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                                          'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
p23 <- ggplot(E2_20_mes_early@meta.data, aes(x = spot, y = EMT, fill = spot)) + 
  geom_boxplot(alpha = 0.5) +
  theme_bw() +
  scale_fill_manual(values = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                               'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray')) +
  theme(axis.text.x = element_blank()) +
  labs(x = 'Spot', y = 'EMT score', fill = 'Spot')

p21 + p23
p22
p24 <- p23 + ggtitle('ITGA3+ITGB4-EMT_early')
p24

mean(filter(E2_20_mes_early@meta.data, spot == 'hotspot')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor1')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor2')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor3')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor4')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'neighbor5')$EMT)
mean(filter(E2_20_mes_early@meta.data, spot == 'else')$EMT)

# 合并箱线图
p8 + p9 + p17 + p24 + plot_layout(guides = 'collect')


###### ITGA3+EMT_early与CAF互作 ######
# 区分spot类型
hotspot <- colnames(E2_20_mes_early)[E2_20_mes_early$type2 == 'ITGA3+'] # hot spot
neighbor1 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 2:8]), 'ID']) 
neighbor2 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 9:24]), 'ID'])
neighbor3 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 25:48]), 'ID'])
neighbor4 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 49:80]), 'ID'])
neighbor5 <- unique(coords[as.numeric(nnmatrix[which(coords$ID %in% hotspot), 81:120]), 'ID']) # neighbor spot

# spot类型分布可视化
E2_20$spot <- ifelse(colnames(E2_20) %in% hotspot, 'hotspot',
                     ifelse(colnames(E2_20) %in% neighbor1, 'neighbor1', 
                            ifelse(colnames(E2_20) %in% neighbor2, 'neighbor2', 
                                   ifelse(colnames(E2_20) %in% neighbor3, 'neighbor3',
                                          ifelse(colnames(E2_20) %in% neighbor4, 'neighbor4',
                                                 ifelse(colnames(E2_20) %in% neighbor5, 'neighbor5', 'else'))))))
SpatialDimPlot(E2_20, group.by = 'spot', image.alpha = 0, pt.size.factor = 3.8, cols = c('hotspot' = 'red','neighbor1' = 'blue', 'neighbor2' = 'green', 
                                                                                         'neighbor3' = 'yellow', 'neighbor4' = 'orange', 'neighbor5' = 'purple', 'else' = 'gray'))
# spot类型之间基因表达对比
tmp_f_subset <- E2_20@meta.data
tmp_f_subset$FN1 <- norm(E2_20@assays$SCT@scale.data['FN1',])
tmp_f_subset$PLAU <- norm(E2_20@assays$SCT@scale.data['PLAU',])

ggplot(tmp_f_subset, aes(x = spot, y = FN1)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)
ggplot(tmp_f_subset, aes(x = spot, y = PLAU)) + 
  geom_violin() +
  stat_summary(
    fun = mean,         
    geom = "crossbar",
    width = 0.3,
    color = "red",
    size = 0.5)