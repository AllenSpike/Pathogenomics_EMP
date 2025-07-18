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

#### 受配体细胞定位 ####
##### 受配体筛选 #####
# 风险值正相关 + 无监督高风险组上调 + 危害性的受配体
load('/data/chenrz/LRinter.dataframe.RData')

df1 <- filter(LRinter.dataframe, corr > 0)
df2 <- filter(LRinter.dataframe, logFC > 0)
df3 <- filter(hr, HR_pvalue < 0.05 & HR > 1)
lr <- str_split(Reduce(intersect, list(df1$LR, df2$LR, rownames(df3))), '_', simplify = T)

##### 细胞大类定位 #####
marker <- FindAllMarkers(scobj, logfc.threshold = 0.25, min.pct = 0.1) # 细胞大类marker

cell <- apply(lr, 1, function(x){
  unlist(lapply(x, function(y){
    res <- marker[marker$gene == y, ]
    res <- arrange(res, desc(avg_log2FC))
    cell <- res[1, 'cluster']
    lfc <- res[1, 'avg_log2FC']
    return(as.character(cell))
  }))
})
lfc <- apply(lr, 1, function(x){
  unlist(lapply(x, function(y){
    res <- marker[marker$gene == y, ]
    res <- arrange(res, desc(avg_log2FC))
    cell <- res[1, 'cluster']
    lfc <- res[1, 'avg_log2FC']
    return(as.numeric(lfc))
  }))
})

lr <- cbind(as.data.frame(lr), as.data.frame(t(cell)), as.data.frame(t(lfc)))
colnames(lr) <- c('L', 'R', 'LCell', 'RCell', 'L_lfc', 'R_lfc')
lr$LR <- paste(lr$L, lr$R, sep = '_')
lr$pathway <- LRinter.dataframe$pw.name[match(lr$LR, LRinter.dataframe$LR)]
lr <- lr[complete.cases(lr),]

sort(table(lr$pathway), decreasing = T)
sort(table(lr$L), decreasing = T)
sort(table(lr$R), decreasing = T)
sort(table(lr$LCell), decreasing = T)
sort(table(lr$RCell), decreasing = T)

## 小提琴图
plots <- apply(lr, 1, function(x){
  p1 <- VlnPlot(scobj, features = x[1], group.by = 'Cell', pt.size = 0)
  p2 <- VlnPlot(scobj, features = x[2], group.by = 'Cell', pt.size = 0)
  p <- p1 + p2 + plot_layout(ncol = 2, guides = 'collect')
  return(p)
})
names(plots) <- lr$LR

## 密度图
scobj$seurat_clusters <- scobj$Cell

ggscplot(object = scobj, features = c(lr$L[1], lr$R[1])) +
  geom_density2d_filled(bins = 10,show.legend = F) +
  geom_scPoint(aes(color = value),size = 0.1) +
  scale_color_gradient(low = "white",high = "#CC0033",name = "Expression") +
  facet_wrap(~gene_name,ncol = 3) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold.italic",size = rel(1)),
        axis.text = element_blank()) +
  scale_fill_manual(values = colorRampPalette(c("white","#336633"))(10))

# 细胞大类汇总表
ddf <- data.frame(Cell = unique(c(lr$LCell, lr$RCell)),
                  Receptor_count = as.numeric(table(lr$RCell)[unique(c(lr$LCell, lr$RCell))]),
                  Ligand_count = as.numeric(table(lr$LCell)[unique(c(lr$LCell, lr$RCell))]))
ddf <- arrange(ddf, desc(Receptor_count))

ggplot() + 
  annotate(geom = "table", 
           x = 1, y = 1, 
           label = list(ddf)) +
  theme_void()

Heatmap(as.matrix(table(lr$LCell, lr$RCell)), border = T, show_column_dend = F, show_row_dend = F,
        rect_gp = gpar(col = "white", lwd = 2),
        heatmap_legend_param = list(title = "Count"))

# 保存细胞大类汇总表
write.csv(lr, file = '/data/chenrz/result_pathbulk/3_lr_result_maincell.csv', row.names = F)


##### 细胞亚群定位 #####
# 细胞亚群marker
marker_epi <- FindAllMarkers(scobj_epi, group.by = 'Cell2')
marker_fibro <- FindAllMarkers(scobj_fibro, group.by = 'Cell2')
marker_endo <- FindAllMarkers(scobj_endo, group.by = 'Cell2')
marker_mye <- FindAllMarkers(scobj_mye, group.by = 'Cell2')
marker_tcell <- FindAllMarkers(scobj_tcell, group.by = 'Cell2')
marker_bcell <- FindAllMarkers(scobj_bcell, group.by = 'Cell2')

marker_ls <- list()
marker_ls[['Epithelia']] <- marker_epi
marker_ls[['Fibroblast']] <- marker_fibro
marker_ls[['Endothelial']] <- marker_endo
marker_ls[['Myeloid']] <- marker_mye
marker_ls[['Tcell']] <- marker_tcell
marker_ls[['Bcell']] <- marker_bcell
marker_ls[['FRC']] <- filter(marker, cluster == 'FRC')
marker_ls[['Pericytes']] <- filter(marker, cluster == 'Pericytes')

cell <- apply(lr, 1, function(x){
  
  L <- x[1]
  R <- x[2]
  Lcell <- x[3]
  Rcell <- x[4]
  Lmarker <- filter(marker_ls[[Lcell]], gene == L) %>% arrange(desc(avg_log2FC))
  Rmarker <- filter(marker_ls[[Rcell]], gene == R) %>% arrange(desc(avg_log2FC))
  
  L_cell <- as.character(Lmarker[1, 'cluster'])
  R_cell <- as.character(Rmarker[1, 'cluster'])
  
  return(c(L_cell, R_cell))
})

lfc <- apply(lr, 1, function(x){
  
  L <- x[1]
  R <- x[2]
  Lcell <- x[3]
  Rcell <- x[4]
  Lmarker <- filter(marker_ls[[Lcell]], gene == L) %>% arrange(desc(avg_log2FC))
  Rmarker <- filter(marker_ls[[Rcell]], gene == R) %>% arrange(desc(avg_log2FC))
  
  L_lfc <- as.numeric(Lmarker[1, 'avg_log2FC'])
  R_lfc <- as.numeric(Rmarker[1, 'avg_log2FC'])
  
  return(c(L_lfc, R_lfc))
  
})

lr <- cbind(lr, as.data.frame(t(cell)), as.data.frame(t(lfc)))
colnames(lr)[c(9, 10, 11, 12)] <- c('LCell_sub', 'RCell_sub', 'Llfc_sub', 'Rlfc_sub')
lr <- lr[complete.cases(lr),]

# 细胞亚群汇总表
ddf <- data.frame(Cell = unique(c(lr$LCell_sub, lr$RCell_sub)),
                  Receptor_count = as.numeric(table(lr$RCell_sub)[unique(c(lr$LCell_sub, lr$RCell_sub))]),
                  Ligand_count = as.numeric(table(lr$LCell_sub)[unique(c(lr$LCell_sub, lr$RCell_sub))]))
ddf <- arrange(ddf, desc(Receptor_count))
ggplot() + 
  annotate(geom = "table", 
           x = 1, y = 1, 
           label = list(ddf)) +
  theme_void()

Heatmap(as.matrix(table(lr$LCell_sub, lr$RCell_sub)), border = T, show_column_dend = F, show_row_dend = F,
        heatmap_legend_param = list(title = "Count"),
        col = circlize::colorRamp2(c(0, 6), c('white' ,'red')), 
        row_title = 'Sender', column_title = 'Receiver')

# 保存细胞亚群汇总表
write.csv(lr, file = '/data/chenrz/result_pathbulk/3_lr_result_subcell.csv', row.names = F)


##### Mes特异性受配体 #####
scobj_merge <- merge(scobj_bcell, y = c(scobj_tcell, scobj_endo, scobj_epi, scobj_fibro, scobj_frc, scobj_mye), 
                     project = "merged_project")
scobj$Cell2 <- scobj_merge$Cell2[match(colnames(scobj), colnames(scobj_merge))]

save(scobj, file = '/data/chenrz/scseq_gse160269/scobj.RData')
rm(scobj_merge)
gc()

# Mes特异性受体
pseudobulk <- AggregateExpression(scobj, group.by = 'Cell2')
pseudobulk <- as.matrix(pseudobulk[["RNA"]])

lr_mes <- filter(lr, RCell_sub == 'Mes')
table(lr_mes$R)

plots <- lapply(c(names(table(lr_mes$R)), 'CACNA1C'), function(x){
  p <- VlnPlot(scobj_epi, features = x, group.by = 'Cell2', pt.size = 0.001,
               col = colorRampPalette(brewer.pal(8, "Set2"))(8)) 
  
  return(p)
})
plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plot_layout(ncol = 5, guides = 'collect')

# Mes特异性配体
table(lr_mes$L)
table(lr_mes$LCell_sub)

# 受配体交互热图
ht <- pseudobulk[, unique(c(lr_mes$LCell_sub, scobj_epi$Cell2))]
ht <- t(scale(t(ht)))
ht_f <- ht[unique(c(lr_mes$L, lr_mes$R)),]

ht_f_r <- as.matrix(ht[unique(lr_mes$R), unique(scobj_epi$Cell2)])
ht_f_l <- as.matrix(ht[unique(lr_mes$L), unique(lr_mes$LCell_sub)])

ht1 <- Heatmap(ht_f_r, show_column_dend = F, show_row_dend = F, border = T,
               row_names_side = 'left',
               heatmap_legend_param = list(title = "Receptor Exp"))
ht2 <- Heatmap(ht_f_l, show_column_dend = F, show_row_dend = F, border = T,
               heatmap_legend_param = list(title = "Ligand Exp"))
draw(ht2, heatmap_legend_side = 'left')
draw(ht1)


#### 靶点细胞定位 ####
target <- read.csv('/data/chenrz/target_enriched_chembl_targets.csv')

##### 靶点筛选 #####
# scges > 0 & padj < 0.05
target <- filter(target, padj < 0.05 & score > 0)

##### 细胞大类定位 #####
cell <- unlist(lapply(target$target, function(x){
  res <- marker[marker$gene == x, ]
  res <- arrange(res, desc(avg_log2FC))
  cell <- res[1, 'cluster']
  lfc <- res[1, 'avg_log2FC']
  return(as.character(cell))
}))

lfc <- unlist(lapply(target$target, function(x){
  res <- marker[marker$gene == x, ]
  res <- arrange(res, desc(avg_log2FC))
  cell <- res[1, 'cluster']
  lfc <- res[1, 'avg_log2FC']
  return(as.numeric(lfc))
}))

target <- cbind(as.data.frame(target), cell, lfc)

# 保存细胞定位汇总表
write.csv(target, file = '/data/chenrz/result_pathbulk/3_target_result_subcell.csv', row.names = F)

##### 上皮细胞亚群定位 #####
cell_epi_sub <- unlist(lapply(target$target, function(x){
  res <- marker_epi[marker_epi$gene == x, ]
  res <- arrange(res, desc(avg_log2FC))
  cell <- res[1, 'cluster']
  lfc <- res[1, 'avg_log2FC']
  return(as.character(cell))
}))

lfc_epi_sub <- unlist(lapply(target$target, function(x){
  res <- marker_epi[marker_epi$gene == x, ]
  res <- arrange(res, desc(avg_log2FC))
  cell <- res[1, 'cluster']
  lfc <- res[1, 'avg_log2FC']
  return(as.numeric(lfc))
}))

target <- cbind(target, cell_epi_sub, lfc_epi_sub)
