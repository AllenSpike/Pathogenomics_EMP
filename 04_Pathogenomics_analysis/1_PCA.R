## WSI病理特征PCA
## 20240926

library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)
library(maftools)
library(estimate)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(patchwork)
library(survival)
library(survminer)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(multiOmicsViz)
library(DESeq2)
library(xCell)
library(ComplexHeatmap)
library(enrichplot)
library(export)

agg_type = c('mean', 'sd', 'skewness', 'kurtosis') # 聚合类型
path_output = './data/patch/xj_4096/' 
path_output2 = './data/patch/tcga_4096/'
path_output3 = './data/patch/gdph_supp_4096/'

#### XJ特征矩阵 ####
feat_xj <- t(do.call(rbind, lapply(agg_type, function(x){
  s <- read.xlsx(paste0(path_output,"/feat_agg.xlsx"), rowNames = T, sheet = x)[1:252,] # 选择前252个特征
  rownames(s) <- paste0(rownames(s), "_", x)
  return(s)
}
)))
colnames(feat_xj) <- gsub('GraphAvg_', 'GraphAvg\\.', colnames(feat_xj)) 
feat_xj <- apply(feat_xj, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

# 标准化
feat_xj_norm <- scale(feat_xj)
feat_xj_norm <- apply(feat_xj_norm, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

#### TCGA特征矩阵 ####
feat_tcga <- t(do.call(rbind, lapply(agg_type, function(x){
  s <- read.xlsx(paste0(path_output2,"/feat_agg.xlsx"), rowNames = T, sheet = x)[1:252,] # 选择前522个特征
  rownames(s) <- paste0(rownames(s), "_", x)
  return(s)
}
)))
rownames(feat_tcga) <- substr(rownames(feat_tcga), 1, 12)
colnames(feat_tcga) <- gsub('GraphAvg_', 'GraphAvg\\.', colnames(feat_tcga))
feat_tcga <- apply(feat_tcga, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

# 标准化
feat_tcga_norm <- scale(feat_tcga)
feat_tcga_norm <- apply(feat_tcga_norm, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na


#### GDPH特征矩阵 ####
feat_gdph <- t(do.call(rbind, lapply(agg_type, function(x){
  s <- read.xlsx(paste0(path_output3,"/feat_agg.xlsx"), rowNames = T, sheet = x)[1:252,] # 选择前252个特征
  rownames(s) <- paste0(rownames(s), "_", x)
  return(s)
}
)))
rownames(feat_gdph) <- substr(rownames(feat_gdph), 1, 7)
colnames(feat_gdph) <- gsub('GraphAvg_', 'GraphAvg\\.', colnames(feat_gdph))
feat_gdph <- apply(feat_gdph, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

# 标准化
feat_gdph_norm <- scale(feat_gdph)
feat_gdph_norm <- apply(feat_gdph_norm, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

#### PCA ####
feat <- rbind(feat_xj_norm, feat_tcga_norm, feat_gdph_norm)
pca <- prcomp(feat, scale = F, center = F)
pca_df <- as.data.frame(pca$x)
summ <- summary(pca)
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

pdf('./fig1.pca.pdf', width = 6, height = 5)
ggplot(data = pca_df, aes(x = PC1,y = PC2, color = c(rep('XJ', 126), rep('TCGA', 51), rep('GDPH', 48))))+
  geom_point(size = 1.5)+
  labs(x = xlab, y = ylab , color = "Group", title = "4096_standardize")+
  guides(fill = "none") +
  stat_ellipse(level = 0.9) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
dev.off()

rm(list = ls())
gc()
