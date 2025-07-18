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

#### 靶点生存依赖验证 ####
model <- read.csv('/data/chenrz/resource/Model.csv')
screen <- read.csv('/data/chenrz/resource/ScreenSequenceMap.csv')
effect <- fread('/data/chenrz/resource/ScreenGeneEffect.csv')

model_f <- filter(model, DepmapModelType == 'ESCC')
screen_f <- filter(screen, ModelID %in% model_f$ModelID)
effect_f <- filter(effect, V1 %in% screen_f$ScreenID)
colnames(effect_f) <- gsub(' \\(.*', '', colnames(effect_f))

gene_vul <- effect_f[, c('V1', 'CACNA1C', 'IGF1R', 'PLXNA1', 'ITGA3', 'ITGB4')]
cellname <- model_f[match(screen_f$ModelID[match(gene_vul$V1, screen_f$ScreenID)],model_f$ModelID), 'StrippedCellLineName']
gene_vul$V1 <- cellname
gene_vul <- aggregate(. ~ V1, data = gene_vul, FUN = median, na.rm = TRUE)

mtp_potential_f$X1 <- gsub('_.*', '',mtp_potential_f$X1)
mtp_potential_f <- left_join(mtp_potential_f, gene_vul, by = c('X1' = 'V1'))

ggplot(mtp_potential_f, aes(x = State, y = CACNA1C, fill = State)) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width = 0.2, linewidth= 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 3) +
  ggpubr::geom_signif(
    comparisons = list(
      c(1, 2),
      c(2, 3),
      c(1, 3)
    ),
    test = "wilcox.test",
    step_increase = 0.1) +
  labs(x = 'EMT state', y = 'CACNA1C dependency', fill = 'EMT state') +
  scale_fill_manual(values = c('EMT_early' = "#b6dd49", 'EMT_stable' = "#f9e762", 'EMT_diff' = "#74cc64")) + 
  theme_bw()

####
## CACNA1C是EMT中后期细胞的特异性靶点，敲除该基因可以降低细胞周期、提升角质化，促使肿瘤细胞走向终末分化
####

#### 受配体/靶点功能验证 ####
info_sig <- fread('/data/chenrz/resource/GSE92742_Broad_LINCS_sig_info.txt') %>% filter(pert_type %in% c('trt_sh', 'trt_xpr'))
info_gene <- fread('/data/chenrz/resource/GSE92742_Broad_LINCS_gene_info.txt')
load('/data/chenrz/resource/lincs_knockdown_escc_pr.RData')

gene_knock <- lincs[, c('CACNA1C', 'IGF1R', 'PLXNA1', 'ITGA3', 'ITGB4')]
epi_f <- intersect(epi, rownames(gene_knock))
mes_f <- intersect(mes, rownames(gene_knock))

# 敲除对上皮/间充质基因的影响
ht <- gene_knock[c(epi_f, mes_f),]
dm <- dist(ht[epi_f,], method = "euclidean")
hc_1 <- hclust(dm)
dm <- dist(ht[mes_f,], method = "euclidean")
hc_2 <- hclust(dm)
ord <- c(epi_f[hc_1$order], mes_f[hc_2$order])

Heatmap(ht[ord,], 
        row_title = 'Marker gene', 
        column_title = 'Knock-down gene',
        row_split = c(rep(1, 12), rep(2, 14)), 
        border = T,
        cluster_rows = F,
        show_column_dend = F,
        heatmap_legend_param = list(title = "LogFC"))

# 敲除对上皮/间充质表型的影响
anno_manual <- rbind(data.frame(term = 'Mesenchymal', gene = mes_f), 
                     data.frame(term = 'Epithelia', gene = epi_f))
res_GSEA_knock <- apply(gene_knock, 2, function(x){
  x <- sort(x, decreasing = T)
  res <- list(
    manual = GSEA(geneList = x, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_manual)
  )
  return(res)
})
lapply(res_GSEA_knock, function(x){x[["manual"]]@result})

# 敲除对其他通路的影响
anno_reactome <- read.gmt('/data/chenrz/resource/reactome.v2024.1.Hs.symbols.gmt')

res_GSEA_knock2 <- apply(gene_knock, 2, function(x){
  x <- sort(x, decreasing = T)
  res <- GSEA(geneList = x, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_reactome)
  return(res)
})

kf <- res_GSEA_knock2$CACNA1C
kf@result$Description <- str_to_title(gsub('_', ' ', gsub('REACTOME_', '', kf@result$Description)))

gseaNb(object = kf,
       curveCol = jjAnno::useMyCol('paired', 2),
       subPlot = 2,
       termWidth = 35,
       addGene = c('CSTA', 'CAPN1', 'KRT14', 'KNTC1', 'PSMD11', 'PSMA6'),
       geneSetID = c('REACTOME_KERATINIZATION', 
                     'REACTOME_CELL_CYCLE_CHECKPOINTS'))

# 同阶段受体之间的拮抗作用
cor(gene_knock, method = 'spearman')

# ITGA3 ITGB4
p1 <- ggplot(data = gene_knock[,c('ITGA3', 'ITGB4')], aes(x = ITGA3, y = ITGB4)) + 
  geom_pointdensity(aes(x = ITGA3,y = ITGB4), size = 0.1) +  
  scale_color_viridis() + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() + 
  theme(legend.position = "none") +
  ggpubr::stat_cor(method = "spearman")
p2 <- ggplot(data = gene_knock[epi_f, c('ITGA3', 'ITGB4')], aes(x = ITGA3, y = ITGB4)) + 
  geom_point(aes(x = ITGA3,y = ITGB4), color = "#d20000", size = 2) +  
  scale_color_viridis() + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  ggpubr::stat_cor(method = "spearman")
p3 <- ggplot(data = gene_knock[mes_f, c('ITGA3', 'ITGB4')], aes(x = ITGA3, y = ITGB4)) + 
  geom_point(aes(x = ITGA3,y = ITGB4), color = "#5abf22", size = 2) +  
  scale_color_viridis() + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  ggpubr::stat_cor(method = "spearman")

p1 + (p2/p3)

# PLXNA1 IGF1R
p1 <- ggplot(data = gene_knock[,c('PLXNA1', 'IGF1R')], aes(x = PLXNA1, y = IGF1R)) + 
  geom_pointdensity(aes(x = PLXNA1,y = IGF1R), size = 0.1) +  
  scale_color_viridis() + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() + 
  theme(legend.position = "none") +
  ggpubr::stat_cor(method = "spearman")
p2 <- ggplot(data = gene_knock[epi_f, c('PLXNA1', 'IGF1R')], aes(x = PLXNA1, y = IGF1R)) + 
  geom_point(aes(x = PLXNA1,y = IGF1R), color = "#d20000", size = 2) +  
  scale_color_viridis() + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  ggpubr::stat_cor(method = "spearman")
p3 <- ggplot(data = gene_knock[mes_f, c('PLXNA1', 'IGF1R')], aes(x = PLXNA1, y = IGF1R)) + 
  geom_point(aes(x = PLXNA1, y = IGF1R), color = "#5abf22", size = 2) +  
  scale_color_viridis() + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  theme_bw() +
  ggpubr::stat_cor(method = "spearman")

p1 + (p2/p3)