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

##### EMT阶段干性能力对比 #####
cytotrace2_res <- cytotrace2(scobj_epi_mes, 
                             is_seurat = TRUE, 
                             slot_type = "counts", 
                             species = 'human')

annotation <- data.frame(phenotype = cytotrace2_res@meta.data$EMT_state) %>% 
  set_rownames(., colnames(cytotrace2_res))

plots <- plotData(cytotrace2_result = cytotrace2_res, 
                  annotation = annotation, 
                  is_seurat = TRUE)

p1 <- plots$CytoTRACE2_UMAP
p2 <- plots$CytoTRACE2_Potency_UMAP
p3 <- plots$CytoTRACE2_Relative_UMAP
p4 <- plots$Phenotype_UMAP
p5 <- plots$CytoTRACE2_Boxplot_byPheno

(p1 | p2) / (p3 | p5)

# 干性卡值
max(filter(cytotrace2_res@meta.data, CytoTRACE2_Potency == 'Differentiated')$CytoTRACE2_Score) # 0.165
max(filter(cytotrace2_res@meta.data, CytoTRACE2_Potency == 'Unipotent')$CytoTRACE2_Score) # 0.333
max(filter(cytotrace2_res@meta.data, CytoTRACE2_Potency == 'Oligopotent')$CytoTRACE2_Score) # 0.499
max(filter(cytotrace2_res@meta.data, CytoTRACE2_Potency == 'Multipotent')$CytoTRACE2_Score) # 0.662

scobj_epi_mes$potency <- cytotrace2_res$CytoTRACE2_Score
scobj_epi_mes$EMT_state <- factor(scobj_epi_mes$EMT_state, levels = c('EMT_early', 'EMT_stable', 'EMT_diff'))

p5 <- ggplot(scobj_epi_mes@meta.data, aes(x = EMT_state, y = potency, fill = EMT_state)) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width = 0.2, linewidth= 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 1.25, color = "white") +
  geom_hline(yintercept = 0.1652585, linetype = 'dashed') +
  geom_hline(yintercept = 0.333304, linetype = 'dashed') +
  geom_hline(yintercept = 0.4998168, linetype = 'dashed') +
  geom_hline(yintercept = 0.6618212, linetype = 'dashed') +
  ggpubr::geom_signif(
    comparisons = list(
      c(1, 2),
      c(2, 3),
      c(1, 3)
    ),
    test = "wilcox.test",
    step_increase = 0.1) +
  labs(x = 'EMT state', y = 'Diffrentiation potential', fill = 'EMT state') +
  scale_fill_manual(values = c('EMT_early' = "#b6dd49", 'EMT_stable' = "#f9e762", 'EMT_diff' = "#74cc64")) + 
  theme_bw()

##### EMT阶段增殖能力对比 ##### 
scobj_epi_mes <- CellCycleScoring(scobj_epi_mes, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
Idents(scobj_epi_mes) <- 'EMT_state'

scobj_epi_mes$Cycle.score <- rowMeans(scobj_epi_mes@meta.data[, c('S.Score', 'G2M.Score')])

ht <- scobj_epi_mes@meta.data
ht$EMT_state <- factor(ht$EMT_state, levels = c('EMT_early', 'EMT_stable', 'EMT_diff'))
p6 <- ggplot(ht, aes(x = EMT_state, y = Cycle.score, fill = EMT_state)) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width = 0.2, linewidth= 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 1.25, color = "white") +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ggpubr::geom_signif(
    comparisons = list(
      c(1, 2),
      c(2, 3),
      c(1, 3)
    ),
    test = "wilcox.test",
    step_increase = 0.1) +
  labs(x = 'EMT state', y = 'Cycle score', fill = 'EMT state') +
  scale_fill_manual(values = c('EMT_early' = "#b6dd49", 'EMT_stable' = "#f9e762", 'EMT_diff' = "#74cc64")) + 
  theme_bw()

##### EMT阶段转移能力对比 #####
ncbi <- fread('/data/chenrz/resource/Homo_sapiens.gene_info', check.names = F)
ncbi$gene_id <- str_extract_all(string = ncbi$dbXrefs, pattern = 'ENSG[0-9]+', simplify = T)[,1]

exp_ccle <- as.data.frame(fread('/data/chenrz/resource/CCLE_RNAseq_rsem_genes_tpm_20180929.txt.gz'))
exp_ccle$gene_id <- substr(exp_ccle$gene_id, 1, 15) 
exp_ccle <- left_join(exp_ccle, ncbi[,c(3, 17)], by = 'gene_id') %>%
  select('gene_id', 'transcript_ids', 'Symbol', everything())

mtp_annot <- read.xlsx('/data/chenrz/resource/MetMap cell line annotation.xlsx', sheet = 1)
mtp_potential <- do.call(rbind, lapply(1:5, function(i) read.xlsx('/data/chenrz/resource/MetMap 500 met potential.xlsx', sheet = i)))
mtp_potential$target_tissue <- c(rep('brain', 488), rep('lung', 488), rep('liver', 488), rep('bone', 488), rep('kidney', 488))
mtp_potential$group <- ifelse(mtp_potential$mean <= -4, 'non-metastatic',
                              ifelse(mtp_potential$mean < -2, 'weakly metastatic', 
                                     'metastatic')) # 转移潜力

# 食管鳞癌原发灶细胞系(18)
cellline <- filter(mtp_annot,
                   tissue == 'OESOPHAGUS' &
                     cancer_subtype == 'squamous cell carcinoma' &
                     corrected_site_of_origin == 'primary')

exp_ccle_f <- exp_ccle[, c("Symbol", cellline$CCLE_name)]
exp_ccle_f <- exp_ccle_f[complete.cases(exp_ccle_f),]  
rownames(exp_ccle_f) <- exp_ccle_f$Symbol
exp_ccle_f <- exp_ccle_f[,-1] # 表达谱

mtp_potential <- filter(mtp_potential, X1 %in% colnames(exp_ccle))
exp_ccle <- exp_ccle[, c("gene_id", "transcript_ids", "Symbol", mtp_potential$X1)]

mtp_potential_f <- filter(mtp_potential, X1 %in% cellline$CCLE_name) 
exp_ccle_f <- exp_ccle_f[,mtp_potential_f$X1]

# 细胞整合
exp_emt <- as.matrix(GetAssayData(scobj_epi_mes, slot = 'counts'))

gene <- intersect(rownames(exp_emt), rownames(exp_ccle_f)) # 交集基因
exp_merge <- log(cbind(exp_emt[gene,], exp_ccle_f[gene,]) + 1, 2) # 对数化
batch <- c(rep('Patient', ncol(exp_emt)), rep('Cell line', ncol(exp_ccle_f))) # 批次

irlba_res <- prcomp(t(exp_merge)) 
factoextra::fviz_pca_ind(irlba_res,
                         col.ind = batch, 
                         palette = c("#00AFBB",  "#FC4E07"),
                         addEllipses = FALSE, 
                         legend.title = "Groups",
                         repel=FALSE, 
                         geom="point") +
  theme_bw() # 去批次之前的pca结果

zero.rows.lst <- lapply(unique(batch), function(batch_level) {
  if (sum(batch == batch_level) > 1) {
    return(which(apply(exp_merge[, batch == batch_level], 1,
                       function(x) {
                         var(x) == 0
                       })))
  } else {
    return(which(rep(1, 3) == 2))
  }    
})
zero.rows <- Reduce(union, zero.rows.lst)
keep.rows <- setdiff(1:nrow(exp_merge), zero.rows)
exp_merge2 <- exp_merge[keep.rows,] # 去掉方差为的0的基因

exp_merge_rmbatch <- sva::ComBat(dat = exp_merge2, batch = batch) # 去批次

irlba_res_combat <- prcomp(t(exp_merge_rmbatch)) 
factoextra::fviz_pca_ind(irlba_res_combat,
                         col.ind = batch, 
                         palette = c("#00AFBB",  "#FC4E07"),
                         addEllipses = FALSE, 
                         legend.title = "Groups",
                         repel=FALSE, 
                         geom="point")+
  theme_bw() # 去批次之后的pca结果

vardata <- apply(exp_merge_rmbatch, 1, var)
exp_merge_rmbatch2 <- cbind(exp_merge_rmbatch, vardata = vardata)
exp_merge_rmbatch2 <- exp_merge_rmbatch2[order(exp_merge_rmbatch2[,"vardata"], decreasing=T),]
exp_merge_rmbatch_f <- exp_merge_rmbatch2[1:round((nrow(exp_merge_rmbatch2)*30)/100), -ncol(exp_merge_rmbatch2)] # top30%最大变异基因

pca_for_components <- prcomp(t(exp_merge_rmbatch_f))
factoextra::fviz_eig(pca_for_components, ncp = 50, addlabels = T, ggtheme = theme())  # 确定主成分数目

# KNN
ref <- t(pca_for_components$x[1:ncol(exp_emt), 1:19])
que <- t(pca_for_components$x[(ncol(exp_emt)+1):ncol(exp_merge_rmbatch_f), 1:19])
neighbors <- FNN::get.knnx(data = t(ref), query = t(que), k = 50, algorithm = "kd_tree") # 50之后稳定出现三个峰

# 计算拟时序和细胞亚群
id <- apply(neighbors$nn.index, 1, function(x) colnames(ref)[x])
time <- apply(id, 2, function(x) sccds$Pseudotime[colnames(sccds) %in% x])  # 拟时序
cell <- apply(id, 2, function(x) sccds$celltype[colnames(sccds) %in% x]) # 细胞亚群
colnames(cell) <- colnames(time) <- colnames(que)

mtp_potential_f$Pseudotime <- apply(time, 2, function(x) median(x))[mtp_potential_f$X1]
mtp_potential_f$Cell <- apply(cell, 2, function(x) names(sort(table(x), decreasing = T))[1])[mtp_potential_f$X1]

ggplot(mtp_potential_f, aes(x = Pseudotime)) + geom_density()

# 计算EMT分值
anno_ls <- list(anno_hallmark_ls[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
                anno_hallmark_ls[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]],
                anno_hallmark_ls[["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]])
anno_ls[[1]]@setName <- 'EMT'
anno_ls[[2]]@setName <- 'Epithelia'
anno_ls[[3]]@setName <- 'Mesenchymal'
anno_ls[[2]]@geneIds <- epi
anno_ls[[3]]@geneIds <- mes
anno_ls <- GeneSetCollection(anno_ls)
params <- gsvaParam(exprData = as.matrix(exp_ccle_f), 
                    geneSets = anno_ls)
gsva_ccle <- as.data.frame(t(gsva(params)))
gsva_ccle$EMT <- gsva_ccle$Mesenchymal - gsva_ccle$Epithelia

mtp_potential_f$EMT <- gsva_ccle$EMT[match(mtp_potential_f$X1, rownames(gsva_ccle))]
mtp_potential_f$EPI <- gsva_ccle$Epithelia[match(mtp_potential_f$X1, rownames(gsva_ccle))]
mtp_potential_f$MES <- gsva_ccle$Mesenchymal[match(mtp_potential_f$X1, rownames(gsva_ccle))]

cor(mtp_potential_f$EMT, mtp_potential_f$mean) # -0.2228758
cor(mtp_potential_f$EMT, mtp_potential_f$Pseudotime) # 0.624107

# 计算迁移能力
mtp_potential_f$State <- ifelse(mtp_potential_f$Pseudotime <= cutoff_emt1_2, 'EMT_early',
                                ifelse(mtp_potential_f$Pseudotime <= cutoff_emt2_3, 'EMT_stable', 
                                       'EMT_diff'))
table(mtp_potential_f$State)

ggplot(mtp_potential_f, aes(x = Pseudotime, y = mean, colour = State)) +
  geom_point(size = 1) +
  geom_smooth(data = mtp_potential_f, aes(x = Pseudotime, y = mean), method = "lm", linetype = "dashed", color = 'black', se = TRUE) +
  theme_bw()

mtp_potential_f$group <- factor(mtp_potential_f$group, levels = c('non-metastatic', 'weakly metastatic', 'metastatic'))

p7 <- ggplot(mtp_potential_f, aes(x = group, y = EMT, fill = group)) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width = 0.2, linewidth= 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 3) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ggpubr::geom_signif(
    comparisons = list(
      c(1, 2),
      c(2, 3),
      c(1, 3)
    ),
    test = "wilcox.test",
    step_increase = 0.1) +
  labs(x = 'Metastatic group', y = 'EMT score', fill = 'Metastatic group') +
  theme_bw()

mtp_potential_f$State <- factor(mtp_potential_f$State, levels = c('EMT_early', 'EMT_stable', 'EMT_diff'))
p8 <- ggplot(mtp_potential_f, aes(x = State, y = mean, fill = State)) + 
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width = 0.2, linewidth= 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 3) +
  geom_hline(yintercept = -2, linetype = 'dashed') +
  geom_hline(yintercept = -4, linetype = 'dashed') +
  ggpubr::geom_signif(
    comparisons = list(
      c(1, 2),
      c(2, 3),
      c(1, 3)
    ),
    test = "wilcox.test",
    step_increase = 0.1) +
  labs(x = 'EMT state', y = 'Metastatic potential', fill = 'EMT state') +
  scale_fill_manual(values = c('EMT_early' = "#b6dd49", 'EMT_stable' = "#f9e762", 'EMT_diff' = "#74cc64")) + 
  theme_bw()

(p7 + p8)/(p5 + p6)

####
## ESCC的EMT并非与迁移潜力成正比
## EMT-early具有最强的转移能力
## EMT-stable具有最强的分化能力
## EMT-diff具有最强的侵袭能力
####

##### EMT阶段依赖基因 #####
library(ClusterGVis)
library(org.Hs.eg.db)
source('./script/BEAM_crz.R')

BEAM_res <- BEAM_2(sccds, branch_point = 1, cores = 2, progenitor_method = 'duplicate')
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

# genes <- row.names(subset(BEAM_res, qval < 1e-110))
genes <- row.names(subset(BEAM_res, qval < 1e-4))
gene_anno <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = c("GO", "ONTOLOGY"))

gene_tf <- unique(filter(gene_anno, GO == "GO:0003700")[,1]) # 转录因子 
gene_psf <- intersect(genes, c('OVOL1', 'OVOL2', 'OVOL3', 'NUMB', 'NUMBL', 'NFATC1', 'GRHL2', 'TGFBI')) # 表型稳定因子
gene_epi <- intersect(genes, c(epi, genes[grepl('KRT', genes)])) # 角化基因
gene_mes <- intersect(genes, c(mes, genes[grepl('^COL', genes)])) # 间充质基因

plot_genes_branched_heatmap_2(sccds[genes,],
                              branch_point = 1,
                              num_clusters = 3,
                              cores = 2,
                              use_gene_short_name = T,
                              show_rownames = T)



##### Mes受体EMT阶段定位 #####
plots <- lapply(c('IGF1R', 'ITGA3', 'ITGB4', 'PLXNA1', 'CACNA1C'), function(x){
  test <- sccds[,which(GetAssayData(scobj_epi_mes)[x,] != 0)]
  p <- monocle::plot_genes_in_pseudotime(test[x,], color_by = "celltype")+
    scale_color_manual(values = colorRampPalette(brewer.pal(7, "Set2"))(11)) +
    labs(color = 'Mes subclusters')
  return(p)
})

plots[[1]] + plots[[4]] + plots[[2]] + plots[[3]] + plot_layout(ncol = 2, guides = 'collect')
plots[[5]]


##### Mes靶点EMT阶段定位 #####
scobj_epi_mes$EMT_state <- df$EMT_state
scobj_epi_mes_f <- scobj_epi_mes[,which(GetAssayData(scobj_epi_mes)['CACNA1C',] != 0)]

Idents(scobj_epi_mes_f) <- 'EMT_state'
scobj_epi_mes_f$EMT_state <- as.factor(scobj_epi_mes_f$EMT_state)
levels(scobj_epi_mes_f) <- c('EMT_early', 'EMT_stable', 'EMT_diff')
VlnPlot_scCustom(scobj_epi_mes_f, 
                 features = 'CACNA1C', 
                 plot_boxplot = T,
                 colors_use = c('#b6dd49', '#f9e762', '#74cc64')) +
  ggpubr::geom_signif(
    comparisons = list(
      c('Mes_early', 'Mes_stable'),
      c('Mes_stable', 'Mes_diff'),
      c('Mes_early', 'Mes_diff')
    ),
    test = "wilcox.test",
    step_increase = 0.1) +
  labs(x = 'EMT_state')