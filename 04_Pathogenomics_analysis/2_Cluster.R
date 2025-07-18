## 病理特征无监督聚类

dyn.load('/home/ljc/anaconda3/lib/libhdf5_hl.so.200') # 设置hdf5r的配置文件

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
library(ggrepel)
library(easier)
library(Seurat)
library(SeuratDisk)
library(ggExtra)
library(pbapply)
library(factoextra)
library(ggVolcano)

path_output = './data/patch/xj_4096/'
agg_type = c('mean', 'sd', 'skewness', 'kurtosis')


# 数据读入、分值计算、数据汇总 -------------------------------------------------
#### 1、读入数据 ####
# 背景文件
ncbi <- read.delim('./data/seq/Homo_sapiens.gene_info', check.names = F) # 基因文件
genelocate <- read.table(system.file("extdata","genelocate.txt",package="multiOmicsViz"), 
                         header = T, sep = '\t',stringsAsFactors=FALSE, check.names=FALSE) # 基因位置文件
gene_esca <- read.xlsx('./data/seq/resource/esca_related_gene_summary.xlsx') # esca相关基因: 来自DisGeNet(https://www.disgenet.org/)
gene_driver <- c('KRT5', 'CDH10', 'LILRB3', 'YEATS2', 'CASP8',
                 'TP53', 'ZNF750', 'NOTCH1', 'FAT1', 'NFE2L2') # ESCC驱动基因（来自2020-CellRes-Whole-genome sequencing of 508 patients identifies key molecular features associated with poor prognosis in esophageal squamous cell carcinoma）
# load('./data/seq/preprocess_result/lincs_knockdown_escc.RData') # escc基因敲除谱

anno_go <- read.gmt('./data/seq/resource/go.v2023.2.Hs.entrez.gmt')
anno_hallmark <- read.gmt('./data/seq/resource/hallmark.v2023.2.Hs.entrez.gmt')
anno_kegg <- read.gmt('./data/seq/resource/kegg_legacy.v2023.2.Hs.entrez.gmt')
anno_wiki <- read.gmt('./data/seq/resource/wikipathways.v2023.2.Hs.entrez.gmt')
anno_reactome <- read.gmt('./data/seq/resource/reactome.v2023.2.Hs.entrez.gmt')

# 元信息表
meta <- read.delim('./data/seq/preprocess_result/meta.txt')
meta_t <- filter(meta, Group == 'Tumour') # 肿瘤样本
meta_n <- filter(meta, Group == 'Non-tumour') # 癌旁样本

# 随访表
follow_up <- fread('./data/seq/preprocess_result/follow_up.txt')

# snv
load('./data/seq/preprocess_result/output_maf.RData')

# cnv
cnv <- read.delim('./data/seq/preprocess_result/all_data_by_genes.txt')
cnv$Gene.Symbol <- gsub('\\|chr\\d', '', cnv$Gene.Symbol)
laml.gistic <- readGistic(gisticAllLesionsFile = './data/seq/preprocess_result/all_lesions.conf_90.txt', 
                          gisticAmpGenesFile = './data/seq/preprocess_result/amp_genes.conf_90.txt', 
                          gisticDelGenesFile = './data/seq/preprocess_result/del_genes.conf_90.txt', 
                          gisticScoresFile = './data/seq/preprocess_result/scores.gistic')

# rna
exp_count <- read.csv('./data/seq/preprocess_result/exp_count.csv') 
exp_tpm <- read.csv('./data/seq/preprocess_result/exp_tpm.csv') 
exp_vsd <- read.csv('./data/seq/preprocess_result/exp_vsd.csv')
rownames(exp_count) <- rownames(exp_tpm) <- rownames(exp_vsd) <- exp_count$X
exp_count <- exp_count[,-1]
exp_tpm <- exp_tpm[,-1]
exp_vsd <- exp_vsd[,-1]

exp_count_t <- exp_count[,meta_t$Tumor_Sample_Barcode]
exp_tpm_t <- exp_tpm[,meta_t$Tumor_Sample_Barcode]
exp_vsd_t <- exp_vsd[,meta_t$Tumor_Sample_Barcode] # 肿瘤样本

exp_count_n <- exp_count[,meta_n$Tumor_Sample_Barcode]
exp_tpm_n <- exp_tpm[,meta_n$Tumor_Sample_Barcode]
exp_vsd_n <- exp_vsd[,meta_n$Tumor_Sample_Barcode] # 癌旁样本


#### 2、计算分值 ####
# 计算肿瘤突变负荷
tmbScore <- tmb(maf = laml)

# 计算T细胞免疫功能障碍评分
# exp_tpm_t_tide <- t(apply(exp_tpm_t, 1, function(x){x - mean(x, na.rm = T)}))
# write.table(exp_tpm_t_tide, './data/seq/preprocess_result/exp_tpm_t_tide.txt', sep = '\t', quote = F, row.names = T, col.names = T)
tideScore <- read.csv('./data/seq/preprocess_result/tide.csv')[,c(1,4,11,12)] # 来自http://tide.dfci.harvard.edu/
colnames(tideScore)[1] <- 'Tumor_Sample_Barcode'

# 计算溶细胞活性
caScore <- as.data.frame(apply(exp_tpm_t[c('GZMA', 'PRF1'),], 2, function(x){mean(x, na.rm = T)}))
caScore$Tumor_Sample_Barcode <- rownames(caScore)
colnames(caScore) <- c('CytolyticActivity', 'Tumor_Sample_Barcode')

# 计算肿瘤纯度/免疫评分/基质评分
estimate = function(dat,pro) {
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")   ## platform
  scores=read.table(output.ds,skip = 2,header = T,check.names = F)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  library(stringr)
  rownames(scores)=str_replace_all(rownames(scores),'[.]','-') # 这里TCGA样本名里面的-变成.了，进行恢复
  write.csv(scores,file="Stromal_Immune_ESTIMATE.Score.csv")
  return(scores)
}

esScore <- estimate(exp_tpm_t, 'escc')
esScore <- read.csv('./Stromal_Immune_ESTIMATE.Score.csv', row.names = 1)
esScore$Tumor_Sample_Barcode <- rownames(esScore)


### 3、数据汇总 ####
meta <- read.delim('./data/seq/preprocess_result/meta.txt')

# 将上述信息汇总到元信息表
meta <- meta %>% 
  left_join(follow_up[,-1], riskScore, by = 'Pathological_Barcode') %>%
  left_join(esScore, by = 'Tumor_Sample_Barcode') %>%
  left_join(tideScore, by = 'Tumor_Sample_Barcode') %>%
  left_join(tmbScore, by = 'Tumor_Sample_Barcode') %>%
  left_join(caScore, by = 'Tumor_Sample_Barcode')

meta$Esophgus_Location[meta$Esophgus_Location == '0'] <- 'middle'  # 有一个记录错误，先将其修正为middle
meta$Pathological_Barcode <- as.integer(meta$Pathological_Barcode)

meta_t <- filter(meta, Group == 'Tumour') # 肿瘤样本
meta_n <- filter(meta, Group == 'Non-tumour') # 癌旁样本

status_xj <- as.data.frame(follow_up[, c(11,9)]) # 随访信息
rownames(status_xj) <- follow_up$Pathological_Barcode


# analysis ----------------------------------------------------------------
#### 1、病理特征无监督聚类 ####
# 特征矩阵
feat_xj <- t(do.call(rbind, lapply(agg_type, function(x){
  s <- read.xlsx(paste0(path_output,"/feat_agg.xlsx"), rowNames = T, sheet = x)[1:252,] # 选择前252个特征
  rownames(s) <- paste0(rownames(s), "_", x)
  return(s)
}
)))
colnames(feat_xj) <- gsub('GraphAvg_', 'GraphAvg\\.', colnames(feat_xj)) 
feat_xj <- feat_xj[rownames(status_xj),] # 有随访信息的样本
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
# feat_xj_norm <- feat_xj[,-which(colSums(feat_xj_norm, na.rm = T) == 0)] # 排除全为NA的特征

# pca
pca <- prcomp(feat_xj_norm, scale = F, center = F)
pca_data <- pca$x[, 1:2] # 提取前2个主成分

# hclust
sorted_contributions <- sort(abs(pca$rotation[,1]), decreasing = TRUE)
top_features <- names(sorted_contributions[1:32]) # 提取对PC1贡献最大的前50个特征

heat <- feat_xj_norm[,top_features]
heat[heat > 2] <- 2
heat[heat < -2] <- -2
distance_matrix <- dist(heat)
hclust_result <- hclust(distance_matrix, method = "complete")
plot(hclust_result)

group <- cutree(hclust_result, k = 2)
table(group)
group <- ifelse(group == 1, 2, 1) # 修正聚类名称

# pca
pca_df <- as.data.frame(pca_data)
pca_df$Cluster <- as.factor(group)
summ <- summary(pca)
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

pdf('./fig1.cluster_pca_xj.pdf', width = 7, height = 6)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 1.5) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c('1' = '#f14166', '2' = '#6a75c2')) +
  labs(x = xlab, y = ylab , title = "") +
  theme_bw()
dev.off()

# 绘制KM曲线
sf <- as.data.frame(group)
sf$case_submitter_id <- rownames(sf)
colnames(sf)[1] <- "Cluster"

status <- status_xj
status$case_submitter_id <- rownames(status)
status <- left_join(status, sf, by = "case_submitter_id")
colnames(status)[1:2] <- c('time', 'status')
fit <- survfit(Surv(time, status) ~ Cluster, data = status)

pdf('./fig1.cluster_km_xj.pdf', width = 7, height = 7)
ggsurvplot(fit = fit, 
           data = status,
           palette = c('#f14166', '#6a75c2'),
           pval = T,
           risk.table = T,
           surv.median.line = "hv", #添加中位生存曲线
           title = "Overall survival", #标题
           ylab = "Cumulative survival (percentage)", xlab = " Time (Months)", #更改横纵坐标
           censor.shape = 124,
           censor.size = 4,
           conf.int = FALSE
)
dev.off()


#### 2、无监督组之间的差异病理特征 ####
df <- as.data.frame(group)
colnames(df) <- 'Cluster'
dm <- dist(heat[rownames(df)[df$Cluster == '1'],], method = "euclidean")
hc_1 <- hclust(dm)
dm <- dist(heat[rownames(df)[df$Cluster == '2'],], method = "euclidean")
hc_2 <- hclust(dm)

ord <- c(rownames(df)[df$Cluster == '1'][hc_1$order], 
         rownames(df)[df$Cluster == '2'][hc_2$order]
) 
meta_t <- meta_t[match(ord, meta_t$Pathological_Barcode),] 
meta_t$Cluster <- df$Cluster[match(ord, rownames(df))] # 样本排序

col_annotation <- HeatmapAnnotation(df = data.frame(Cluster = meta_t$Cluster, row.names = meta_t$Pathological_Barcode),
                                    show_annotation_name = F,
                                    col = list(Cluster = c('1' = '#f14166', '2' = '#6a75c2')))
heat <- heat[ord,]
row_colors <- str_split(colnames(heat), '_', simplify = T)[,1]
row_colors <- ifelse(row_colors == 'T-T', '#f14166',
                     ifelse(row_colors == 'I-I', '#6a75c2',
                            ifelse(row_colors == 'S-S', '#f9a411', 
                                   ifelse(row_colors == 'T-I', '#ae5b94',
                                          ifelse(row_colors == 'T-S', '#f5733c',
                                                 ifelse(row_colors == 'I-S', '#b28d6a', NA))))))
ht <- Heatmap(t(heat),
              name = 'Feature score',
              col = colorRampPalette(c("#5f176c", "#009c9a","#f3e700"))(50),
              top_annotation = col_annotation,
              heatmap_legend_param = list(
                direction = "horizontal",   # 图例横向排列
                legend_width = unit(4, "cm"), # 图例宽度
                title = "Feature score" # 图例标题
              ),
              use_raster = T,
              show_row_names = T,
              show_column_names = F, 
              show_row_dend = T, 
              cluster_rows = T,
              cluster_columns = T,  
              show_column_dend = T,
              row_names_side = 'left', 
              row_dend_side = 'right', 
              row_names_gp = gpar(fontsize = 6, col = row_colors))

lgd <- Legend(labels = c('I-I', 'S-S', 'T-I', 'T-S', 'I-S'),
              title = 'Interaction',
              title_position = 'topleft',
              legend_gp = gpar(fill = c('#6a75c2', '#f9a411', '#ae5b94', '#f5733c', '#b28d6a')),
              gap = unit(2, "cm"),
              ncol = 1)

pdf('./fig1.cluster_heatmap_xj.pdf', width = 7, height = 6)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_list = list(lgd)) 
dev.off()

##### top feature intensity #####
top_feature_intensity <- stack(apply(heat, 2, function(x){
  cluster1 <- sum(x[rownames(df)[df$Cluster == '1']])
  cluster2 <- sum(x[rownames(df)[df$Cluster == '2']])
  return(cluster1 - cluster2)
}))

# 保存intensity结果
wb <- createWorkbook() # 创建工作簿
addWorksheet(wb, 'Intensity_xj') # 加入工作表
writeData(wb, 'Intensity_xj', top_feature_intensity, rowNames = F, colNames = T) # 写入数据
saveWorkbook(wb, paste0("./fig1.feature_intensity_result.xlsx")) # 保存工作簿


#### 3、无监督组之间的差异表达基因 ####
meta_t$group <- group
grouplist <- factor(ifelse(meta_t$group == 1, 'case', 'control'), levels = c('case', 'control')) 
grouplist <- relevel(grouplist, ref = 'control') # 选择group2为参考组(预后较好)
colData <- data.frame(row.names = meta_t$Tumor_Sample_Barcode, grouplist = grouplist)
dds <- DESeqDataSetFromMatrix(countData = round(exp_count_t)[,meta_t$Tumor_Sample_Barcode],
                              colData = colData,
                              design = ~ grouplist)
dds <- dds[rowSums(counts(dds)) > 1 ,]
dds <- DESeq(dds)
res <- data.frame(results(dds))
res$gene <- rownames(res)

res <- add_regulate(res, log2FC_name = "log2FoldChange", fdr_name = "padj",log2FC = 1, fdr = 0.05)

#### Vocalno ####
ggvolcano(res,
          x = "log2FoldChange",
          y = "padj",
          label = "gene",
          label_number = 20,
          output = FALSE)+
  theme(
    plot.title = element_blank(),   
    plot.caption = element_blank(),  
    legend.position = "none"
  )


#### ORA ####
res_f <- filter(res, pvalue < 0.05 & log2FoldChange > 0.6) # 上调基因
gene <- filter(ncbi, Symbol %in% rownames(res_f))$GeneID
res_ora <- list(
  go = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_go),
  hallmark = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_hallmark),
  kegg = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_kegg),
  wiki = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_wiki),
  reactome = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_reactome)
)

# 可视化
df <- res_ora$go@result
df <- filter(df, ID %in% c('GOBP_MESENCHYME_DEVELOPMENT',
                           'GOCC_COLLAGEN_CONTAINING_EXTRACELLULAR_MATRIX',
                           'GOCC_CONTRACTILE_FIBER',
                           'GOBP_MUSCLE_SYSTEM_PROCESS',
                           'GOBP_MUSCLE_CELL_DIFFERENTIATION',
                           'GOBP_EPITHELIAL_CELL_PROLIFERATION')) # EMT相关GO term

df$v <- -log10(df$pvalue)
df$Description <- str_to_title(gsub('_', ' ', gsub('^GO\\w\\w_', '', df$Description)))
df$Description <- factor(factor(df$Description), levels = rev(df$Description))
ggplot(data = df,
       aes(x = Description, y = v)) + 
  geom_bar(stat = "identity", position = "dodge", fill = '#f68a85', width = 0.9)+ 
  coord_flip()+
  theme_minimal()+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) +
  labs(x = "GO term", y = "-log10(p-value)")

# 保存ORA富集结果
wb <- loadWorkbook('./fig1.enrichment_result.xlsx')
addWorksheet(wb, 'ORA_xj')
writeData(wb, 'ORA_xj', res_ora$go@result, rowNames = F, colNames = T)
saveWorkbook(wb, paste0("./fig1.enrichment_result.xlsx"), overwrite = T) # 保存工作簿

#### GSEA ####
res_s <- left_join(res, ncbi[,2:3], by = c('gene' = 'Symbol')) %>% arrange(desc(log2FoldChange))
val <- res_s$log2FoldChange
names(val) <- res_s$GeneID
val <- val[!is.na(names(val))]

res_GSEA <- list(
  hallmark = GSEA(geneList = val, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_hallmark),
  go = GSEA(geneList = val, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_go)
)

kf <- res_GSEA$hallmark
kf@result$Description <- str_to_title(gsub('_', ' ', gsub('HALLMARK_', '', kf@result$Description)))

pdf('./fig1.cluster_gsea_xj.pdf', width = 7, height = 5)
gseaplot2(kf,
          c('HALLMARK_MYOGENESIS', 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
            'HALLMARK_ALLOGRAFT_REJECTION', 'HALLMARK_INTERFERON_GAMMA_RESPONSE', 'HALLMARK_INTERFERON_ALPHA_RESPONSE'),
          color = c('#3f8cc2', '#f0464e', '#2c797d', '#ff9a35', '#1db648','#6a75c2'),
          pvalue_table = F,
          base_size = 12,
          ES_geom = "line", 
          subplots = 1:2, 
          rel_heights = c(0.7, 0.3)
)
dev.off()

# 保存GSEA富集结果
wb <- loadWorkbook('./fig1.enrichment_result.xlsx')
addWorksheet(wb, 'GSEA_xj')
writeData(wb, 'GSEA_xj', res_GSEA$hallmark@result, rowNames = F, colNames = T)
saveWorkbook(wb, paste0("./fig1.enrichment_result.xlsx"), overwrite = T) # 保存工作簿

rm(list = ls())
gc()
