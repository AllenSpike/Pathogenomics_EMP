## 病理特征无监督聚类表征

library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)
library(moments)
library(clusterProfiler)
library(glmnet)
library(survival)
library(pbapply)
library(ComplexHeatmap)
library(circlize)
library(ggthemes)
library(ggplot2)
library(patchwork)
library(enrichplot)

agg_type = c('mean', 'sd', 'skewness', 'kurtosis') # 聚合类型
path_feat = './data/patch/xj_4096/feature/' # 特征路径
path_output = './data/patch/xj_4096/' # 结果路径


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

#### XJ临床信息 ####
follow_up <- fread('./data/seq/preprocess_result/follow_up.txt') # 随访表
meta <- read.delim('./data/seq/preprocess_result/meta.txt') # 肿瘤元信息表
meta_t <- filter(meta, Group == 'Tumour') # 肿瘤样本
meta_n <- filter(meta, Group == 'Non-tumour') # 癌旁样本
status_xj <- as.data.frame(follow_up[, c(11,9)]) # 随访信息
colnames(status_xj) <- c('time', 'status')
rownames(status_xj) <- follow_up$Pathological_Barcode

feat_xj <- feat_xj[rownames(status_xj),]
feat_xj_norm <- feat_xj_norm[rownames(status_xj),] # 有随访信息的样本

#### 单因素cox特征筛选 ####
hr <- apply(feat_xj_norm, 2, function(x, clin){
  fit <- coxph(Surv(Overall_Survival_Months, Live_Status) ~ x, data = clin)
  fitSummary <- summary(fit)
  res <- c(fitSummary$conf.int[,"exp(coef)"], # 风险比（HR）
           fitSummary$coefficients[,"Pr(>|z|)"], # HR p值
           fitSummary$conf.int[,"lower .95"], # HR 95%置信区间低值
           fitSummary$conf.int[,"upper .95"],  # HR 95%置信区间高值
           fitSummary$concordance[["C"]], # 一致性指数（C-index）
           fitSummary$concordance[["se(C)"]], # C-index 标准差
           fitSummary$concordance[["C"]] - 1.96 * fitSummary$concordance[["se(C)"]], # C-index 95%置信区间低值
           fitSummary$concordance[["C"]] + 1.96 * fitSummary$concordance[["se(C)"]] # C-index 95%置信区间高值
  )
  names(res) <- c('HR', 'HR_pvalue', 'HR_L95CI', 'HR_H95CI',
                  'C-index', 'C-index_SE', 'C-index_L95CI', 'C-index_H95CI')
  return(res)
  return(res)
}, clin = follow_up) 

hr <- as.data.frame(t(hr))
hr$coef <- log(hr$HR)

# 保存hr结果
wb <- createWorkbook() # 创建工作簿
addWorksheet(wb, 'HR_xj') # 加入工作表
writeData(wb, 'HR_xj', hr, rowNames = T, colNames = T) # 写入数据
saveWorkbook(wb, paste0("./2_unicox_result.xlsx")) # 保存工作簿

#### 显著预后特征图 ####
# 筛选显著预后特征
hr_f <- filter(hr, HR_pvalue < 0.05)
hr_f$interaction <- str_split(rownames(hr_f), '_', simplify = T)[,1]
hr_f$color <- ifelse(hr_f$interaction == 'T-T', '#f14166',
                     ifelse(hr_f$interaction == 'I-I', '#6a75c2',
                            ifelse(hr_f$interaction == 'S-S', '#f9a411', 
                                   ifelse(hr_f$interaction == 'T-I', '#ae5b94',
                                          ifelse(hr_f$interaction == 'T-S', '#f5733c',
                                                 ifelse(hr_f$interaction == 'I-S', '#b28d6a', NA))))))
hr_f$type <- unlist(lapply(strsplit(rownames(hr_f), "_"), function(x) tail(x, 1)))
hr_f <- arrange(hr_f, desc(coef))
hr_f$varname <- rownames(hr_f)

# 柱状图
p1 <- ggplot(hr_f, aes(x = reorder(rownames(hr_f), coef), 
                       y = coef, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "right", 
        axis.text.y = element_text(size = 7, 
                                   color = rev(hr_f$color)
        )) +
  labs(title = "", x = "Feature", y = "Univariate cox regression coefficient", fill = 'Aggregation Type') +
  scale_fill_manual(values = c("mean" = "#6e2f9f", "sd" = "#ff00ff", "skewness" = "#0000ff", "kurtosis" = "#008281"))

p2 <- ggplot(hr_f, aes(x = reorder(rownames(hr_f), coef), 
                       y = coef, color = interaction)) + 
  geom_point() + 
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
  labs(title = "", x = "", y = "", color = 'Interaction') + 
  scale_color_manual(values = c("T-T" = "#f14166", "I-I" = "#6a75c2", "S-S" = "#f9a411", "T-I" = "#ae5b94", 'T-S' = '#f5733c','I-S' = '#b28d6a')) + 
  theme(legend.position = c(0.5, 0), legend.direction = 'horizontal')

pdf('./fig1.unicox_significant_barplot.pdf', width = 7, height = 9.5)
p1 + p2 + plot_layout(heights = c(1, 0))
dev.off()

# 森林图
hr_f <- hr_f[order(hr_f$varname),]
p1 <- ggplot(hr_f, aes(y = varname, x = HR)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 7.5, alpha = 0.2, fill = "#6a75c2") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 7.5, ymax = 22.5, alpha = 0.2, fill = '#b28d6a') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 22.5, ymax = 37.5, alpha = 0.2, fill = '#f9a411') +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 37.5, ymax = 49.5, alpha = 0.2, fill = "#ae5b94") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 49.5, ymax = 65.5, alpha = 0.2, fill = '#f5733c') +
  geom_point(size = 1.5, color = "orange") +
  geom_errorbarh(aes(xmax = HR_H95CI, xmin = HR_L95CI), size= 0.7, height = 0.3, colour = "orange") +
  geom_vline(aes(xintercept = 1), color="gray", linetype="dashed", size = 0.7)+
  ylab('Feature')+
  xlab('Hazard ratio (HR)')+
  theme_few()+
  theme(axis.text.y = element_text(size = 9, color = "black"))+
  theme(axis.text.x = element_text(size = 9, color = "black"))+
  theme(title=element_text(size = 12)) +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(l = 30, t = 10, r = 30),
        axis.text.x = element_text(angle = 45, hjust = 1, 
                                   size = 5, 
                                   color = hr_f$color)) 

p2 <- ggplot(hr_f, aes(x = reorder(rownames(hr_f), coef), 
                       y = coef, color = interaction)) + 
  geom_point() + 
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
  labs(title = "", x = "", y = "", color = 'Interaction') + 
  scale_color_manual(values = c("T-T" = "#f14166", "I-I" = "#6a75c2", "S-S" = "#f9a411", "T-I" = "#ae5b94", 'T-S' = '#f5733c','I-S' = '#b28d6a')) + 
  theme(legend.position = c(0.5, 0), legend.direction = 'horizontal')

pdf('./2_unicox_significant_forest.pdf', width = 12, height = 7)
p1 + p2 + plot_layout(heights = c(1, 0))
dev.off()


#### 显著预后特征数值分布图 ####
ggplot(hr_f, aes(x = coef, y = interaction, fill = interaction)) +
  ggridges::geom_density_ridges_gradient() +
  theme_bw() +
  theme(legend.position = "right") +
  labs(title = "", x = "Coefficient", y = "Interaction", fill = 'Interaction') +
  scale_fill_manual(values = c("T-T" = "#f14166", "I-I" = "#6a75c2", "S-S" = "#f9a411", "T-I" = "#ae5b94", 'T-S' = '#f5733c','I-S' = '#b28d6a'))


#### 显著预后特征相关性热图 ####
feat_xj_norm <- feat_xj_norm[,order(colnames(feat_xj_norm))]
feat_cor <- cor(feat_xj_norm, use = "pairwise.complete.obs", method = 'pearson')
is.na(feat_cor) <- 0 # 将NA填充为0
# feat_cor <- abs(feat_cor)
feat_cor <- feat_cor[hr_f$varname, hr_f$varname]

index <- unlist(lapply(unique(str_split(colnames(feat_cor), '_', simplify = T)[,1]), function(type){
  feat_cor_f <- feat_cor[str_detect(colnames(feat_cor), type), str_detect(colnames(feat_cor), type)]
  hclust_result <- hclust(dist(feat_cor_f), method = "complete") # 层次聚类
  index <- colnames(feat_cor_f)[hclust_result$order]
  return(index)
})) # 在每种交互类型内部进行聚类
feat_cor <- feat_cor[index, index] 

row_annotation <- data.frame(
  Interaction = str_split(colnames(feat_cor), '_', simplify = T)[,1],
  check.names = F
)
rownames(row_annotation) <- colnames(feat_cor)
col_annotation <- row_annotation

col_annotation <- HeatmapAnnotation(df = col_annotation,
                                    show_annotation_name = F,
                                    annotation_legend_param = list(Interaction = list(direction = "vertical")),
                                    col = list(Interaction = c("T-T" = "#f14166", "I-I" = "#6a75c2", "S-S" = "#f9a411", "T-I" = "#ae5b94", 'T-S' = '#f5733c','I-S' = '#b28d6a')))
row_annotation <- rowAnnotation(df = row_annotation,
                                show_annotation_name = F,
                                annotation_legend_param = list(Interaction = list(direction = "vertical")),
                                col = list(Interaction = c("T-T" = "#f14166", "I-I" = "#6a75c2", "S-S" = "#f9a411", "T-I" = "#ae5b94", 'T-S' = '#f5733c','I-S' = '#b28d6a')))
ht <- Heatmap(abs(feat_cor),
              name = 'Absolute Pearson R',
              col = colorRamp2(c(0, 0.5, 1), c('white' ,'#ffe888', '#920024')),
              use_raster = T,
              show_row_names = F,
              show_column_names = F,
              cluster_rows = F,
              cluster_columns = F,
              left_annotation = row_annotation,
              top_annotation = col_annotation,
              heatmap_legend_param = list(direction = 'vertical'))

pdf('./fig1.histomic_correlation.pdf', width = 7, height = 5)
draw(ht, merge_legend = T, heatmap_legend_side = 'right', annotation_legend_side = "right")
dev.off()

#### 细胞交互富集分析 ####
## HR富集
val <- hr$coef
names(val) <- rownames(hr)
val <- sort(val, decreasing = T) 

inter_df <- data.frame(
  term = paste(str_split(rownames(hr), '_', simplify = T)[,1],
               unlist(lapply(str_split(rownames(hr), '_'), function(x) tail(x, 1))),
               sep = '_'),
  # term = str_split(rownames(hr), '_', simplify = T)[,1],
  gene = rownames(hr)
)

inter_res <- GSEA(geneList = val,
                  TERM2GENE = inter_df,
                  minGSSize = 1,
                  maxGSSize = Inf,
                  pvalueCutoff = 1,
                  eps = 0,
                  verbose = T) 
View(inter_res@result)

pdf('./supp1_enrichment_significant_HR.pdf', width = 7, height = 7)
gseaplot2(inter_res,
          c('I-S_sd', 'S-S_sd', 'I-I_sd',
            'T-S_kurtosis', 'I-S_kurtosis', 'T-I_kurtosis'),
          color = c('#3f8cc2', '#f0464e', '#2c797d', '#ff9a35', '#1db648','#6a75c2'),
          # inter_res@result$ID,
          pvalue_table = F,
          base_size = 12,
          ES_geom = "line", 
          subplots = 1:3, 
          rel_heights = c(0.6, 0.2, 0.2)
)
dev.off()

plot <- inter_res@result
plot$label <- substr(plot$ID, 1, 3)

p1 <- ggplot(plot, aes(x = NES, y = reorder(ID, NES), fill = label)) + 
  geom_col() +  
  scale_fill_manual(values = c('T-T' = '#f14166', 'I-I' = '#6a75c2', 'S-S' = '#f9a411',
                               'T-I' = '#ae5b94', 'T-S' = '#f5733c', 'I-S' = '#b28d6a')) +
  theme_bw() + 
  labs(y = 'Feature type', fill = 'Interaction type')


## 变异富集
res_diff <- apply(feat_xj, 2, function(x){
  df <- data.frame(
    sd = sd(x, na.rm = T), # 标准差
    vc = abs(sd(x, na.rm = T) / mean(x, na.rm = T)), # 变异系数
    pval_kruskal = kruskal.test(as.list(x))[['p.value']] # kruskal test p值
  )
  return(df)
})
res_diff <- as.data.frame(do.call(rbind, res_diff))

val <- res_diff$vc
names(val) <- rownames(res_diff)
val <- sort(val, decreasing = T) 

inter_df <- data.frame(
  term = paste(str_split(rownames(res_diff), '_', simplify = T)[,1], 
               unlist(lapply(str_split(rownames(res_diff), '_'), function(x) tail(x, 1))), 
               sep = '_'),
  gene = rownames(res_diff)
)

inter_res2 <- GSEA(geneList = val,
                   TERM2GENE = inter_df,
                   minGSSize = 1,
                   maxGSSize = Inf,
                   pvalueCutoff = 1,
                   scoreType = "pos",
                   eps = 0,
                   verbose = T,
                   nPermSimple = 10000) 
View(inter_res2@result)

pdf('./supp1_feature_enrichment_VC.pdf', width = 7, height = 7)
gseaplot2(inter_res2, 
          c('S-S_skewness', 'T-S_skewness', 'T-T_skewness'), 
          color = c('#ff9a35', '#3f8cc2', '#f0464e'),
          pvalue_table = F,
          base_size = 12,
          ES_geom = "line", 
          subplots = 1:3, 
          rel_heights = c(0.6, 0.2, 0.2)
)
dev.off()

plot2 <- inter_res2@result
plot2$label <- substr(plot2$ID, 1, 3)

p2 <- ggplot(plot2, aes(x = NES, y = reorder(ID, NES), fill = label)) + 
  geom_col() +  
  scale_fill_manual(values = c('T-T' = '#f14166', 'I-I' = '#6a75c2', 'S-S' = '#f9a411',
                               'T-I' = '#ae5b94', 'T-S' = '#f5733c', 'I-S' = '#b28d6a')) +
  theme_bw() + 
  labs(y = 'Feature type', fill = 'Interaction type')


p1 + p2 + plot_layout(guides = 'collect')

# 保存细胞交互富集结果
wb <- createWorkbook() # 创建工作簿
addWorksheet(wb, 'HR_xj')
writeData(wb, 'HR_xj', inter_res, rowNames = F, colNames = T)
addWorksheet(wb, 'VC_xj')
writeData(wb, 'VC_xj', inter_res2, rowNames = F, colNames = T)
saveWorkbook(wb, paste0("./2_enrichment_result.xlsx")) # 保存工作簿

####
rm(list = ls())
gc()
