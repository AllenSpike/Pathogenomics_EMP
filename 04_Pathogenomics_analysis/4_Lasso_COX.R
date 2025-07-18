## 病理特征预后

library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)
library(survival)
library(survminer)
library(pbapply)
library(ggplot2)
library(SummarizedExperiment)
library(glmnet)
library(caret)
library(clusterProfiler)
library(patchwork)
library(circlize)
library(ComplexHeatmap)

setwd('/home/ljc/0_project/0_ESCC/')

path_output = './data/patch/xj_4096/'
path_output2 = './data/patch/tcga_4096/'
path_output3 = './data/patch/gdph_supp_4096/'

agg_type = c('mean', 'sd', 'skewness', 'kurtosis')


# 数据读入 --------------------------------------------------------------------
#### XJ临床信息 ####
follow_up <- fread('./data/seq/preprocess_result/follow_up.txt') # 随访表
meta <- read.delim('./data/seq/preprocess_result/meta.txt') # 肿瘤元信息表
meta_t <- filter(meta, Group == 'Tumour') # 肿瘤样本
meta_n <- filter(meta, Group == 'Non-tumour') # 癌旁样本

status_xj <- as.data.frame(follow_up[, c(11,9)]) # 随访信息
colnames(status_xj) <- c('time', 'status')
rownames(status_xj) <- follow_up$Pathological_Barcode

#### XJ特征矩阵 ####
feat_xj <- t(do.call(rbind, lapply(agg_type, function(x){
  s <- read.xlsx(paste0(path_output,"/feat_agg.xlsx"), rowNames = T, sheet = x)[1:252,] # 选择前252个特征
  rownames(s) <- paste0(rownames(s), "_", x)
  return(s)
}
)))
colnames(feat_xj) <- gsub('GraphAvg_', 'GraphAvg\\.', colnames(feat_xj)) 
feat_xj <- feat_xj[rownames(status_xj),] # 有随访信息的样本

# 标准化
feat_xj_norm <- scale(feat_xj)
feat_xj_norm <- apply(feat_xj_norm, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

#### XJ分子信息 ####
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

## gsva
ssgsea_xj <- gsva(
  expr = log2(as.matrix(exp_tpm_t) + 1),
  gset.idx = anno_hallmark,
  method = "gsva",
  kcdf = "Gaussian",
  verbose = T
)

clin_xj <- meta_t
clin_xj <- cbind(clin_xj, t(ssgsea_xj[, clin_xj$Tumor_Sample_Barcode]))
rownames(clin_xj) <- clin_xj$Pathological_Barcode
clin_xj <- clin_xj[rownames(feat_xj),]
clin_xj <- clin_xj[,5:ncol(clin_xj)]
colnames(clin_xj)[1] <- 'samp'

#### TCGA临床信息 ####
clinical <- fread('./data/seq/preprocess_result/TCGA-ESCA/clinical.tsv')
clinical$Overall_Survival_Months <- ifelse(clinical$vital_status == 'Dead',
                                           as.numeric(clinical$days_to_death)/30,
                                           as.numeric(clinical$days_to_last_follow_up)/30)
clinical$Live_Status <- ifelse(clinical$vital_status == 'Alive', 0, 1)
clinical <- clinical[!duplicated(clinical$case_submitter_id), ]


#### TCGA特征矩阵 ####
feat_tcga <- t(do.call(rbind, lapply(agg_type, function(x){
  s <- read.xlsx(paste0(path_output2,"/feat_agg.xlsx"), rowNames = T, sheet = x)[1:252,] # 选择前252个特征
  rownames(s) <- paste0(rownames(s), "_", x)
  return(s)
}
)))
rownames(feat_tcga) <- substr(rownames(feat_tcga), 1, 12)
colnames(feat_tcga) <- gsub('GraphAvg_', 'GraphAvg\\.', colnames(feat_tcga))
clinical_f <- clinical[clinical$case_submitter_id %in% rownames(feat_tcga), ] # 有随访信息的样本

feat_tcga <- apply(feat_tcga, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

status_tcga <- as.data.frame(clinical_f[, c(159, 160)])
colnames(status_tcga) <- c('time', 'status')
rownames(status_tcga) <- clinical_f$case_submitter_id

# 标准化
feat_tcga_norm <- scale(feat_tcga)
feat_tcga_norm <- apply(feat_tcga_norm, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na


#### Gdph临床信息 ####
gdph <- read.xlsx('./data/wsi/gdph/GDPH_临床信息_单纯手术_20241106.xlsx')
gdph <- gdph[!duplicated(gdph$病理号),] # 去重
gdph$`生存状态（修正版）` <- ifelse(gdph$`生存状态（修正版）` == 'alive', 0, 1)

status_gdph <- as.data.frame(gdph[, c(114, 106)])
colnames(status_gdph) <- c('time', 'status')
rownames(status_gdph) <- gdph$病理号


#### Gdph特征矩阵 ####
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


# lassocox --------------------------------------------------------------
cor_matrix <- cor(feat_xj_norm, use = "pairwise.complete.obs", method = 'pearson')
high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9)
feature <- setdiff(colnames(feat_xj), colnames(feat_xj)[high_cor_vars]) # 低相关特征

df <- data.frame()
cvs <- list()

for (i in seq(1:25)) {
  
  set.seed(i) # 遍历随机种子, 调整交叉验证的训练集和验证集划分
  cvs[[as.character(i)]] <- list()
  
  for (a in seq(0.2, 1, length=9)) { # 遍历alpha值
    
    # lasso-cox
    cv <- cv.glmnet(feat_xj[,feature],
                    Surv(status_xj$time, status_xj$status), 
                    family = "cox", 
                    type.measure = "C", 
                    nfolds = 5, 
                    alpha = a,
                    standardize = TRUE)
    
    coefficient <- coef(cv, s = "lambda.1se")
    Active.index <- which(as.numeric(coefficient) != 0)
    Active.coefficient <- as.numeric(coefficient)[Active.index]
    sig_mult_cox <- rownames(coefficient)[Active.index]
    
    # XJ预后测试
    riskScore <- predict(cv, newx = feat_xj[,feature], type = "link", s = "lambda.1se")
    riskScore <- as.data.frame(riskScore)
    colnames(riskScore)[1] <- "riskScore"
    riskScore$sample <- as.integer(rownames(riskScore))
    riskScore_cli <- inner_join(riskScore, follow_up, by = c('sample' = "Pathological_Barcode")) # 风险值
    cutoff <- median(riskScore_cli$riskScore, na.rm = T) # 风险组卡值
    riskScore_cli$riskGroup <- ifelse(riskScore_cli$riskScore > cutoff, "High","Low") # 风险组
    if (length(unique(riskScore_cli$riskGroup)) == 1) {
      pvalue <- 1
    } else {
      pvalue <- survdiff(Surv(Overall_Survival_Months, Live_Status) ~ riskGroup, data = riskScore_cli)[["pvalue"]] 
    } # 预后p值
    
    # TCGA队列预后测试
    riskScore_tcga <- predict(cv, newx = feat_tcga[,feature], type = "link", s = "lambda.1se")
    riskScore_tcga <- as.data.frame(riskScore_tcga)
    colnames(riskScore_tcga)[1] <- "riskScore"
    riskScore_tcga$case_submitter_id <- rownames(riskScore_tcga)
    riskScore_cli_tcga <- inner_join(riskScore_tcga, clinical_f, by = "case_submitter_id") # 风险值
    riskScore_cli_tcga$riskGroup <- ifelse(riskScore_cli_tcga$riskScore > cutoff, "High","Low") # 风险组
    if (length(unique(riskScore_cli_tcga$riskGroup)) == 1) {
      pvalue_tcga <- 1
    } else {
      pvalue_tcga <- survdiff(Surv(Overall_Survival_Months, Live_Status) ~ riskGroup, data = riskScore_cli_tcga)[["pvalue"]]
    } # 预后p值
    
    # GDPH队列预后测试
    riskScore_gdph <- predict(cv, newx = feat_gdph[,feature], type = "link", s = "lambda.1se") # 修订
    riskScore_gdph <- as.data.frame(riskScore_gdph)
    colnames(riskScore_gdph)[1] <- "riskScore"
    riskScore_gdph$病理号 <- as.numeric(rownames(riskScore_gdph))
    riskScore_cli_gdph <- inner_join(riskScore_gdph, gdph[, c(5, 114, 106)], by = "病理号") # 风险值
    riskScore_cli_gdph$riskGroup <- ifelse(riskScore_cli_gdph$riskScore > cutoff, "High","Low") # 风险组
    if (length(unique(riskScore_cli_gdph$riskGroup)) == 1) {
      pvalue_gdph <- 1
    } else {
      pvalue_gdph <- survdiff(Surv(OS_month, `生存状态（修正版）`) ~ riskGroup, data = riskScore_cli_gdph)[["pvalue"]]
    } # 预后p值
    
    df <- rbind(df, data.frame(i, a, length(sig_mult_cox), pvalue, pvalue_tcga, pvalue_gdph))
    cvs[[as.character(i)]][[as.character(a)]] <- cv # 保存模型
    print(paste(i, a, length(sig_mult_cox), pvalue, pvalue_tcga, pvalue_gdph, collapse = ' \\t '))
    
  }
}

## 保存结果
save(cvs, df, file = './fig3.lassocox_result_lambda1se.RData')


# Analysis --------------------------------------------------------------------
load('./fig3.lassocox_result_lambda1se.RData')

df_f <- filter(df, pvalue < 0.05 & pvalue_tcga < 0.05 & pvalue_gdph < 0.05)
cv_final <- cvs[['2']][['0.4']] # 选择模型

#### 多因素COX ####
coefficient <- coef(cv_final, s = "lambda.1se")
Active.index <- which(as.numeric(coefficient) != 0)
Active.coefficient <- as.numeric(coefficient)[Active.index]
sig_mult_cox <- rownames(coefficient)[Active.index]

multicox <- coxph(Surv(status_xj$time, status_xj$status) ~ ., data = as.data.frame(feat_xj_norm[,sig_mult_cox]))
summ <- summary(multicox)
summ

# 柱状图
ds <- stack(fit$coefficients)
ds$ind <- gsub('`', '', ds$ind)
colnames(ds) <- c('coef', 'feature')
ds$type <- unlist(lapply(strsplit(ds$feature, "_"), function(x) tail(x, 1)))
ds$interaction <- substr(ds$feature, 1, 3)
ds$color <- ifelse(ds$interaction == 'T-T', '#f14166',
                   ifelse(ds$interaction == 'I-I', '#6a75c2',
                          ifelse(ds$interaction == 'S-S', '#f9a411', 
                                 ifelse(ds$interaction == 'T-I', '#ae5b94',
                                        ifelse(ds$interaction == 'T-S', '#f5733c',
                                               ifelse(ds$interaction == 'I-S', '#b28d6a', NA))))))
ds <- ds[order(ds$coef),]

p1 <- ggplot(ds, aes(x = reorder(feature, coef),
                     y = coef, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 7,
                                   color = ds$color
        )) +
  labs(title = "", x = "Feature", y = "Multiple cox regression coefficient", fill = 'Aggregation Type') +
  scale_fill_manual(values = c("mean" = "#6e2f9f", "sd" = "#ff00ff", "skewness" = "#0000ff", "kurtosis" = "#008281"))

p2 <- ggplot(ds, aes(x = reorder(feature, coef),
                     y = coef, color = interaction)) +
  geom_point() +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  labs(title = "", x = "", y = "", color = 'Interaction') +
  scale_color_manual(values = c("T-T" = "#f14166", "I-I" = "#6a75c2", "S-S" = "#f9a411", "T-I" = "#ae5b94", 'T-S' = '#f5733c','I-S' = '#b28d6a')) +
  theme(legend.position = c(0.5, 0), legend.direction = 'horizontal')

pdf('./fig3.multiple_cox_coefficient.pdf', width = 6, height = 5)
p1 + p2 + plot_layout(heights = c(1, 0))
dev.off()


#### KM ####
# XJ
riskScore <- predict(cv_final, newx = feat_xj[,feature], type = "link", s = "lambda.1se")
riskScore <- as.data.frame(riskScore)
colnames(riskScore)[1] <- "riskScore"
riskScore$sample <- as.integer(rownames(riskScore))
riskScore_cli <- inner_join(riskScore, follow_up, by = c('sample' = "Pathological_Barcode")) # 风险值
cutoff <- median(riskScore_cli$riskScore, na.rm = T) # 风险组卡值
riskScore_cli$RiskGroup <- ifelse(riskScore_cli$riskScore > cutoff, "High","Low")

fit <- survfit(Surv(Overall_Survival_Months, Live_Status) ~ RiskGroup, data = riskScore_cli)

pdf('./fig3.lassocox_km_xj', width = 7, height = 7)
ggsurvplot(fit = fit, 
           data = riskScore_cli,
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

# 保存样本分类，使用clam训练模型
dh <- data.frame(
  case_id = riskScore_cli$sample,
  slide_id = riskScore_cli$sample,
  label = riskScore_cli$RiskGroup
)

write.csv(dh, file = '/home/ljc/0_project/0_ESCC/0_CLAM/dataset_csv/xj_tumor_vs_normal.csv', row.names = F)

# TCGA
riskScore_tcga <- predict(cv_final, newx = feat_tcga[,feature], type = "link", s = "lambda.1se")
riskScore_tcga <- as.data.frame(riskScore_tcga)
colnames(riskScore_tcga)[1] <- "riskScore"
riskScore_tcga$case_submitter_id <- rownames(riskScore_tcga)
riskScore_cli_tcga <- inner_join(riskScore_tcga, clinical_f, by = "case_submitter_id") # 风险值
riskScore_cli_tcga$RiskGroup <- ifelse(riskScore_cli_tcga$riskScore > cutoff, "High","Low") # 风险组

fit <- survfit(Surv(Overall_Survival_Months, Live_Status) ~ RiskGroup, data = riskScore_cli_tcga)

pdf('./fig3.lassocox_km_tcga', width = 7, height = 7)
ggsurvplot(fit = fit, 
           data = riskScore_cli_tcga,
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

# EMT分值和风险值的相关性(加入TCGA?)
dg <- data.frame(
  EMT = clin_xj$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  RiskScore = riskScore_cli$riskScore,
  RiskGroup = riskScore_cli$RiskGroup
)

pdf('./fig3.corr_EMT_RiskScore.pdf', width = 7, height = 5)
ggplot(aes(RiskScore, EMT), data = dg) + 
  geom_point(aes(color = RiskGroup)) + 
  scale_color_manual(values = c('#f14166', '#6a75c2')) +
  geom_smooth(method = 'lm', level = 0.95) + 
  stat_cor(method = 'spearman', color = 'black') +
  theme_minimal()
dev.off()

# GDPH
riskScore_gdph <- predict(cv_final, newx = feat_gdph[,feature], type = "link", s = "lambda.1se")
riskScore_gdph <- as.data.frame(riskScore_gdph)
colnames(riskScore_gdph)[1] <- "riskScore"
riskScore_gdph$病理号 <- as.numeric(rownames(riskScore_gdph))
riskScore_cli_gdph <- inner_join(riskScore_gdph, gdph[, c(5, 114, 106)], by = "病理号") # 风险值
riskScore_cli_gdph$RiskGroup <- ifelse(riskScore_cli_gdph$riskScore > cutoff, "High","Low") # 风险组

fit <- survfit(Surv(OS_month, `生存状态（修正版）`) ~ RiskGroup, data = riskScore_cli_gdph)

pdf('./fig3.lassocox_km_gdph', width = 7, height = 7)
ggsurvplot(fit = fit, 
           data = riskScore_cli_gdph,
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


#### 注释热图 ####
tn <- read.xlsx('./data/seq/TN分期.xlsx', sheet = 1) # TNM分期表
feat_f <- as.data.frame(feat_xj_norm[,sig_mult_cox])
feat_f$Pathological_Barcode <- as.integer(rownames(feat_f))
feat_f2 <- data.frame(
  `I-I` = rowMeans(feat_f[,grep('I-I', colnames(feat_f))]),
  `I-S` = rowMeans(feat_f[,grep('I-S', colnames(feat_f))]),
  `S-S` = rowMeans(feat_f[,grep('S-S', colnames(feat_f))]),
  `T-I` = rowMeans(feat_f[,grep('T-I', colnames(feat_f))]),
  `T-S` = rowMeans(feat_f[,grep('T-S', colnames(feat_f))]),
  Pathological_Barcode = feat_f$Pathological_Barcode,
  check.names = F
)

meta <- read.delim('./data/seq/preprocess_result/meta.txt') %>% 
  filter(Group == 'Tumour') %>% 
  left_join(follow_up[,-1],  by = 'Pathological_Barcode') %>%
  left_join(riskScore, by = c('Pathological_Barcode' = 'sample')) %>%
  left_join(tn, by = c('Tumor_Sample_Barcode' = '测序')) %>%
  left_join(feat_f, by = 'Pathological_Barcode') %>%
  left_join(feat_f2, by = 'Pathological_Barcode')

meta$Esophgus_Location[meta$Esophgus_Location == '0'] <- 'middle'  # 有一个记录错误，先将其修正为middle
meta[is.na(meta)] <- NA # 缺失记录
meta <- arrange(meta, riskScore) # 按照风险值排序

annotation_df <- data.frame(
  RiskScore = meta$riskScore,
  `I-I` = meta$`I-I`,
  `I-S` = meta$`I-S`,
  `S-S` = meta$`S-S`,
  `T-I` = meta$`T-I`,
  `T-S` = meta$`T-S`,
  Age = meta$Age,
  Sex = factor(ifelse(meta$Sex == 'male', 'Male', 'Female')),
  Nation = factor(ifelse(meta$Nation == '汉', 'Han', 'Wei')),
  Smoking = factor(ifelse(meta$Smoking_Status == 0, 'Non-smoker', 'Smoker')),
  Drinking = factor(ifelse(meta$Drinking_Status == 0, 'Non-drinker', 'Drinker'), levels = c('Non-drinker', 'Drinker')),
  Location = factor(str_to_title(meta$Esophgus_Location), levels = c('Lower', 'Middle-Lower', 'Middle', 'Upper')), 
  Stage = factor(meta$Stage),
  `T stage` = factor(meta$T.stage),
  `N stage` = factor(meta$N.stage),
  Subtype = factor(meta$Subtype),
  check.names = F)

annotation_df <- arrange(annotation_df, RiskScore) # 按照风险值排序

# 注释条颜色
annotation_colors = list(
  RiskScore = colorRamp2(c(min(annotation_df$RiskScore), median(annotation_df$RiskScore), max(annotation_df$RiskScore)), c("#999999", "#DDDDDD", "purple")),
  `I-I` = colorRamp2(c(min(annotation_df$`I-I`), mean(annotation_df$`I-I`), max(annotation_df$`I-I`)), c("#999999", "#DDDDDD", "purple")),
  `I-S` = colorRamp2(c(min(annotation_df$`I-S`), mean(annotation_df$`I-S`), max(annotation_df$`I-S`)), c("#999999", "#DDDDDD", "purple")),
  `S-S` = colorRamp2(c(min(annotation_df$`S-S`), mean(annotation_df$`S-S`), max(annotation_df$`S-S`)), c("#999999", "#DDDDDD", "purple")),
  `T-I` = colorRamp2(c(min(annotation_df$`T-I`), mean(annotation_df$`T-I`), max(annotation_df$`T-I`)), c("#999999", "#DDDDDD", "purple")),
  `T-S` = colorRamp2(c(min(annotation_df$`T-S`), mean(annotation_df$`T-S`), max(annotation_df$`T-S`)), c("#999999", "#DDDDDD", "purple")),
  Age = colorRamp2(c(min(annotation_df$Age), mean(annotation_df$Age), max(annotation_df$Age)), c("#999999", "#DDDDDD", "purple")),
  Sex = c("Female" = "#FF6666", "Male" = "#6666FF"),
  Nation = c("Han" = "#FF6666", "Wei" = "#6666FF"),
  Smoking = c("Non-smoker" = "#FF6666", "Smoker" = "#6666FF"),
  Drinking = c("Non-drinker" = "#FF6666", "Drinker" = "#6666FF"),
  Location = c('Lower' = "#f9a411", 'Middle-Lower' = "#6a75c2", 'Middle' = "#f14166", 'Upper' = '#009f81'),
  Stage = c('Stage I' = "#f9a411", 'Stage II' = "#6a75c2", 'Stage III' = "#f14166", 'Stage IV' = '#009f81'),
  `T stage` = c('T1' = "#f9a411", 'T2' = "#6a75c2", 'T3' = "#f14166", 'T4' = '#009f81'),
  `N stage` = c('N0' = "#f9a411", 'N1' = "#6a75c2", 'N2' = "#f14166", 'N3' = '#009f81'),
  Subtype = c("Subtype1" = "#f9a411", "Subtype2" = "#6a75c2", "Subtype3" = "#f14166")
)

ha = HeatmapAnnotation(df = annotation_df, 
                       col = annotation_colors, 
                       gp = gpar(col = "white"), # 白色边框
                       show_annotation_name = TRUE, 
                       annotation_name_side = "left", 
                       annotation_legend_param = lapply(annotation_colors, function(x) list(ncol = 1, direction = "horizontal")) # 分类变量的legend垂直放，连续变量legend水平放
)

##
heat <- t(meta[,20:37])
ht <- Heatmap(heat, 
              show_row_names = FALSE,                                                                          
              show_column_names = FALSE, 
              show_heatmap_legend = TRUE, 
              cluster_rows = T, 
              cluster_columns = F,
              top_annotation = ha,
              heatmap_legend_param = list(
                direction = "horizontal",   # 图例横向排列
                legend_width = unit(4, "cm"), # 图例宽度
                title = "Feature score" # 图例标题
              ))

pdf('./fig3.lassocox_heatmap.pdf', width = 16, height = 7)
draw(ht, heatmap_legend_side = "bottom")
dev.off()


rm(list = ls())
gc()
