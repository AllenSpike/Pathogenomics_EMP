### 病理无监督分层差异基因/细胞/受配体/转录因子
### 20241129

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
library(GSVA)
library(xCell)
library(ComplexHeatmap)
library(enrichplot)
library(export)
library(ggrepel)
library(ggExtra)
library(pbapply)
library(glmnet)
library(dynamicTreeCut)
library(ggVolcano)
library(gghalves)
library(GseaVis)
library(BulkSignalR)
library(pbmcapply)
library(sva)
library(caret)
library(UpSetR)
library(ggVennDiagram)
library(ggpmisc)
library(fmsb)
library(octad)
library(doParallel)

path_output = './data/patch/xj_4096/'
path_output2 = './data/patch/tcga_4096/'
path_output3 = './data/patch/gdph_supp_4096/'

agg_type = c('mean', 'sd', 'skewness', 'kurtosis')


# 数据汇总 ---------------------------------------------------------------------
# 背景文件
ncbi <- read.delim('./data/seq/resource/Homo_sapiens.gene_info', check.names = F) # 基因文件
genelocate <- read.table(system.file("extdata","genelocate.txt",package="multiOmicsViz"), 
                         header = T, sep = '\t',stringsAsFactors=FALSE, check.names=FALSE) # 基因位置文件
gene_esca <- read.xlsx('./data/seq/resource/esca_related_gene_summary.xlsx') # ESCA相关基因: 来自DisGeNet(https://www.disgenet.org/)
gene_driver <- c('KRT5', 'CDH10', 'LILRB3', 'YEATS2', 'CASP8',
                 'TP53', 'ZNF750', 'NOTCH1', 'FAT1', 'NFE2L2') # ESCC驱动基因（来自2020-CellRes-Whole-genome sequencing of 508 patients identifies key molecular features associated with poor prognosis in esophageal squamous cell carcinoma）

anno_go <- read.gmt('./data/seq/resource/go.v2024.1.Hs.symbols.gmt')
anno_hallmark <- read.gmt('./data/seq/resource/hallmark.v2023.2.Hs.symbols.gmt')
anno_kegg <- read.gmt('./data/seq/resource/kegg_legacy.v2024.1.Hs.symbols.gmt')
anno_reactome <- read.gmt('./data/seq/resource/reactome.v2024.1.Hs.symbols.gmt')

#### XJ临床信息 ####
follow_up <- fread('./data/seq/preprocess_result/follow_up.txt') # 随访表
meta_xj <- read.delim('./data/seq/preprocess_result/meta.txt') # 肿瘤元信息表
meta_t <- filter(meta_xj, Group == 'Tumour') # 肿瘤样本
meta_n <- filter(meta_xj, Group == 'Non-tumour') # 癌旁样本
status_xj <- as.data.frame(follow_up[, c(11,9)]) # 随访信息
colnames(status_xj) <- c('time', 'status')
rownames(status_xj) <- follow_up$Pathological_Barcode

#### XJ病理特征 ####
feat_xj <- t(do.call(rbind, lapply(agg_type, function(x){
  s <- read.xlsx(paste0(path_output,"/feat_agg.xlsx"), rowNames = T, sheet = x)[1:252,] # 选择前252个特征
  rownames(s) <- paste0(rownames(s), "_", x)
  return(s)
}
)))
colnames(feat_xj) <- gsub('GraphAvg_', 'GraphAvg\\.', colnames(feat_xj))
feat_xj <- feat_xj[rownames(status_xj),] # 有随访信息的样本
feat_xj <- apply(feat_xj, 2, function(x){
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

# 标准化
feat_xj_norm <- scale(feat_xj)
feat_xj_norm <- apply(feat_xj_norm, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na


#### XJ分子信息 ####
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

#### TCGA临床信息 ####
clinical <- fread('./data/seq/preprocess_result/TCGA-ESCA/clinical.tsv')
clinical$Overall_Survival_Months <- ifelse(clinical$vital_status == 'Dead',
                                           as.numeric(clinical$days_to_death)/30,
                                           as.numeric(clinical$days_to_last_follow_up)/30)
clinical$Live_Status <- ifelse(clinical$vital_status == 'Alive', 0, 1)
clinical <- clinical[!duplicated(clinical$case_submitter_id), ]

#### TCGA病理特征 ####
feat_tcga <- t(do.call(rbind, lapply(agg_type, function(x){
  s <- read.xlsx(paste0(path_output2,"/feat_agg.xlsx"), rowNames = T, sheet = x)[1:252,] # 选择前252个特征
  rownames(s) <- paste0(rownames(s), "_", x)
  return(s)
}
)))
feat_tcga <- feat_tcga[setdiff(rownames(feat_tcga), c("TCGA-IG-A3I8-01Z-00-DX1.00D0CD6B-06E7-4F3B-BD24-2E18DE926045",
                                                      "TCGA-IG-A3QL-01Z-00-DX1.D7C07E36-F35B-4369-A215-159E0085767C",
                                                      "TCGA-IG-A6QS-01Z-00-DX1.F8F804CA-97FA-4941-BB0A-1E651424E009",
                                                      "TCGA-IG-A97H-01Z-00-DX1.DB5FAE75-4F0D-4D25-BE6F-928CA784EEF7",
                                                      "TCGA-IG-A97I-01Z-00-DX1.3C27B8C0-A346-4B05-98EF-6470C91F57BC",
                                                      "TCGA-L5-A4OM-01Z-00-DX1.09D4AD58-654C-495D-BD03-4973977E50EE",
                                                      "TCGA-LN-A9FO-01Z-00-DX1.E8183919-40B9-4CF2-BC78-00B8C3E6C2C3")),] # 取没有污渍的单纯手术样本(52)
rownames(feat_tcga) <- substr(rownames(feat_tcga), 1, 12)
colnames(feat_tcga) <- gsub('GraphAvg_', 'GraphAvg\\.', colnames(feat_tcga))

clinical_f <- clinical[clinical$case_submitter_id %in% rownames(feat_tcga), ] # 有随访信息的样本
status_tcga <- as.data.frame(clinical_f[, c(159, 160)])
colnames(status_tcga) <- c('time', 'status')
rownames(status_tcga) <- clinical_f$case_submitter_id
status_tcga["TCGA-L5-A43J", "time"] <- 0.01 # 该样本的生存时间为0，会导致后续计算出错，因此将其生存时间设为0.01
feat_tcga <- apply(feat_tcga, 2, function(x){
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

# 标准化
feat_tcga_norm <- scale(feat_tcga)
feat_tcga_norm <- apply(feat_tcga_norm, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na

#### TCGA分子信息 ####
# SNV
load('./data/seq/preprocess_result/TCGA-ESCA_SNP.RData')
laml_tcga <- read.maf(data)

laml_tcga@data$Tumor_Sample_Barcode <- substr(laml_tcga@data$Tumor_Sample_Barcode, 1, 12)
laml_tcga@variants.per.sample$Tumor_Sample_Barcode <- substr(laml_tcga@variants.per.sample$Tumor_Sample_Barcode, 1, 12)
laml_tcga@variant.type.summary$Tumor_Sample_Barcode <- substr(laml_tcga@variant.type.summary$Tumor_Sample_Barcode, 1, 12)
laml_tcga@variant.classification.summary$Tumor_Sample_Barcode <- substr(laml_tcga@variant.classification.summary$Tumor_Sample_Barcode, 1, 12)
laml_tcga@maf.silent$Tumor_Sample_Barcode <- substr(laml_tcga@maf.silent$Tumor_Sample_Barcode, 1, 12)
laml_tcga@clinical.data$Tumor_Sample_Barcode <- substr(laml_tcga@clinical.data$Tumor_Sample_Barcode, 1, 12)

# CNV
cnv_tcga <- read.delim('./data/seq/preprocess_result/TCGA-ESCC_CNV_all/all_data_by_genes.txt')
cnv_tcga$Gene.Symbol <- gsub('\\|chr\\d', '', cnv_tcga$Gene.Symbol)
colnames(cnv_tcga) <- substr(colnames(cnv_tcga), 1, 12)
colnames(cnv_tcga)[4:196] <- gsub('\\.', '-', colnames(cnv_tcga)[4:196])
cnv_tcga <- cnv_tcga[,!duplicated(colnames(cnv_tcga))]

laml.gistic_tcga <- readGistic(gisticAllLesionsFile = './data/seq/preprocess_result/TCGA-ESCC_CNV_all/all_lesions.conf_90.txt',
                               gisticAmpGenesFile = './data/seq/preprocess_result/TCGA-ESCC_CNV_all/amp_genes.conf_90.txt',
                               gisticDelGenesFile = './data/seq/preprocess_result/TCGA-ESCC_CNV_all/del_genes.conf_90.txt',
                               gisticScoresFile = './data/seq/preprocess_result/TCGA-ESCC_CNV_all/scores.gistic')


# RNA
load('./data/seq/preprocess_result/TCGA-ESCA_RNA.RData')
tcga_escc <- data[,grepl('Squamous', as.data.frame(colData(data))$primary_diagnosis)] # 取escc表达谱
clinical_2 <- as.data.frame(colData(tcga_escc)) # 取escc临床信息
exp_tcga_count <- assay(tcga_escc, "stranded_second")[, clinical_2$sample_type == 'Primary Tumor'] # 取原发灶
exp_tcga_tpm <- assay(tcga_escc, "tpm_unstrand")[, clinical_2$sample_type == 'Primary Tumor']
rownames(exp_tcga_tpm) <- rownames(exp_tcga_count) <- rowData(tcga_escc)$gene_name # 修改基因名
colnames(exp_tcga_tpm) <- colnames(exp_tcga_count) <- substr(colnames(exp_tcga_count), 1, 12)

samp <- intersect(rownames(feat_tcga), colnames(exp_tcga_count)) # 51
exp_tcga_count <- exp_tcga_count[,samp]
exp_tcga_tpm <- exp_tcga_tpm[,samp]

clinical_f <- clinical_f[clinical_f$case_submitter_id %in% samp, ] # 50
feat_tcga <- feat_tcga[samp,]
feat_tcga_norm <- feat_tcga_norm[samp,]

#### GDPH临床信息 ####
gdph <- read.xlsx('./data/wsi/gdph/GDPH_临床信息_单纯手术_20241106.xlsx')
gdph <- gdph[!duplicated(gdph$病理号),] # 去重
gdph$`生存状态（修正版）` <- ifelse(gdph$`生存状态（修正版）` == 'alive', 0, 1)
status_gdph <- as.data.frame(gdph[, c(114, 106)])
colnames(status_gdph) <- c('time', 'status')
rownames(status_gdph) <- gdph$病理号

#### GDPH病理特征 ####
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

feat_gdph <- feat_gdph[!duplicated(rownames(feat_gdph)),] # 1935862有两个，去重
samp <- intersect(rownames(feat_gdph), rownames(status_gdph)) # 共有样本
status_gdph <- status_gdph[samp,]
feat_gdph <- feat_gdph[samp,]
# feat_gdph_dl <- feat_gdph_dl[samp,]

# 标准化
feat_gdph_norm <- scale(feat_gdph)
feat_gdph_norm <- apply(feat_gdph_norm, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = T)
  return(x)
}) # 每列均值填充na


# 病理基因组--------------------------------------------------------------------
cor_matrix <- cor(feat_xj, use = "pairwise.complete.obs", method = 'pearson')
high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9)
feature <- setdiff(colnames(feat_xj), colnames(feat_xj)[high_cor_vars])

#### LassoCox风险值 ####
df <- data.frame()
cvs <- list()

set.seed(42)

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

# 选择模型
df_f <- filter(df, pvalue < 0.05 & pvalue_tcga < 0.05 & pvalue_gdph < 0.05)
cv_final <- cvs[['4']][['0.3']]
plot(cv_final)

save(cv_final, file = './supp1_lasso_cox_model.RData')
load('./supp1_lasso_cox.RData')

# 风险值
risk_xj <- as.data.frame(predict(cv_final, newx = feat_xj[,feature], type = "link", s = "lambda.1se"))
risk_gdph <- as.data.frame(predict(cv_final, newx = feat_gdph[,feature], type = "link", s = "lambda.1se"))
risk_tcga <- as.data.frame(predict(cv_final, newx = feat_tcga[,feature], type = "link", s = "lambda.1se"))

colnames(risk_gdph) <- colnames(risk_tcga) <- colnames(risk_xj) <- "riskScore"
risk_cutoff <- median(risk_xj$riskScore) # cutoff

# 模型特征
rs <- coef(cv_final, s = "lambda.1se")
coef <- rs[which(as.numeric(rs) != 0)]
feature <- rownames(rs)[which(as.numeric(rs) != 0)]

##### KM ##### 
# XJ
risk_xj$riskGroup <- ifelse(risk_xj$riskScore > risk_cutoff, "High","Low") # 风险组
status_xj$riskGroup <- risk_xj$riskGroup[match(rownames(status_xj), rownames(risk_xj))]
fit <- survfit(Surv(time, status) ~ riskGroup, data = status_xj)
p1 <- ggsurvplot(fit = fit, 
                 data = status_xj,
                 palette = c('#f14166', '#6a75c2'),
                 pval = T,
                 risk.table = T,
                 surv.median.line = "hv", #添加中位生存曲线
                 title = "XJ", #标题
                 ylab = "Cumulative survival (percentage)", xlab = " Time (Months)", #更改横纵坐标
                 censor.shape = 124,
                 censor.size = 4,
                 conf.int = FALSE
)

tumorDiff <- survdiff(Surv(time, status) ~ riskGroup, data = status_xj)
D1 <- tumorDiff$obs[1]
D2 <- tumorDiff$obs[2]
E1 <- tumorDiff$exp[1]
E2 <- tumorDiff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
SE_lnHR = sqrt(1/E1 + 1/E2)
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95

# tcga
risk_tcga$riskGroup <- ifelse(risk_tcga$riskScore > risk_cutoff, "High","Low") # 风险组
status_tcga$riskGroup <- risk_tcga$riskGroup[match(rownames(status_tcga), rownames(risk_tcga))]
fit <- survfit(Surv(time, status) ~ riskGroup, data = status_tcga)
p2 <- ggsurvplot(fit = fit, 
           data = status_tcga,
           palette = c('#f14166', '#6a75c2'),
           pval = T,
           risk.table = T,
           surv.median.line = "hv", #添加中位生存曲线
           title = "TCGA", #标题
           ylab = "Cumulative survival (percentage)", xlab = " Time (Months)", #更改横纵坐标
           censor.shape = 124,
           censor.size = 4,
           conf.int = FALSE
)

tumorDiff <- survdiff(Surv(time, status) ~ riskGroup, data = status_tcga)
D1 <- tumorDiff$obs[1]
D2 <- tumorDiff$obs[2]
E1 <- tumorDiff$exp[1]
E2 <- tumorDiff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
SE_lnHR = sqrt(1/E1 + 1/E2)
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95

# gdph
risk_gdph$riskGroup <- ifelse(risk_gdph$riskScore > risk_cutoff, "High","Low") # 风险组
status_gdph$riskGroup <- risk_gdph$riskGroup[match(rownames(status_gdph), rownames(risk_gdph))]
fit <- survfit(Surv(time, status) ~ riskGroup, data = status_gdph)
p3 <- ggsurvplot(fit = fit, 
           data = status_gdph,
           palette = c('#f14166', '#6a75c2'),
           pval = T,
           risk.table = T,
           surv.median.line = "hv", #添加中位生存曲线
           title = "GDPH", #标题
           ylab = "Cumulative survival (percentage)", xlab = " Time (Months)", #更改横纵坐标
           censor.shape = 124,
           censor.size = 4,
           conf.int = FALSE
)

tumorDiff <- survdiff(Surv(time, status) ~ riskGroup, data = status_gdph)
D1 <- tumorDiff$obs[1]
D2 <- tumorDiff$obs[2]
E1 <- tumorDiff$exp[1]
E2 <- tumorDiff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
SE_lnHR = sqrt(1/E1 + 1/E2)
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95

splots <- list()
splots[[1]] <- p1
splots[[2]] <- p2
splots[[3]] <- p3
plots <- arrange_ggsurvplots(splots, print = T, ncol = 3, nrow = 1)
ggsave("./supp1_lasso_cox_KM.pdf", plots, width = 15, height = 6)



####### 注释热图 ###### 
tn <- read.xlsx('./data/seq/TN分期.xlsx', sheet = 1) # TNM分期表
heat <- as.data.frame(feat_xj_norm[,feature])
heat$Pathological_Barcode <- as.integer(rownames(heat))
heat2 <- data.frame(
  `I-I` = rowMeans(heat[,grep('I-I', colnames(heat))]),
  `S-S` = rowMeans(heat[,grep('S-S', colnames(heat))]),
  `T-I` = heat[,grep('T-I', colnames(heat))],
  `T-S` = heat[,grep('T-S', colnames(heat))],
  Pathological_Barcode = heat$Pathological_Barcode,
  check.names = F
)

risk_xj$Pathological_Barcode <- as.integer(rownames(risk_xj))
meta_ht <- read.delim('./data/seq/preprocess_result/meta.txt') %>% 
  filter(Group == 'Tumour') %>% 
  left_join(follow_up[,-1],  by = 'Pathological_Barcode') %>%
  left_join(risk_xj, by = c('Pathological_Barcode')) %>%
  left_join(tn, by = c('Tumor_Sample_Barcode' = '测序')) %>%
  left_join(heat, by = 'Pathological_Barcode') %>%
  left_join(heat2, by = 'Pathological_Barcode')

meta_ht$Esophgus_Location[meta_ht$Esophgus_Location == '0'] <- 'middle'  # 有一个记录错误，先将其修正为middle
meta_ht[is.na(meta_ht)] <- NA # 缺失记录
meta_ht <- arrange(meta_ht, riskScore) # 按照风险值排序

annotation_df <- data.frame(
  RiskScore = meta_ht$riskScore,
  `I-I` = meta_ht$`I-I`,
  `S-S` = meta_ht$`S-S`,
  `T-I` = meta_ht$`T-I`,
  `T-S` = meta_ht$`T-S`,
  Age = meta_ht$Age,
  Sex = factor(ifelse(meta_ht$Sex == 'male', 'Male', 'Female')),
  Nation = factor(ifelse(meta_ht$Nation == '汉', 'Han', 'Wei')),
  Smoking = factor(ifelse(meta_ht$Smoking_Status == 0, 'Non-smoker', 'Smoker')),
  Drinking = factor(ifelse(meta_ht$Drinking_Status == 0, 'Non-drinker', 'Drinker'), levels = c('Non-drinker', 'Drinker')),
  Location = factor(str_to_title(meta_ht$Esophgus_Location), levels = c('Lower', 'Middle-Lower', 'Middle', 'Upper')), 
  Stage = factor(meta_ht$Stage),
  `T stage` = factor(meta_ht$T.stage),
  `N stage` = factor(meta_ht$N.stage),
  Subtype = factor(meta_ht$Subtype),
  check.names = F)

annotation_df <- arrange(annotation_df, RiskScore) # 按照风险值排序

# 注释条颜色
annotation_colors = list(
  RiskScore = colorRamp2(c(min(annotation_df$RiskScore), median(annotation_df$RiskScore), max(annotation_df$RiskScore)), c("#999999", "#DDDDDD", "purple")),
  `I-I` = colorRamp2(c(min(annotation_df$`I-I`), mean(annotation_df$`I-I`), max(annotation_df$`I-I`)), c("#999999", "#DDDDDD", "purple")),
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
heat <- t(meta_ht[,21:27])
ht <- Heatmap(heat, 
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 7), row_names_side = "left",
              show_column_names = FALSE, show_row_dend = F,
              show_heatmap_legend = TRUE, 
              cluster_rows = T, 
              cluster_columns = F,
              top_annotation = ha,
              heatmap_legend_param = list(
                direction = "horizontal",   # 图例横向排列
                legend_width = unit(4, "cm"), # 图例宽度
                title = "Feature score" # 图例标题
              ))

draw(ht, heatmap_legend_side = "bottom")


# #### 队列测试 ####
# # top_features_ls <- list()
# # 
# ## xj
# name <- 'XJ'
# feat <- feat_xj
# status <- status_xj
# pca <- prcomp(feat, scale = T, center = T)
# pca_data <- pca$x[, 1:10] # 提取前n个主成分
# # 
# # ## tcga
# # name <- 'TCGA'
# # feat <- feat_tcga
# # status <- status_tcga
# # pca <- prcomp(feat, scale = F, center = T)
# # pca_data <- pca$x[, 1:10] # 提取前n个主成分
# # 
# # ## gdph
# # name <- 'GDPH'
# # feat <- feat_gdph
# # status <- status_gdph
# # pca <- prcomp(feat, scale = T, center = T)
# # pca_data <- pca$x[, 1:10] # 提取前n个主成分
# # 
# # ## XJ + TCGA
# # name <- 'XJ + TCGA'
# # feat <- rbind(feat_xj, feat_tcga)
# # feat <- t(ComBat(t(feat), batch = c(rep(1, nrow(feat_xj)), rep(2, nrow(feat_tcga))))) # 去批次
# # feat_norm <- scale(feat)
# # status <- rbind(status_xj, status_tcga)
# # pca <- prcomp(feat, scale = T, center = T)
# # pca_data <- pca$x[, 1:10] # 提取前n个主成分
# # 
# # ## XJ + GDPH
# # name <- 'XJ + GDPH'
# # feat <- rbind(feat_xj, feat_gdph)
# # feat <- t(ComBat(t(feat), batch = c(rep(1, nrow(feat_xj)), rep(2, nrow(feat_gdph)))))
# # feat_norm <- scale(feat)
# # status <- rbind(status_xj, status_gdph)
# # pca <- prcomp(feat, scale = T, center = T)
# # pca_data <- pca$x[, 1:10]
# # 
# # ## TCGA + GDPH
# # name <- 'TCGA + GDPH'
# # feat <- rbind(feat_tcga, feat_gdph)
# # feat <- t(ComBat(t(feat), batch = c(rep(1, nrow(feat_tcga)), rep(2, nrow(feat_gdph)))))
# # feat_norm <- scale(feat)
# # status <- rbind(status_tcga, status_gdph)
# # pca <- prcomp(feat, scale = T, center = T)
# # pca_data <- pca$x[, 1:10]
# # 
# # ## XJ + TCGA + GDPH
# # name <- 'XJ + TCGA + GDPH'
# # feat <- rbind(feat_xj, feat_tcga, feat_gdph)
# # feat <- t(ComBat(t(feat), batch = c(rep(1, nrow(feat_xj)), rep(2, nrow(feat_tcga)), rep(3, nrow(feat_gdph)))))
# # feat_norm <- scale(feat)
# # status <- rbind(status_xj, status_tcga, status_gdph)
# # pca <- prcomp(feat, scale = T, center = T)
# # pca_data <- pca$x[, 1:10]
# 
# # kmeans
# group1 <- kmeans(pca_data, centers = 2, nstart = 10, iter.max = 10000)$cluster
# table(group1)
# group1 <- ifelse(group1 == 1, 2, 1) # 修正聚类名称
# 
# # pca
# pca_df <- as.data.frame(pca_data)
# pca_df$Cluster <- as.factor(group1)
# pca_df$Data <- ifelse(rownames(pca_df) %in% rownames(feat_xj), 'XJ', 'TCGA')
# summ <- summary(pca)
# xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
# ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
# ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, shape = Data)) +
#   geom_point(size = 1.5) +
#   stat_ellipse(level = 0.95) +
#   scale_color_manual(values = c('1' = '#f14166', '2' = '#6a75c2')) +
#   labs(x = xlab, y = ylab , title = "") +
#   theme_bw()
# 
# # 绘制KM曲线
# sf <- pca_df[,c('Cluster', 'Data')]
# sf$case_submitter_id <- rownames(pca_df)
# 
# status$case_submitter_id <- rownames(status)
# status <- left_join(status, sf, by = "case_submitter_id")
# colnames(status)[1:2] <- c('time', 'status')
# fit <- survfit(Surv(time, status) ~ Cluster, data = status)
# 
# ggsurvplot(fit = fit,
#            data = status,
#            palette = c('#f14166', '#6a75c2'),
#            pval = T,
#            risk.table = T,
#            surv.median.line = "hv", #添加中位生存曲线
#            title = name, #标题
#            ylab = "Cumulative survival (percentage)", xlab = " Time (Months)", #更改横纵坐标
#            censor.shape = 124,
#            censor.size = 4,
#            conf.int = FALSE
# )
# 
# # 提取对前n个主成分贡献最大的前m个特征 
# sorted_contributions <- sort(rowMeans(pca$rotation[,1:10]), decreasing = TRUE) 
# top_features <- sorted_contributions
# top_features_ls[[name]] <- top_features
# 
# # 共识特征 - 贡献排名前25%的特征
# upset(fromList(lapply(top_features_ls, function(x) names(x[x > quantile(x)[4]]))),
#       nsets = 7, set_size.show = FALSE,
#       order.by = "freq")
# 
# Reduce(intersect, lapply(top_features_ls, function(x) names(x[x > quantile(x)[4]]))) # 交集
  

#### 合并XJ/TCGA病理特征无监督聚类 ####
feat <- rbind(feat_xj, feat_tcga)
feat <- t(ComBat(t(feat), batch = c(rep(1, nrow(feat_xj)), rep(2, nrow(feat_tcga))))) # 去批次
feat_norm <- scale(feat)

# pca
pca <- prcomp(feat, scale = T, center = T)
plot(pca$sdev[1:50])
cumsum(pca$sdev)/sum(pca$sdev)
pca_data <- pca$x[, 1:10] # 提取前n个主成分

sorted_contributions <- sort(rowMeans(pca$rotation[,1:10]), decreasing = TRUE) 
top_features <- sorted_contributions[1:50]

# kmeans
group1 <- kmeans(pca_data, centers = 2, nstart = 1000, iter.max = 10000)$cluster
table(group1)
group1 <- ifelse(group1 == 1, 2, 1) # 修正聚类名称

# pca
pca_df <- as.data.frame(pca_data)
pca_df$Cluster <- as.factor(group1)
pca_df$Data <- ifelse(rownames(pca_df) %in% rownames(feat_xj), 'XJ', 'TCGA')

summ <- summary(pca)
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, shape = Data)) +
  geom_point(size = 1.5) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c('1' = '#f14166', '2' = '#6a75c2')) +
  labs(x = xlab, y = ylab , title = "") +
  theme_bw()

# 绘制KM曲线
sf <- pca_df[,c('Cluster', 'Data')]
sf$case_submitter_id <- rownames(pca_df)

status <- rbind(status_xj, status_tcga)
status$case_submitter_id <- rownames(status)
status <- left_join(status, sf, by = "case_submitter_id")
colnames(status)[1:2] <- c('time', 'status')
fit <- survfit(Surv(time, status) ~ Cluster, data = status)

ggsurvplot(fit = fit, 
           data = status,
           palette = c('#f14166', '#6a75c2'),
           pval = T,
           risk.table = T,
           surv.median.line = "hv", #添加中位生存曲线
           title = "XJ + TCGA", #标题
           ylab = "Cumulative survival (percentage)", xlab = " Time (Months)", #更改横纵坐标
           censor.shape = 124,
           censor.size = 4,
           conf.int = FALSE
)

meta <- rbind(clinical_f[,c(1,2)], meta_t[, c(4,5)], use.names = FALSE)
rownames(meta) <- meta$case_submitter_id
meta$Cluster <- as.character(group1[match(meta$case_submitter_id, names(group1))]) # 样本排序
meta$Data <- ifelse(rownames(meta) %in% rownames(feat_xj), 'XJ', 'TCGA')

tumorDiff <- survdiff(Surv(time, status) ~ riskGroup, data = status)
D1 <- tumorDiff$obs[1]
D2 <- tumorDiff$obs[2]
E1 <- tumorDiff$exp[1]
E2 <- tumorDiff$exp[2]
HR <- (D1/D2)/(E1/E2)
HR
SE_lnHR = sqrt(1/E1 + 1/E2)
L = log(HR)
lower <- exp(L - 1.96*SE_lnHR)
upper <- exp(L + 1.96*SE_lnHR)
ci95 <- c(lower=lower, upper=upper)
ci95

#### 无监督组的代表性病理特征 ####
heat <- feat_norm[,names(top_features[1:50])]

df <- as.data.frame(group1)
colnames(df) <- 'Cluster'
df$Cluster <- as.character(df$Cluster)
df$Data <- ifelse(rownames(df) %in% rownames(feat_xj), 'XJ', 'TCGA')

dm <- dist(heat[rownames(df)[df$Cluster == '1'],], method = "euclidean")
hc_1 <- hclust(dm)
dm <- dist(heat[rownames(df)[df$Cluster == '2'],], method = "euclidean")
hc_2 <- hclust(dm)

ord <- c(rownames(df)[df$Cluster == '1'][hc_1$order], 
         rownames(df)[df$Cluster == '2'][hc_2$order]) 
heat <- heat[ord,]
col_annotation <- HeatmapAnnotation(df = df[ord,],
                                    show_annotation_name = F,
                                    col = list(Cluster = c('1' = '#f14166', '2' = '#6a75c2'), 
                                               Data = c('XJ' = '#f9a411', 'TCGA' = '#ae5b94')))
row_colors <- str_split(colnames(heat), '_', simplify = T)[,1]
row_colors <- ifelse(row_colors == 'T-T', '#f14166',
                     ifelse(row_colors == 'I-I', '#6a75c2',
                            ifelse(row_colors == 'S-S', '#f9a411', 
                                   ifelse(row_colors == 'T-I', '#ae5b94',
                                          ifelse(row_colors == 'T-S', '#f5733c',
                                                 ifelse(row_colors == 'I-S', '#b28d6a', NA))))))
ht <- Heatmap(t(heat),
              name = 'Feature score',
              col = circlize::colorRamp2(c(-2, 0, 2), c("#5f176c", "#009c9a","#f3e700")),
              top_annotation = col_annotation,
              heatmap_legend_param = list(
                direction = "horizontal",   # 图例横向排列
                legend_width = unit(4, "cm"), # 图例宽度
                title = "Feature score" # 图例标题
              ),
              use_raster = F,
              show_row_names = T,
              show_column_names = F, 
              show_row_dend = T, 
              cluster_rows = T,
              cluster_columns = F,  
              show_column_dend = F,
              row_names_side = 'left', 
              row_dend_side = 'right', 
              row_names_gp = gpar(fontsize = 6, col = row_colors))

lgd <- Legend(labels = c('T-T', 'I-I', 'S-S', 'T-I', 'T-S', 'I-S'),
              title = 'Interaction',
              title_position = 'topleft',
              legend_gp = gpar(fill = c('#f14166', '#6a75c2', '#f9a411', '#ae5b94', '#f5733c', '#b28d6a')),
              gap = unit(2, "cm"),
              ncol = 1)

draw(ht, heatmap_legend_side = "bottom", annotation_legend_list = list(lgd)) 


#### 无监督聚类差异SNV和CNV ####
meta <- as.data.frame(meta)
meta$sample <- ifelse(meta$Data == 'XJ', meta$case_id, meta$case_submitter_id)
meta$Tumor_Sample_Barcode <- meta$sample

# SNV
laml_tcga@data <- filter(laml_tcga@data, Tumor_Sample_Barcode %in% meta$sample)
laml_tcga@variants.per.sample <- filter(laml_tcga@variants.per.sample, Tumor_Sample_Barcode %in% meta$sample)
laml_tcga@variant.type.summary <- filter(laml_tcga@variant.type.summary, Tumor_Sample_Barcode %in% meta$sample)
laml_tcga@variant.classification.summary <- filter(laml_tcga@variant.classification.summary, Tumor_Sample_Barcode %in% meta$sample)
laml_tcga@maf.silent <- filter(laml_tcga@maf.silent, Tumor_Sample_Barcode %in% meta$sample)
laml_tcga@clinical.data <- filter(laml_tcga@clinical.data, Tumor_Sample_Barcode %in% meta$sample)

laml_combined <- laml

laml_combined@data <- data.table(rbind(as.data.frame(laml@data)[,intersect(colnames(laml_tcga@data), colnames(laml@data))], 
                                       as.data.frame(laml_tcga@data)[,intersect(colnames(laml_tcga@data), colnames(laml@data))]))
laml_combined@variants.per.sample <- data.table(rbind(as.data.frame(laml@variants.per.sample)[,intersect(colnames(laml_tcga@variants.per.sample), colnames(laml@variants.per.sample))], 
                                                      as.data.frame(laml_tcga@variants.per.sample)[,intersect(colnames(laml_tcga@variants.per.sample), colnames(laml@variants.per.sample))]))
laml_combined@variant.type.summary <- data.table(rbind(as.data.frame(laml@variant.type.summary)[,intersect(colnames(laml_tcga@variant.type.summary), colnames(laml@variant.type.summary))], 
                                                       as.data.frame(laml_tcga@variant.type.summary)[,intersect(colnames(laml_tcga@variant.type.summary), colnames(laml@variant.type.summary))]))
laml_combined@variant.classification.summary <- data.table(rbind(as.data.frame(laml@variant.classification.summary)[,intersect(colnames(laml_tcga@variant.classification.summary), colnames(laml@variant.classification.summary))], 
                                                                 as.data.frame(laml_tcga@variant.classification.summary)[,intersect(colnames(laml_tcga@variant.classification.summary), colnames(laml@variant.classification.summary))]))
laml_combined@gene.summary <- data.table(rbind(as.data.frame(laml@gene.summary)[,intersect(colnames(laml_tcga@gene.summary), colnames(laml@gene.summary))], 
                                               as.data.frame(laml_tcga@gene.summary)[,intersect(colnames(laml_tcga@gene.summary), colnames(laml@gene.summary))]))
laml_combined@summary <- data.table(rbind(as.data.frame(laml@summary)[,intersect(colnames(laml_tcga@summary), colnames(laml@summary))], 
                                          as.data.frame(laml_tcga@summary)[,intersect(colnames(laml_tcga@summary), colnames(laml@summary))]))
laml_combined@maf.silent <- data.table(rbind(as.data.frame(laml@maf.silent)[,intersect(colnames(laml_tcga@maf.silent), colnames(laml@maf.silent))], 
                                             as.data.frame(laml_tcga@maf.silent)[,intersect(colnames(laml_tcga@maf.silent), colnames(laml@maf.silent))]))

laml_combined@clinical.data <- data.table(meta)
laml_combined_1 <- subsetMaf(laml_combined, clinQuery = "Cluster == 1")
laml_combined_2 <- subsetMaf(laml_combined, clinQuery = "Cluster == 2")

res_mut <- mafCompare(m1 = laml_combined_1, 
                      m2 = laml_combined_2, 
                      m1Name = 'Cluster1', 
                      m2Name = 'Cluster2',
                      minMut = 3)

genes_mut <- filter(res_mut$results, pval < 0.05)$Hugo_Symbol

####
## Integrated cohort of esophageal squamous cell cancer reveals genomic features underlying clinical characteristics
## TTN/ZNF429：由于巨大编码长度，突变率高并不代表其具有生物学意义
## PKHD1L1：在ESCC年轻患者中突变率较高；在甲状腺癌和肺腺癌中突变具有预后价值（PKHD1L1 blocks the malignant behavior of lung adenocarcinoma cells and restricts tumor growth by regulating CBX7）
## EHBP1：参与内吞囊泡运输
## LRRC56：在转移性样本中具有异常表达，并且在乳腺癌中可通过调控RhoA/ROKs通路，影响EMT相关蛋白的表达（Leucine-rich repeat-containing 56 promotes breast cancer progression via modulation of the RhoA/ROCKs signaling axis）
####

laml_combined@clinical.data$Cluster <- paste0('Cluster', laml_combined@clinical.data$Cluster)
laml_combined@clinical.data$Cluster <- as.factor(laml_combined@clinical.data$Cluster)
fabcolors <- c('#f14166', '#6a75c2')
names(fabcolors) <- c('Cluster1', 'Cluster2')
fabcolors <- list(Cluster = fabcolors)

oncoplot(maf = laml_combined, 
         genes = genes_mut,
         clinicalFeatures = 'Cluster', 
         annotationColor = fabcolors,
         bgCol = 'gray90',
         draw_titv = TRUE,
         sortByAnnotation = TRUE,
         showTitle = FALSE,
         showPct = FALSE,
         sortByMutation = FALSE)

tmbScore <- tmb(maf = laml_combined)
meta$TMB <- tmbScore$total_perMB_log[match(meta$Tumor_Sample_Barcode, tmbScore$Tumor_Sample_Barcode)]

ggplot(meta, aes(x = Cluster, y = TMB, fill = Cluster)) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width = 0.2, linewidth= 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 1.25, color = "white") +
  theme_bw() +
  ggpubr::geom_signif(comparisons = list(c(1,2))) +
  scale_fill_manual(values = c("1" = '#f14166', "2" = '#6a75c2'))

# CNV
cnv <- laml.gistic@cnv.summary
cnv$FGA <- cnv$total/4207

cnv_tcga <- laml.gistic_tcga@cnv.summary
cnv_tcga$FGA <- cnv_tcga$total/10539
cnv_tcga$Tumor_Sample_Barcode <- substr(cnv_tcga$Tumor_Sample_Barcode, 1, 12)
cnv_tcga <- cnv_tcga[!duplicated(cnv_tcga$Tumor_Sample_Barcode),]

cnv_combined <- rbind(cnv, cnv_tcga)
meta$FGA <- cnv_combined$FGA[match(meta$sample, cnv_combined$Tumor_Sample_Barcode)]

ggplot(meta, aes(x = Cluster, y = FGA, fill = Cluster)) +
  geom_half_violin(side = "r", color = NA, alpha = 0.35) +
  geom_half_boxplot(side = "r", errorbar.draw = T, width = 0.2, linewidth= 0.5) +
  geom_half_point_panel(side = "l", shape = 21, size = 1.25, color = "white") +
  theme_bw() +
  ggpubr::geom_signif(comparisons = list(c(1,2))) +
  scale_fill_manual(values = c("1" = '#f14166', "2" = '#6a75c2'))


#### 无监督聚类差异表达基因 ####
gene <- intersect(rownames(exp_count_t), rownames(exp_tcga_count))
exp <- as.matrix(cbind(exp_count_t[gene,], exp_tcga_count[gene,]))
exp <- ComBat_seq(exp, batch = c(rep(1, ncol(exp_count_t)), rep(2, ncol(exp_tcga_count))))
exp <- round(exp)

meta <- meta[match(colnames(exp), meta$sample),]

grouplist <- factor(ifelse(meta$Cluster == 1, 'case', 'control'), levels = c('case', 'control')) 
grouplist <- relevel(grouplist, ref = 'control') # 选择group2为参考组(预后较好)
colData <- data.frame(row.names = meta$sample, grouplist = grouplist)
exp <- exp[,meta$sample]

dds <- DESeqDataSetFromMatrix(countData = exp,
                              colData = colData,
                              design = ~ grouplist)
dds <- dds[rowSums(counts(dds)) > 1 ,]
dds <- DESeq(dds)
res <- data.frame(results(dds))
res$gene <- rownames(res)
res <- add_regulate(res, log2FC_name = "log2FoldChange", fdr_name = "padj", log2FC = 1, fdr = 0.5)
res$regulate <- ifelse(res$regulate == 'Up', 'Cluster1', 
                       ifelse(res$regulate == 'Down', 'Cluster2', res$regulate))
res <- arrange(res, desc(log2FoldChange))

save(res, file = './signature.RData')

# Volcano
ggvolcano(res,
          x = "log2FoldChange",
          y = "padj",
          label = "gene",
          FDR_cut = 0.5, 
          log2FC_cut = 1,
          fills = c("#f14166", "#6a75c2", "#999999"),
          colors = c("#f14166", "#6a75c2", "#999999"),
          custom_label = c(head(res$gene, 10), tail(res$gene, 10)),
          output = FALSE) +
  theme(
    plot.title = element_blank(),   
    plot.caption = element_blank(),  
    legend.position = "right"
  ) +
  labs(title = "", x = "log2FC", y = "-log10(p-value)") +
  guides(
    fill = guide_legend(title = "Cluster"),
    color = guide_legend(title = "Cluster")
  )

# ORA
res_f <- filter(res, pvalue < 0.05 & log2FoldChange > 1)
gene <- res_f$gene
res_ora <- list(
  go = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_go),
  hallmark = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_hallmark)
)

df1 <- res_ora$go@result
df1 <- arrange(df1, desc(Count))
df1 <- df1[1:10,]
df1$v <- -log10(df1$pvalue)
df1 <- arrange(df1, desc(v))

res_f <- filter(res, pvalue < 0.05 & log2FoldChange < -1)
gene <- res_f$gene
res_ora <- list(
  go = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_go),
  hallmark = enricher(gene, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_hallmark)
)

df2 <- res_ora$go@result
df2 <- arrange(df2, desc(Count))
df2 <- df2[1:10,]
df2$v <- -log10(df2$pvalue)
df2 <- arrange(df2, desc(v))

df <- rbind(df1, df2)
df$type <- c(rep('Cluster1', nrow(df1)), rep('Cluster2', nrow(df2)))

df$Description <- str_to_title(gsub('_', ' ', gsub('^GO\\w\\w_', '', df$Description)))
df$Description <- factor(factor(df$Description), levels = rev(df$Description))

ggplot(df, aes(x = Description, y = v, size = v, fill = type)) +
  geom_point(shape = 21, stroke = 0.5) +  
  scale_fill_manual(values = c("Cluster1" = '#f14166', "Cluster2" = '#6a75c2'), name = 'Cluster') + 
  scale_size(range = c(5, 10), name = "-log10(p-value)") +
  labs(x = "GO term", y = "-log10(p-value)") +
  theme_bw() +
  coord_flip()

# GESA
val <- res$log2FoldChange
names(val) <- res$gene

res_GSEA <- list(
  hallmark = GSEA(geneList = val, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_hallmark),
  go = GSEA(geneList = val, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_go)
)

kf <- res_GSEA$hallmark
kf@result$Description <- str_to_title(gsub('_', ' ', gsub('HALLMARK_', '', kf@result$Description)))

gseaNb(object = kf,
       curveCol = jjAnno::useMyCol('paired', 7),
       subPlot = 2,
       termWidth = 35,
       addGene = c('FN1', 'SNAI2', 'VIM'),
       geneSetID = c('HALLMARK_P53_PATHWAY', 
                     'HALLMARK_TGF_BETA_SIGNALING', 
                     'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 
                     'HALLMARK_HYPOXIA',
                     'HALLMARK_INFLAMMATORY_RESPONSE', 
                     'HALLMARK_IL6_JAK_STAT3_SIGNALING', 
                     'HALLMARK_ALLOGRAFT_REJECTION'))


#### 无监督聚类差异TME ####
library(xCell)
library(MCPcounter)
library(immunedeconv)

# # mcpcounter所需文件
# genes <- read.delim('./data/seq/resource/genes.txt',
#                     sep="\t", stringsAsFactors=FALSE, header=TRUE, 
#                     colClasses="character",check.names=FALSE) # 标志基因
# probesets <- read.delim('./data/seq/resource/probesets.txt',
#                         sep="\t", stringsAsFactors=FALSE, header=TRUE, 
#                         colClasses="character",check.names=FALSE) # 探针名

# cibersort所需文件
source("/home/ljc/0_project/0_ESCC/data/seq/resource/CIBERSORT.R")
sig <- read.table("/home/ljc/0_project/0_ESCC/data/seq/resource/LM22.txt", ,row.names=1, header=T, sep="\t",check.names=F)

##
gene <- intersect(rownames(exp_tpm_t), rownames(exp_tcga_tpm))
exp_tpm <- as.matrix(cbind(exp_tpm_t[gene,], exp_tcga_tpm[gene,]))
exp_tpm <- ComBat(exp_tpm, batch = c(rep(1, ncol(exp_tpm_t)), rep(2, ncol(exp_tcga_tpm)))) # 去批次

# tme_mcp <- MCPcounter.estimate(exp_tpm, 
#                                featuresType = 'HUGO_symbols', 
#                                probesets = probesets, 
#                                genes = genes) # mcpcounter 

tme_xcell <- xCellAnalysis(exp_tpm, 
                           rnaseq = TRUE, 
                           scale = FALSE,
                           parallel.sz = 2, 
                           parallel.type = "SOCK") # xcell

tme_cibersort <- CIBERSORT(mixture_file = exp_tpm, 
                           sig_matrix = sig, 
                           perm = 100) # cibersort

## boxplot
cibersort_tme <- t(tme_cibersort[,1:22])
cibersort_tme <-reshape2::melt(cibersort_tme)
cibersort_tme$cluster <- meta$Cluster[match(cibersort_tme$Var2, meta$sample)]

ggplot(cibersort_tme, aes(x = Var1, y = value, fill = cluster)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c('1' = '#f14166', '2' = '#6a75c2')) + 
  labs(x = 'Cell', y = 'Proportion', fill = 'Cluster')


## heatmap
x <- sort(apply(tme_xcell, 1, function(x){
  a <- x[filter(meta, Cluster == 1)$sample]
  b <- x[filter(meta, Cluster == 2)$sample]
  logfc <- log2(mean(a) / mean(b))
  return(logfc)
}))

heat <- tme_xcell[c(names(head(x, 10)), names(tail(x, 10))), 
                  meta$sample]

dm <- dist(t(heat[,meta$sample[meta$Cluster == '1']]), method = "euclidean")
hc_1 <- hclust(dm)
dm <- dist(t(heat[,meta$sample[meta$Cluster == '2']]), method = "euclidean")
hc_2 <- hclust(dm)
ord <- c(meta$sample[meta$Cluster == '1'][hc_1$order],
         meta$sample[meta$Cluster == '2'][hc_2$order])
df <- data.frame(Cluster = meta$Cluster, Data = meta$Data, row.names = meta$sample)
col_annotation <- HeatmapAnnotation(df = df[ord,],
                                    show_annotation_name = F,
                                    col = list(Cluster = c('1' = '#f14166', '2' = '#6a75c2'), 
                                               Data = c('XJ' = '#f9a411', 'TCGA' = '#ae5b94')))
heat <- heat[,ord]
ht <- Heatmap(t(scale(t(heat))),
              col = circlize::colorRamp2(c(-2, 0, 2), c("#377EB8","#F0F0F0","#E41A1C")),
              # col = circlize::colorRamp2(c(0, 1.5, 3), c("#377EB8","#F0F0F0","#E41A1C")),
              top_annotation = col_annotation,
              heatmap_legend_param = list(
                direction = "horizontal",
                legend_width = unit(4, "cm"),
                title = "TME score"
              ),
              use_raster = T,
              show_row_dend = F, 
              show_column_dend = F, 
              show_column_names = F, 
              show_row_names = T, 
              row_names_gp = gpar(fontsize = 8),
              border = T,
              column_split = c(rep('Cluster1', ncol(heat[,meta$sample[meta$Cluster == '1']])),
                               rep('Cluster2', ncol(heat[,meta$sample[meta$Cluster == '2']]))),
              row_split = c(rep(1, 10),
                            rep(2, 10)),
              cluster_columns = T,
              cluster_rows = T)

draw(ht, heatmap_legend_side = "bottom")


#### 无监督聚类差异受配体 ####
# # diff mode
# exp.comp <- as.BSRDataModelComp(prepareDataset(exp[,meta$sample]))
# colA <- which(meta$Cluster == 1)
# colB <- which(meta$Cluster == 2)
# stats <- res[,c(5, 2)]
# colnames(stats) <- c('pval', 'logFC')
# 
# gene <- intersect(rownames(exp.comp@ncounts), rownames(stats))
# 
# stats <- stats[gene,]
# stats$expr <- rowMeans(exp[gene, colA])
# bsrcc <- defineClusterComp(exp.comp, colA, colB, stats)
# exp.comp <- addClusterComp(exp.comp, bsrcc, "random.example")
# bsrinf.comp <- initialInference(exp.comp,
#                                 max.pval = 1,
#                                 min.logFC = 0.1,
#                                 min.t.logFC = 0.1,
#                                 min.pw.size = 1,
#                                 neg.receptors = TRUE,
#                                 pos.targets = FALSE,
#                                 neg.targets = FALSE,
#                                 use.full.network = TRUE,
#                                 with.complex = TRUE,
#                                 "random.example")
# bsrinf.comp.redBP <- reduceToBestPathway(bsrinf.comp)
# LRinter.dataframe <- LRinter(bsrinf.comp.redBP)
# LRinter.dataframe$LR <- paste0(LRinter.dataframe$L, '_', LRinter.dataframe$R)
# bsrsig.redBP <- getLRGeneSignatures(bsrinf.comp.redBP, qval.thres = 1)
# scoresLR <- scoreLRGeneSignatures(exp.comp,
#                                   bsrsig.redBP,
#                                   LR.weight = 0.5,
#                                   name.by.pathway = FALSE,
#                                   robust = FALSE,
#                                   abs.z.score = FALSE)
# rownames(scoresLR) <- gsub(' / ', '_', gsub('[{}]', '', rownames(scoresLR)))
# 
# scoresLR_abs <- scoreLRGeneSignatures(exp.comp,
#                                       bsrsig.redBP,
#                                       name.by.pathway = FALSE,
#                                       robust = FALSE,
#                                       abs.z.score = TRUE)
# rownames(scoresLR_abs) <- gsub(' / ', '_', gsub('[{}]', '', rownames(scoresLR_abs)))
# 
# logfc <- sort(apply(scoresLR_abs, 1, function(x){
#   a <- x[filter(meta, Cluster == 1)$sample]
#   b <- x[filter(meta, Cluster == 2)$sample]
#   logfc <- log2(mean(a) / mean(b))
#   return(logfc)
# }), decreasing = T)
# 
# LRinter.dataframe$logFC <- logfc[match(LRinter.dataframe$LR, names(logfc))]
# LRinter.dataframe_f <- filter(LRinter.dataframe, LR.pval < 0.05) %>% arrange(rank)
# scoresLR_f <- scoresLR[intersect(rownames(scoresLR), LRinter.dataframe_f$LR),]
# scoresLR_abs_f <- scoresLR_abs[intersect(rownames(scoresLR_abs), LRinter.dataframe_f$LR),]

# abs mode
n.proc <- 2
cl <- makeCluster(n.proc)
registerDoParallel(cl)
exp.abs <- prepareDataset(exp[,meta$sample])
exp.abs <- learnParameters(exp.abs,
                           plot.folder = "/home/ljc/0_project/0_ESCC/man",
                           filename = "sdc") # 很耗时
save(exp.abs, file = './exp.abs.RData')
stopCluster(cl)

bsrinf <- initialInference(exp.abs,
                           min.cor = 0.1,
                           rank.p = 0.55,
                           use.full.network = TRUE,
                           with.complex = FALSE)
bsrinf.redBP <- reduceToBestPathway(bsrinf)
LRinter.dataframe <- LRinter(bsrinf.redBP)
LRinter.dataframe$LR <- paste0(LRinter.dataframe$L, '_', LRinter.dataframe$R)
bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres=1)
scoresLR <- scoreLRGeneSignatures(exp.abs,
                                  bsrsig.redBP,
                                  name.by.pathway = FALSE,
                                  robust = FALSE,
                                  abs.z.score = FALSE)
rownames(scoresLR) <- gsub(' / ', '_', gsub('[{}]', '', rownames(scoresLR)))
scoresLR_abs <- scoreLRGeneSignatures(exp.abs,
                                      bsrsig.redBP,
                                      name.by.pathway = FALSE,
                                      robust = FALSE,
                                      abs.z.score = TRUE)
rownames(scoresLR_abs) <- gsub(' / ', '_', gsub('[{}]', '', rownames(scoresLR_abs)))

logfc <- sort(apply(scoresLR_abs, 1, function(x){
  a <- x[filter(meta, Cluster == 1)$sample]
  b <- x[filter(meta, Cluster == 2)$sample]
  logfc <- log2(mean(a) / mean(b))
  return(logfc)
}), decreasing = T)

LRinter.dataframe$logFC <- logfc[match(LRinter.dataframe$LR, names(logfc))]
LRinter.dataframe_f <- filter(LRinter.dataframe, pval < 0.05) %>% arrange(rank)

scoresLR_f <- scoresLR[intersect(rownames(scoresLR), LRinter.dataframe_f$LR),]
scoresLR_abs_f <- scoresLR_abs[intersect(rownames(scoresLR_abs), LRinter.dataframe_f$LR),]

# heatmap
LRinter.dataframe <- arrange(LRinter.dataframe, desc(logFC))
heat <- scoresLR_abs[c(head(filter(LRinter.dataframe_f, LR %in% rownames(scoresLR))$LR, 10),
                       tail(filter(LRinter.dataframe_f, LR %in% rownames(scoresLR))$LR, 10)),]
# heat <- scoresLR[intersect(rownames(scoresLR), LRinter.dataframe_f$LR),]
heat <- heat[,meta$sample]

dm <- dist(t(heat[,meta$sample[meta$Cluster == '1']]), method = "euclidean")
hc_1 <- hclust(dm)
dm <- dist(t(heat[,meta$sample[meta$Cluster == '2']]), method = "euclidean")
hc_2 <- hclust(dm)
ord <- c(meta$sample[meta$Cluster == '1'][hc_1$order],
         meta$sample[meta$Cluster == '2'][hc_2$order])
df <- data.frame(Cluster = meta$Cluster, Data = meta$Data, row.names = meta$sample)
col_annotation <- HeatmapAnnotation(df = df[ord,],
                                    show_annotation_name = F,
                                    col = list(Cluster = c('1' = '#f14166', '2' = '#6a75c2'), 
                                               Data = c('XJ' = '#f9a411', 'TCGA' = '#ae5b94')))
heat <- heat[,ord]

ht <- Heatmap(t(scale(t(heat))), 
              col = circlize::colorRamp2(c(-3, 0, 3), c("#377EB8","#F0F0F0","#E41A1C")),
              # col = circlize::colorRamp2(c(0, 1.5, 3), c("#377EB8","#F0F0F0","#E41A1C")),
              top_annotation = col_annotation,
              heatmap_legend_param = list(
                direction = "horizontal",
                legend_width = unit(4, "cm"),
                title = "Ligand-Receptor score"
              ),
              use_raster = T,
              show_row_dend = F, 
              show_column_dend = F, 
              show_column_names = F, 
              show_row_names = T, 
              row_names_gp = gpar(fontsize = 8),
              border = T,
              column_split = c(rep('Cluster1', ncol(heat[,meta$sample[meta$Cluster == '1']])),
                               rep('Cluster2', ncol(heat[,meta$sample[meta$Cluster == '2']]))),
              row_split = c(rep(1, 10),
                            rep(2, 10)),
              cluster_columns = T,
              cluster_rows = T)

draw(ht, heatmap_legend_side = "bottom")

# HR
meta$time <- status$time[match(meta$case_submitter_id, status$case_submitter_id)]
meta$status <- status$status[match(meta$case_submitter_id, status$case_submitter_id)]

hr <- apply(scoresLR_f, 1, function(x, clin){
  fit <- coxph(Surv(time, status) ~ x, data = meta)
  fitSummary <- summary(fit)
  res <- c(fitSummary$conf.int[,"exp(coef)"], 
           fitSummary$coefficients[,"Pr(>|z|)"], 
           fitSummary$conf.int[,"lower .95"], 
           fitSummary$conf.int[,"upper .95"],  
           fitSummary$concordance[["C"]], 
           fitSummary$concordance[["se(C)"]], 
           fitSummary$concordance[["C"]] - 1.96 * fitSummary$concordance[["se(C)"]], 
           fitSummary$concordance[["C"]] + 1.96 * fitSummary$concordance[["se(C)"]]
  )
  names(res) <- c('HR', 'HR_pvalue', 'HR_L95CI', 'HR_H95CI',
                  'C-index', 'C-index_SE', 'C-index_L95CI', 'C-index_H95CI')
  return(res)
  return(res)
}, clin = follow_up) 

hr <- as.data.frame(t(hr))
hr$coef <- log(hr$HR)
hr <- arrange(hr, desc(HR))
hr$pw.name <- LRinter.dataframe$pw.name[match(rownames(hr), LRinter.dataframe$LR)]
hr$pw.name <- factor(hr$pw.name, levels = rev(unique(hr$pw.name)))
hr$varname <- rownames(hr)
hr$varname <- factor(hr$varname, levels = rev(hr$varname))

hr_f <- rbind(head(hr, 10), tail(hr, 10))
hr_f$Type <- ifelse(hr_f$HR > 1, 'Harmful', 'Protective')

ggplot(hr_f, aes(y = varname, x = HR, colour = Type)) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmax = HR_H95CI, xmin = HR_L95CI), size= 0.7, height = 0.3) +
  geom_vline(aes(xintercept = 1), color="gray", linetype="dashed", size = 0.7)+
  ylab('Ligand-Receptor pair')+
  xlab('Hazard ratio (HR)')+
  ggthemes::theme_few()+
  theme(axis.text.y = element_text(size = 9, color = "black"))+
  theme(axis.text.x = element_text(size = 9, color = "black"))+
  theme(title=element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.y = element_text(size = 7, 
                                   color = hr_f$color))

df <- group_by(hr, pw.name) %>% 
  summarise(HR = mean(HR)) %>% 
  arrange(desc(HR)) %>% 
  dplyr::select(pw.name)
hr$pw.name <- factor(hr$pw.name, levels = rev(unique(df$pw.name)))

hr_f1 <- filter(hr, 
                pw.name %in% c(head(df$pw.name, 20), 
                               tail(df$pw.name, 20)))

ggplot(hr_f1, aes(x = HR, y = pw.name)) +
  geom_segment(aes(x = 0, xend = HR, y = pw.name, yend = pw.name), 
               color = "gray50", linewidth = 0.5) + 
  geom_point(aes(color = HR), size = 3) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 1) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +  
  labs(
    x = "HR",
    y = "Pathway"
  ) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())

# logfc dotplot
LRinter.dataframe_f <- arrange(LRinter.dataframe_f, desc(logFC))
LRinter.dataframe_f$pw.name[LRinter.dataframe_f$pw.name == "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)"] <- 'Regulation of IGF transport and uptake by IGFBPs'

df <- group_by(LRinter.dataframe_f, pw.name) %>% summarise(logFC = mean(logFC)) %>% arrange(desc(logFC)) %>% dplyr::select(pw.name)
LRinter.dataframe_f$pw.name <- factor(LRinter.dataframe_f$pw.name, levels = rev(unique(df$pw.name)))
LRinter.dataframe_f <- LRinter.dataframe_f[complete.cases(LRinter.dataframe_f),]

LRinter.dataframe_f1 <- filter(LRinter.dataframe_f, 
                               pw.name %in% c(head(df$pw.name, 20), 
                                              tail(df$pw.name, 20)))
LRinter.dataframe_f1$LogFC <- LRinter.dataframe_f1$logFC

ggplot(LRinter.dataframe_f1, aes(x = LogFC, y = pw.name)) +
  geom_segment(aes(x = 0, xend = LogFC, y = pw.name, yend = pw.name), 
               color = "gray50", linewidth = 0.5) +  # 横线
  geom_point(aes(color = logFC), size = 3) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  
  labs(
    x = "LogFC",
    y = "Pathway"
  ) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())


#### 无监督聚类差异敏感药物 ####
# 缓存文件路径：/home/ljc/.cache/R/ExperimentHub
load('./3_signature.RData')

# sRGES
colnames(res)[7] <- 'Symbol'
sRGES <- runsRGES(res, max_gene_size=100, permutations=1000, output=TRUE, 
                  outputFolder = './octad') # 药物列表
octadDrugEnrichment(sRGES = sRGES, 
                    outputRank = TRUE,
                    target_type = 'chembl_targets', 
                    outputFolder = './octad') # 药物靶点富集
# /home/ljc/0_project/0_ESCC/octad/enrichFolder/chembl_targets/enriched_chembl_targets.csv

# 逆转疾病特征热图
library(cmapR)

info_sig <- fread('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_sig_info.txt')  %>% filter(pert_type %in% c('trt_cp'))
info_gene <- fread('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_gene_info.txt')
res_lincs <- read.csv('./octad/all__lincs_score.csv')

cid <- filter(res_lincs, pert_iname %in% sRGES_f$pert_iname)$sig_id
rid <- filter(info_gene, pr_gene_symbol %in% res$Symbol)$pr_gene_id

gtx <- parse_gctx('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx', 
                  cid = cid, 
                  rid = rid)
mat_pert <- gtx@mat

# 靶点富集图
eh_dataframe <- data.frame(title = c("cmpd_sets_ChemCluster", 
                                     "cmpd_sets_chembl_targets", "cmpd_sets_mesh"))
row.names(eh_dataframe) <- c("EH7266", "EH7267", "EH7268")
random_gsea_score <- get_ExperimentHub_data("EH7275")
eh_dataframe$object <- row.names(eh_dataframe)
cmpd_sets <- get_ExperimentHub_data((eh_dataframe[eh_dataframe$title == 
                                                    paste0("cmpd_sets_", 'chembl_targets'), "object"]))
cmpdSets <- cmpd_sets$cmpd.sets
names(cmpdSets) <- cmpd_sets$cmpd.set.names
sRGES$rank <- rank(sRGES$sRGES)

top_target <- 'CACNA1C'
cmpdSets[[top_target]]
target_drugs_score <- sRGES$rank[sRGES$pert_iname %in% cmpdSets[[top_target]]]

limma::barcodeplot(sRGES$sRGES, target_drugs_score, 
                   main = top_target, xlab = "sRGES")

# cellline
symbol2ensemble <- bitr(rownames(exp), 
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Hs.eg.db)
exp_ensemble <- exp[symbol2ensemble$SYMBOL,]
rownames(exp_ensemble) <- symbol2ensemble$ENSEMBL[match(rownames(exp_ensemble), 
                                                        symbol2ensemble$SYMBOL)]

cell_line_cluster1 <- computeCellLine(expSet = exp_ensemble, 
                                      case_id = filter(meta, Cluster == 1)$sample,
                                      source = 'expSet', 
                                      outputFolder =  './octad/',
                                      output = TRUE,
                                      LINCS_overlaps = TRUE)
cell_line_cluster1

# sRGEs和IC50的相关系数
source('./script/0_new/topLineEval_crz.R')
corr_srges_ic50 <- topLineEval_crz(topline = 'EFO27', mysRGES = sRGES, outputFolder = './octad')
View(corr_srges_ic50$ic50df)
corr_srges_ic50[["aucgraph"]]

# 细胞系敲除靶点
info_sig <- fread('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_sig_info.txt') %>% filter(pert_type %in% c('trt_sh', 'trt_xpr'))
info_gene <- fread('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_gene_info.txt')

unique(info_sig$pert_iname) # 4371个敲除基因
table(info_sig$cell_id) # 20个细胞系
table(info_sig$pert_type) # 1种处理类型

load('./data/seq/preprocess_result/lincs_knockdown_escc_pr.RData')

gene_knock <- sort(lincs[, 'CACNA1C'], decreasing = T)

res_GSEA_knock <- list(
  hallmark = GSEA(geneList = gene_knock, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_hallmark),
  go = GSEA(geneList = gene_knock, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_go),
  reatome = GSEA(geneList = gene_knock, pvalueCutoff = 1, pAdjustMethod = 'BH', TERM2GENE = anno_reactome)
)


#### 风险值和通路/TME/受配体的相关系数  ####
meta$risk <- rbind(risk_xj[,1:2], risk_tcga)[match(meta$case_submitter_id, rownames(rbind(risk_xj[,1:2], risk_tcga))),1]

# 风险值和通路的相关系数
ptw <- gsva(
  expr = as.matrix(exp_tpm),
  gset.idx = unstack(anno_hallmark, gene ~ term),
  method = "ssgsea",
  kcdf = "Gaussian",
  verbose = T
)
df <- cbind(meta, t(ptw[,meta$sample]))
sort(apply(df[,12:61], 2, function(x){
  cor(df$risk, x, method = 'spearman')
}))

# XJ风险值和通路的相关系数
ptw <- gsva(
  expr = as.matrix(exp_tpm_t),
  gset.idx = unstack(anno_hallmark, gene ~ term),
  method = "ssgsea",
  kcdf = "Gaussian",
  verbose = T
)
df <- cbind(filter(meta, Data == 'XJ'), t(ptw[,filter(meta, Data == 'XJ')$sample]))

corr <- sort(apply(df[,8:57], 2, function(x){
  cor(df$risk, x, method = 'spearman')
}))
corr <- as.data.frame(corr)
corr$variable <- rownames(corr)

corr <- arrange(corr, corr)
corr$variable <- factor(corr$variable, levels = corr$variable)
ggplot(corr, aes(x = variable, y = corr)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(x = 'Pathway', y = 'Corr') + 
  coord_flip() +
  theme(axis.text.y = element_text(size = 7))


# 风险值和TME的相关系数
df <- cbind(meta, t(tme_xcell[,meta$sample]))
sort(apply(df[,12:78], 2, function(x){
  cor(df$risk, x, method = 'spearman')
}))

# 风险值和受配体的相关系数
df <- cbind(meta, t(scoresLR[,meta$sample]))
corr <- sort(apply(df[,12:934], 2, function(x){
  cor(df$risk, x, method = 'spearman')
}))
LRinter.dataframe$corr <- corr[LRinter.dataframe$LR]

# 保存结果
save(hr, LRinter.dataframe, file = './LRinter.dataframe.RData')

# 筛选受配体
# 风险值正相关 + 无监督高风险组上调 + 危害性的受配体
df1 <- filter(LRinter.dataframe, corr > 0)
df2 <- filter(LRinter.dataframe, logFC > 0)
df3 <- filter(hr, HR_pvalue < 0.05 & HR > 1)
lr <- str_split(Reduce(intersect, list(df1$LR, df2$LR, rownames(df3))), '_', simplify = T)

venn_data <- list(
  Risk_pos = df1$LR,
  Cluster1_up = df2$LR,
  HR_harmful = rownamWes(df3)
)

ggVennDiagram(venn_data, label_alpha = 0, set_size = 3.5) +  
  scale_fill_gradient(low = "white", high = "blue") +  # 颜色渐变
  theme(legend.position = "right")

lr <- as.data.frame(lr)
lr$LR <- paste(lr[,1], lr[,2], sep = '_')
lr$pathway <- LRinter.dataframe$pw.name[match(lr$LR, LRinter.dataframe$LR)]
lr <- lr[complete.cases(lr),]
ddf <- data.frame(Pathway = unique(names(table(lr$pathway))),
                  Count = as.numeric(table(lr$pathway)))
ddf <- arrange(ddf, desc(Count))
ggplot() + 
  annotate(geom = "table", 
           x = 1, y = 1, 
           label = list(ddf)) +
  theme_void()

# 保存受配体结果
wb <- createWorkbook() # 创建工作簿
addWorksheet(wb, 'Logfc')
writeData(wb, 'Logfc', LRinter.dataframe, rowNames = F, colNames = T)
addWorksheet(wb, 'HR')
writeData(wb, 'HR', hr, rowNames = F, colNames = T)
addWorksheet(wb, 'Filter')
writeData(wb, 'HR', lr, rowNames = F, colNames = T)
saveWorkbook(wb, paste0("./3_lr_result.xlsx")) # 保存工作簿
