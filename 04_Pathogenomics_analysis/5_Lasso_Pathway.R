## 病理特征预测通路表达

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
library(sva)
library(caret)
library(GSVA)
library(clusterProfiler)
library(doSNOW)
library(foreach)
library(ComplexHeatmap)
library(parallel)
library(metap)

setwd('/home/ljc/0_project/0_ESCC/')

path_output = './data/patch/xj_4096/'
path_output2 = './data/patch/tcga_4096/'
path_output3 = './data/patch/gdph_supp_4096/'

agg_type = c('mean', 'sd', 'skewness', 'kurtosis')


# 数据读入 --------------------------------------------------------------------
#### 背景文件 ####
ncbi <- read.delim('./data/seq/resource/Homo_sapiens.gene_info', check.names = F) # 基因文件
class <- read.xlsx('./data/seq/resource/pathway_class.xlsx', sheet = 1) # 通路类别

anno_hallmark <- read.gmt('./data/seq/resource/hallmark.v2023.2.Hs.symbols.gmt')
anno_kegg <- read.gmt('./data/seq/resource/kegg_legacy.v2024.1.Hs.symbols.gmt')
anno_reactome <- read.gmt('./data/seq/resource/reactome.v2024.1.Hs.symbols.gmt')
anno_cgp <- read.gmt('./data/seq/resource/c2.cgp.v2024.1.Hs.symbols.gmt')
anno_pid <- read.gmt('./data/seq/resource/pid.v2024.1.Hs.symbols.gmt')

anno <- Reduce(rbind, list(anno_hallmark, anno_kegg, anno_reactome, anno_cgp, anno_pid))
anno <- filter(anno, term %in% class$ptw)
anno <- unstack(anno, gene ~ term)


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

#### Gdph临床信息 ####
gdph <- read.xlsx('./data/wsi/gdph/GDPH_临床信息_单纯手术_20241106.xlsx')
gdph <- gdph[!duplicated(gdph$病理号),] # 去重
gdph$`生存状态（修正版）` <- ifelse(gdph$`生存状态（修正版）` == 'alive', 0, 1)
status_gdph <- as.data.frame(gdph[, c(114, 106)])
colnames(status_gdph) <- c('time', 'status')
rownames(status_gdph) <- gdph$病理号


# 去批次 ---------------------------------------------------------------------
gene <- intersect(rownames(exp_tpm_t), rownames(exp_tcga_tpm)) # 交集基因
grouplist <- as.factor(c(rep('XJ', dim(exp_tpm_t)[2]), rep('TCGA', dim(exp_tcga_tpm)[2]))) # 分组
exp <- cbind(exp_tpm_t[gene,], exp_tcga_tpm[gene,])
exp_rmbatch <- ComBat(as.matrix(exp), batch = grouplist)

exp_tpm_t <- exp_rmbatch[,1:dim(exp_tpm_t)[2]]
exp_tcga_tpm <- exp_rmbatch[,(dim(exp_tpm_t)[2]+1):dim(exp_rmbatch)[2]]


# Lasso-Pathway --------------------------------------------------------------
# xj
ssgsea_xj <- gsva(
  expr = exp_tpm_t,
  gset.idx = anno,
  method = "ssgsea",
  kcdf = "Gaussian",
  verbose = T
)

clin_xj <- meta_t
clin_xj <- cbind(clin_xj, t(ssgsea_xj[, clin_xj$Tumor_Sample_Barcode]))
rownames(clin_xj) <- clin_xj$Pathological_Barcode
clin_xj <- clin_xj[rownames(feat_xj),]
clin_xj <- clin_xj[,5:ncol(clin_xj)]
colnames(clin_xj)[1] <- 'samp'

# tcga
ssgsea_tcga <- gsva(
  expr = exp_tcga_tpm,
  gset.idx = anno,
  method = "ssgsea",
  kcdf = "Gaussian",
  verbose = T
)

clin_tcga <- clinical_f
clin_tcga <- cbind(clin_tcga, t(ssgsea_tcga[,clin_tcga$case_submitter_id]))
clin_tcga <- as.data.frame(clin_tcga)
rownames(clin_tcga) <- clin_tcga$case_submitter_id
clin_tcga <- clin_tcga[,c(2, 161:ncol(clin_tcga))]
colnames(clin_tcga)[1] <- 'samp'


#### XJ Test ####
# 并行化
detectCores()
core = 5
cl = makeCluster(core, setup_strategy = "sequential")
registerDoSNOW(cl)

# 进度条
pb <- txtProgressBar(max = length(names(anno_hallmark)), style = 3)
progress <- function(n) setTxtProgressBar(pb, n) # 进度条

# 低相关特征
cor_matrix <- cor(feat_xj_norm, use = "pairwise.complete.obs", method = 'pearson')
high_cor_vars <- findCorrelation(cor_matrix, cutoff = 0.9)
feature <- setdiff(colnames(feat_xj), colnames(feat_xj)[high_cor_vars])

df <- data.frame()
cvFoldIds <- createFolds(rownames(feat_xj), k = 5, returnTrain = T)

set.seed(42)

result <- foreach (ptw = rownames(ssgsea_xj), .packages = c("glmnet", 'survival', 'dplyr', 'stats'),
                   .options.snow = list(progress = progress)) %dopar% {
                     
                     clin <- clin_xj[, c('samp', ptw)]
                     
                     for (index in cvFoldIds) { # 五折交叉验证计算相关系数和p值
                       
                       cv <- cv.glmnet(feat_xj[index, feature],
                                       clin[,ptw][match(rownames(feat_xj)[index], clin$samp)],
                                       type.measure = "deviance",
                                       nfolds = 5,
                                       standardize = TRUE,
                                       intercept = TRUE) # 五折交叉验证确定最佳lambda
                       
                       res <- predict(cv, newx = feat_xj[-index, feature], type = "link", s = "lambda.min") # 测试集
                       res <- as.data.frame(res)
                       colnames(res)[1] <- "riskScore"
                       res$sample <- as.integer(rownames(res))
                       res <- inner_join(res, follow_up, by = c('sample' = "Pathological_Barcode"))
                       res$ptw <- clin[,ptw][match(res$sample, rownames(clin))]
                       corr <- cor.test(res$riskScore, res$ptw, method = 'spearman') # 相关系数和p值
                       
                       df <- rbind(df, data.frame(corr$estimate, 
                                                  corr$p.value))
                     }
                     
                     return(df)
                     
                   } 

# 保存结果
names(result) <- rownames(ssgsea_xj)
result_df1 <- as.data.frame(docall(rbind, lapply(result, function(x) {
  x <- x[complete.cases(x),]
  res <- c(max(x[,1]), min(x[,1]), mean(x[,1]), ifelse(nrow(x) > 1, sump(x[,2])$p, x[,2]))
  names(res) <- c('corr_max', 'corr_min','corr_mean', 'pval_comb')
  return(res)
})))
result_df1$fdr <- p.adjust(result_df1$pval_comb, method = 'fdr')
result_df1$ptw <- rownames(result_df1)
result_df1 <- left_join(result_df1, class, by = 'ptw') %>% arrange(desc(corr_mean))

# 写出结果
wb <- createWorkbook('./fig2.predicted_pathway.xlsx')
addWorksheet(wb, 'xj') # 加入工作表
writeData(wb, 'xj', result_df1, rowNames = F, colNames = T) # 写入数据
saveWorkbook(wb, paste0("./fig2.predicted_pathway.xlsx"), overwrite = T) # 保存工作簿


#### TCGA Test ####
set.seed(90)

df <- data.frame()
cvFoldIds <- createFolds(rownames(feat_tcga), k = 5, returnTrain = T)

result2 <- foreach (ptw = rownames(ssgsea_tcga), .packages = c("glmnet", 'survival', 'dplyr', 'stats'),
                    .options.snow = list(progress = progress)) %dopar% {
                      
                      clin <- clin_tcga[, c('samp', ptw)]
                      
                      for (index in cvFoldIds) { # 5折交叉验证计算相关系数和p值
                        
                        cv <- cv.glmnet(feat_tcga[index, feature],
                                        clin[,ptw][match(rownames(feat_tcga)[index], clin$samp)],
                                        type.measure = "deviance",
                                        nfolds = 5,
                                        standardize = TRUE,
                                        intercept = TRUE) # 5折交叉验证确定最佳lambda
                        
                        res <- predict(cv, newx = feat_tcga[-index, feature], type = "link", s = "lambda.min") # 测试集
                        res <- as.data.frame(res)
                        colnames(res)[1] <- "riskScore"
                        res$sample <- rownames(res)
                        res <- inner_join(res, clinical_f, by = c('sample' = "case_submitter_id"))
                        res$ptw <- clin[,ptw][match(res$sample, rownames(clin))]
                        corr <- cor.test(res$riskScore, res$ptw, method = 'spearman') # 相关系数和p值
                        
                        df <- rbind(df, data.frame(corr$estimate, 
                                                   corr$p.value))
                      }
                      
                      return(df)
                      
                    } 


# 保存结果
names(result2) <- rownames(ssgsea_tcga)
result_df2 <- as.data.frame(docall(rbind, lapply(result2, function(x) {
  x <- x[complete.cases(x),]
  res <- c(max(x[,1]), min(x[,1]), mean(x[,1]), ifelse(nrow(x) > 1, sump(x[,2])$p, x[,2]))
  names(res) <- c('corr_max', 'corr_min','corr_mean', 'pval_comb')
  return(res)
})))
result_df2$fdr <- p.adjust(result_df2$pval_comb, method = 'fdr')
result_df2$ptw <- rownames(result_df2)
result_df2 <- left_join(result_df2, class, by = 'ptw') %>% arrange(desc(corr_mean))

# 写出结果
wb <- loadWorkbook('./fig2.predicted_pathway.xlsx')
addWorksheet(wb, 'tcga') # 加入工作表
writeData(wb, 'tcga', result_df2, rowNames = F, colNames = T) # 写入数据
saveWorkbook(wb, paste0("./fig2.predicted_pathway.xlsx"), overwrite = T) # 保存工作簿

close(pb)
stopCluster(cl)


#### 可视化 ####
result_df1 <- read.xlsx('./3_predicted_pathway.xlsx', sheet = 1)
result_df2 <- read.xlsx('./3_predicted_pathway.xlsx', sheet = 2)

result_df <- merge(filter(result_df1, pval_comb < 0.25)[,c(3,6,7)], result_df2[,c(3,6)], by = 'ptw')
result_df <- result_df[complete.cases(result_df),]
result_df <- result_df[,c(3,1,2,4)]
colnames(result_df) <- c('Category', 'Pathway', 'Corr_XJ', 'Corr_TCGA')
result_df <- result_df[order(result_df$Corr_XJ, decreasing = T),]
result_df <- result_df[order(result_df$Category, decreasing = T),]
result_df$Pathway <- str_to_title(gsub('_', ' ', sub("^[^_]+_", "", result_df$Pathway)))
result_df <- result_df[!duplicated(result_df$Pathway),]

row_annotation <- data.frame(
  Category = result_df$Category,
  check.names = F
)
row_annotation$Category <- ifelse(row_annotation$Category == 'Metastasis and immune response', "Metastasis & Immune", row_annotation$Category)
rownames(row_annotation) <- result_df$Pathway
row_annotation <- rowAnnotation(df = row_annotation,
                                show_annotation_name = F,
                                annotation_legend_param = list(Category = list(direction = "vertical")),
                                col = list(Category = c("Metabolism " = "#f14166", 
                                                        "Metastasis & Immune" = "#6a75c2", 
                                                        "Proliferation " = "#f9a411")))

heat <- result_df[,c(3,4)]
rownames(heat) <- result_df$Pathway

pdf('3_predicted_pathway.pdf', width = 5, height = 10)
Heatmap(scale(as.matrix(heat)), 
        name = 'Z score',
        col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        use_raster = F,
        show_row_names = T,
        show_column_names = T,
        cluster_columns = F, 
        cluster_rows = F,
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 10),
        left_annotation = row_annotation,
        column_names_rot = 45, 
        rect_gp = gpar(col = "white", lwd = 2),
        heatmap_legend_param = list(direction = 'vertical'))
dev.off()

rm(list = ls())
gc()