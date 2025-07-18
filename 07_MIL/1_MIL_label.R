library(mclust)
library(data.table)
library(sva)
library(openxlsx)
library(survival)
library(Seurat)
library(survminer)

#### XJ临床信息 ####
follow_up <- fread('./data/seq/preprocess_result/follow_up.txt') # 随访表
meta_xj <- read.delim('./data/seq/preprocess_result/meta.txt') # 肿瘤元信息表
meta_t <- filter(meta_xj, Group == 'Tumour')

status_xj <- as.data.frame(follow_up[, c(11,9)]) # 随访信息
colnames(status_xj) <- c('time', 'status')
rownames(status_xj) <- follow_up$Tumor_Sample_Barcode

#### TCGA临床信息 ####
clinical <- fread('./data/seq/preprocess_result/TCGA-ESCA/clinical.tsv')
clinical$Overall_Survival_Months <- ifelse(clinical$vital_status == 'Dead',
                                           as.numeric(clinical$days_to_death)/30,
                                           as.numeric(clinical$days_to_last_follow_up)/30)
clinical$Live_Status <- ifelse(clinical$vital_status == 'Alive', 0, 1)
clinical <- clinical[!duplicated(clinical$case_submitter_id), ]

status_tcga <- as.data.frame(clinical[, c(159, 160)])
colnames(status_tcga) <- c('time', 'status')
rownames(status_tcga) <- clinical$case_submitter_id
status_tcga["TCGA-L5-A43J", "time"] <- 0.01

#### 分子表达谱 ####
# XJ
exp_count <- read.csv('./data/seq/preprocess_result/exp_count.csv') 
exp_tpm <- read.csv('./data/seq/preprocess_result/exp_tpm.csv') 
exp_vsd <- read.csv('./data/seq/preprocess_result/exp_vsd.csv')
rownames(exp_count) <- rownames(exp_tpm) <- rownames(exp_vsd) <- exp_count$X
exp_count <- exp_count[,-1]
exp_tpm <- exp_tpm[,-1]
exp_vsd <- exp_vsd[,-1]
exp_count_t <- exp_count[,meta_t$Tumor_Sample_Barcode]
exp_tpm_t <- exp_tpm[,meta_t$Tumor_Sample_Barcode]
exp_vsd_t <- exp_vsd[,meta_t$Tumor_Sample_Barcode]

# TCGA
agg_type = c('mean', 'sd', 'skewness', 'kurtosis')
path_output2 = './data/patch/tcga_4096/'
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

# 合并
gene <- intersect(rownames(exp_count_t), rownames(exp_tcga_count))
exp <- as.matrix(cbind(exp_count_t[gene,], exp_tcga_count[gene,]))
exp <- ComBat_seq(exp, batch = c(rep(1, ncol(exp_count_t)), rep(2, ncol(exp_tcga_count))))
exp <- round(exp) # 175

#### 靶点表达分组 ####
mclust <- Mclust(exp['CACNA1C',], modelNames = 'V')
summary(mclust)
norm = function(x){log2(x - min(x) + 1)}

df <- data.frame(
  CACNA1C = exp['CACNA1C',],
  label = as.character(mclust$classification)
)
df$CACNA1C <- norm(df$CACNA1C)
p1 <- ggplot(df, aes(x = label, y = CACNA1C, fill = label)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = c('#f9a411', '#6a75c2', '#f14166'))
p2 <- ggplot(df, aes(x = CACNA1C, y = label, fill = label)) +
  ggridges::geom_density_ridges_gradient() +
  theme_bw()+
  scale_fill_manual(values = c('#f9a411', '#6a75c2', '#f14166'))
p1/p2


#### 靶点预后价值 ####
df$id <- rownames(df)
status <- rbind(status_xj, status_tcga)
status$id <- rownames(status)
status <- left_join(status, df, by = "id")
colnames(status)[1:2] <- c('time', 'status')

fit <- survfit(Surv(time, status) ~ label, data = status)
ggsurvplot(fit = fit, 
                 data = status,
                 palette = c('#f9a411', '#6a75c2', '#f14166'),
                 pval = T,
                 risk.table = T,
                 surv.median.line = "hv", 
                 title = "XJ + TCGA",
                 ylab = "Cumulative survival (percentage)", xlab = " Time (Months)",
                 censor.shape = 124,
                 censor.size = 4,
                 conf.int = FALSE
)

tumorDiff <- survdiff(Surv(time, status) ~ label, data = status)
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

#### 靶点空间表达 ####
# PT0：1909933-6
PT0 <- Load10X_Spatial(data.dir = "./data/seq/stseq_GHRL23100290/PT0_res/outs/",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", 
                       slice = "slice1")
PT0 <- SCTransform(PT0, assay = 'Spatial')

label <- fread('/home/ljc/0_project/0_ESCC/0_STARCH/pt0/labels_PT0.csv')
position <- fread('/home/ljc/0_project/0_ESCC/0_STARCH/pt0/tissue_positions_list.csv')
position$V1 <- paste0(sprintf("%.1f", position$array_row), 'x', sprintf("%.1f", position$array_col))
position <- left_join(position, label, by = 'V1')
PT0$CNV <- position$V2[match(colnames(PT0), position$barcode)]

PT0_tumor <- subset(PT0, CNV == '0')
SpatialFeaturePlot(PT0_tumor, features = "CACNA1C", crop = TRUE, pt.size.factor = 2) + theme(legend.position = "right")
SpatialFeaturePlot(PT0, features = "CACNA1C", crop = TRUE, pt.size.factor = 2) + theme(legend.position = "right")

# PT1：1930377-13
PT1 <- Load10X_Spatial(data.dir = "./data/seq/stseq_GHRL23100290/PT1_res/outs/",
                       filename = "filtered_feature_bc_matrix.h5",
                       assay = "Spatial", 
                       slice = "slice1")
PT1 <- SCTransform(PT1, assay = 'Spatial')

label <- fread('/home/ljc/0_project/0_ESCC/0_STARCH/pt1/labels_PT1.csv')
position <- fread('/home/ljc/0_project/0_ESCC/0_STARCH/pt1/tissue_positions_list.csv')
position$V1 <- paste0(sprintf("%.1f", position$array_row), 'x', sprintf("%.1f", position$array_col))
position <- left_join(position, label, by = 'V1')
PT1$CNV <- position$V2[match(colnames(PT1), position$barcode)]

PT1_tumor <- subset(PT1, CNV == '0')
SpatialFeaturePlot(PT1_tumor, features = "CACNA1C", crop = TRUE, pt.size.factor = 2, max.cutoff = 1.3) + theme(legend.position = "right")
SpatialFeaturePlot(PT1, features = "CACNA1C", crop = TRUE, pt.size.factor = 2) + theme(legend.position = "right")


#### 切片名称 ####
samp <- read.csv('/home/ljc/0_project/0_ESCC/data/patch/xj_tcga/process_list_autogen.csv')
samp$id <- str_split(samp$slide_id, '\\.', simplify = T)[,1]
samp$id <- ifelse(str_count(samp$id) > 20, str_sub(samp$id, 1, 12), samp$id)

#### 标签构建 ####
df$id <- df$patient <- rownames(df)
df$id <- meta_t[match(df$id, meta_t$Tumor_Sample_Barcode), 'Pathological_Barcode']
df$id <- ifelse(is.na(df$id), df$patient, df$id)
df$slide_id <- samp$slide_id[match(df$id, samp$id)]
write.csv(df, file = '/home/ljc/0_project/0_ESCC/data/patch/xj_tcga/label.csv', row.names = F)

annotation <- read.csv('/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/annotations.csv')
annotation$dataset <- ifelse(substr(annotation$patient, 1, 4) == 'TCGA', 'TCGA', 'XJ')
annotation$category <- df$label[match(annotation$patient, 
                                      str_split(df$slide_id, '\\.[ndpisvs]', simplify = T)[,1])]
annotation <- annotation[complete.cases(annotation), ]
write.csv(annotation, file = '/home/ljc/0_project/0_ESCC/0_slideflow/xj_tcga_gdph/annotations.csv', row.names = F, quote = F)

