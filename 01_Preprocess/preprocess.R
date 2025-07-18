## 确认食管癌元信息、随访、snv数据、转录组数据
## 20240125

# install.packages('pacman')
# pacman::p_load(dplyr, stringr, tidyr, openxlsx, data.table, conflicted, BiocManager， readxl)
# BiocManager::install(c('maftools', 'GenomicFeatures', 'DESeq2', 'clusterProfiler', 'org.Hs.eg.db'))

library(dplyr)
library(stringr)
library(tidyr)
library(openxlsx)
library(readxl)
library(data.table)
library(maftools)
library(GenomicFeatures)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(cmapR)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

## 元信息表
id <- read.xlsx('./data/seq/病理和测序的-匹配.xlsx') # 125个肿瘤样本的病理-测序id匹配表
load('./data/seq/vsd_new_moregene_group_subtype_Set.rdata') # 基因表达谱
meta <- vsd_new_moregene_group_subtype_Set@phenoData@data # 元信息表

meta$id <- rownames(meta)
meta <- left_join(meta, id, by = c('id' = '测序号码')) # 将元信息表和id表结合
colnames(meta) <- c('Group', 'Subtype', 'Subtype_Group', 'Tumor_Sample_Barcode', 'Pathological_Barcode')
table(meta$Group) # 125个肿瘤样本，125个癌旁样本

write.table(meta, file = './data/seq/preprocess_result/meta.txt', row.names = F, quote = F, sep = '\t') # 输出

## 随访
follow_up <- read.csv('./data/seq/阳梦随访2023年11月.csv') # 随访总结表
colnames(follow_up) <- c('Tumor_Sample_Barcode', 'Pathological_Barcode', 'Age', 'Sex', 'Smoking_Status', 'Smoking_Age_Freq', 'Drinking_Status', 'Esophgus_Location', 
                         'Live_Status', 'Overall_Survival_Days', 'Overall_Survival_Months', 'Nation')
write.table(follow_up, file = './data/seq/preprocess_result/follow_up.txt', row.names = F, quote = F, sep = '\t') # 输出


# snv ---------------------------------------------------------------------
## maf文件格式详解：https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
snv <- fread('./data/seq/SNVs_all.txt')

# laml <- read.delim(file = '/Users/spikw/R/workspace/data/preprocess_result/output.maf') # 和yz整理的maf文件进行比较
# intersect(toupper(colnames(laml)), toupper(colnames(snv))) # 10
# setdiff(toupper(colnames(laml)), toupper(colnames(snv))) # laml_yz多了癌旁id、测序深度、AACHANGE（氨基酸变异？）
# setdiff(toupper(colnames(snv)), toupper(colnames(laml))) # snv多了基因id、转录本信息（id、所在链、外显子位置、起始位置）、cDNA变异、氨基酸密码子变异、蛋白氨基酸变异

intersect(unique(snv$Tumor_Sample_Barcode), id$测序号码) # 50
setdiff(unique(snv$Tumor_Sample_Barcode), id$测序号码) # 0
setdiff(id$测序号码, unique(snv$Tumor_Sample_Barcode)) # 75 - snv仅有50个样本，和原文相符

## 输出snv为maf文件
write.table(snv, file = './data/seq/output.maf', row.names = F, quote = F, sep = '\t') # 全部样本

write.table(filter(snv, Tumor_Sample_Barcode %in% filter(meta, Subtype == 'Subtype1')$Tumor_Sample_Barcode),
            file = './data/seq/output_subtyp1.maf', row.names = F, quote = F, sep = '\t') # subtype1
write.table(filter(snv, Tumor_Sample_Barcode %in% filter(meta, Subtype != 'Subtype1')$Tumor_Sample_Barcode),
            file = './data/seq/output_subtyp1_except.maf', row.names = F, quote = F, sep = '\t') # subtype1_except

write.table(filter(snv, Tumor_Sample_Barcode %in% filter(meta, Subtype == 'Subtype2')$Tumor_Sample_Barcode),
            file = './data/seq/output_subtyp2.maf', row.names = F, quote = F, sep = '\t') # subtype2
write.table(filter(snv, Tumor_Sample_Barcode %in% filter(meta, Subtype != 'Subtype2')$Tumor_Sample_Barcode),
            file = './data/seq/output_subtyp2_except.maf', row.names = F, quote = F, sep = '\t') # subtype2_except

write.table(filter(snv, Tumor_Sample_Barcode %in% filter(meta, Subtype == 'Subtype3')$Tumor_Sample_Barcode),
            file = './data/seq/output_subtyp3.maf', row.names = F, quote = F, sep = '\t') # subtype3
write.table(filter(snv, Tumor_Sample_Barcode %in% filter(meta, Subtype != 'Subtype3')$Tumor_Sample_Barcode),
            file = './data/seq/output_subtyp3_except.maf', row.names = F, quote = F, sep = '\t') # subtype3_except

laml <- read.maf('./data/seq/output.maf') # 重新读入为maf对象
laml@clinical.data <- as.data.table(left_join(laml@clinical.data, as.data.table(meta), by = 'Tumor_Sample_Barcode')) # 修改临床信息slot
save(laml, file = './data/seq/preprocess_result/output_maf.RData') # 保存修改后maf对象

laml_sub1 <- read.maf('./data/seq/output_subtyp1.maf')
laml_sub1@clinical.data <- as.data.table(left_join(laml_sub1@clinical.data, as.data.table(meta), by = 'Tumor_Sample_Barcode'))
save(laml_sub1, file = './data/seq/preprocess_result/output_maf_subtype1.RData')

laml_sub2 <- read.maf('./data/seq/output_subtyp2.maf')
laml_sub2@clinical.data <- as.data.table(left_join(laml_sub2@clinical.data, as.data.table(meta), by = 'Tumor_Sample_Barcode'))
save(laml_sub2, file = './data/seq/preprocess_result/output_maf_subtype2.RData')

laml_sub3 <- read.maf('./data/seq/output_subtyp3.maf')
laml_sub3@clinical.data <- as.data.table(left_join(laml_sub3@clinical.data, as.data.table(meta), by = 'Tumor_Sample_Barcode'))
save(laml_sub3, file = './data/seq/preprocess_result/output_maf_subtype3.RData')

laml_sub1_except <- read.maf('./data/seq/output_subtyp1_except.maf')
laml_sub1_except@clinical.data <- as.data.table(left_join(laml_sub1_except@clinical.data, as.data.table(meta), by = 'Tumor_Sample_Barcode'))
save(laml_sub1_except, file = './data/seq/preprocess_result/output_maf_subtype1_except.RData')

laml_sub2_except <- read.maf('./data/seq/output_subtyp2_except.maf')
laml_sub2_except@clinical.data <- as.data.table(left_join(laml_sub2_except@clinical.data, as.data.table(meta), by = 'Tumor_Sample_Barcode'))
save(laml_sub2_except, file = './data/seq/preprocess_result/output_maf_subtype2_except.RData')

laml_sub3_except <- read.maf('./data/seq/output_subtyp3_except.maf')
laml_sub3_except@clinical.data <- as.data.table(left_join(laml_sub3_except@clinical.data, as.data.table(meta), by = 'Tumor_Sample_Barcode'))
save(laml_sub3_except, file = './data/seq/preprocess_result/output_maf_subtype3_except.RData')


# cnv --------------------------------------------------------------------------
####
## 4个必要文件: scores.gistic、all_lesionsall_lesions.conf_90.txt, amp_genes.conf_90.txt, del_genes.conf_90.txt
## 1个基因层面的汇总文件: all_data_by_genes.txt
####
cnv <- read.table('./data/seq/preprocess_result/all_data_by_genes.txt')
laml.gistic <- readGistic(gisticAllLesionsFile = './data/seq/preprocess_result/all_lesions.conf_90.txt', 
                          gisticAmpGenesFile = './data/seq/preprocess_result/amp_genes.conf_90.txt', 
                          gisticDelGenesFile = './data/seq/preprocess_result/del_genes.conf_90.txt', 
                          gisticScoresFile = './data/seq/preprocess_result/scores.gistic')

getSampleSummary(laml.gistic)
getGeneSummary(laml.gistic)
getCytobandSummary(laml.gistic)

## maftools可视化
gisticChromPlot(gistic = laml.gistic, markBands = NULL) # 在线性基因组上，可视化扩增/缺失区域的gi分值
gisticBubblePlot(gistic = laml.gistic, markBands = NULL) # 可视化扩增/缺失区域的fdr 
gisticOncoPlot(gistic = laml.gistic, 
               clinicalData = getClinicalData(laml),
               clinicalFeatures = 'Subtype',
               sortByAnnotation = TRUE, 
               top = 28) # cnv瀑布图


# rnaseq -----------------------------------------------------------------------
## 读入数据
load('./data/seq/RNAseq_txi.RData')
exp_count <- txi_salmon$counts # count 
exp_tpm <- txi_salmon$abundance # tpm
exp_vsd <- as.data.frame(vsd_new_moregene_group_subtype_Set@assayData[["exprs"]]) # vsd - 基于count进行方差稳定变换后的表达值;用于无监督聚类

## 检查vsd数据
grouplist <- as.factor(meta$group)
colData <- data.frame(row.names = rownames(meta), grouplist = grouplist)
dds <- DESeqDataSetFromMatrix(countData = round(exp_count)[,rownames(meta)],
                              colData = colData,
                              design = ~ grouplist)
vsd <- vst(dds, blind=FALSE)
View(assay(vsd))
View(exp_vsd) # 该结果和vsd一致

## 取交集基因,输出三种类型的rnaseq表达谱
write.csv(exp_vsd, file = './data/seq/preprocess_result/exp_vsd.csv')
write.csv(exp_count[rownames(exp_vsd),], file = './data/seq/preprocess_result/exp_count.csv')
write.csv(exp_tpm[rownames(exp_vsd),], file = './data/seq/preprocess_result/exp_tpm.csv')

rm(list = ls())
gc()


# LINCS -------------------------------------------------------------------
calc_reverscore = function(signature, lincs_signatures) {
  
  cmap_exp_sig <- Rfast::colRanks(-1 * lincs_signatures,
                                  method = "max")
  names.list <- list(rownames(lincs_signatures), colnames(lincs_signatures))
  dimnames(cmap_exp_sig) <- names.list
  
  cmap_score_ultimate = function(sig_up, sig_down, drug_signature) {
    num_genes <- length(drug_signature)
    ks_up <- 0
    ks_down <- 0
    connectivity_score <- 0
    
    drug_signature = rank(drug_signature)
    up_tags_rank = drug_signature[as.vector(sig_up)]
    down_tags_rank = drug_signature[as.vector(sig_down)]
    up_tags_position = sort(up_tags_rank)
    down_tags_position = sort(down_tags_rank)
    num_tags_up <- length(up_tags_position)
    num_tags_down <- length(down_tags_position)
    
    if (num_tags_up > 1) {
      a_up <- 0
      b_up <- 0
      a_up <- max(sapply(seq_len(num_tags_up), function(j) {
        j/num_tags_up - up_tags_position[j]/num_genes
      }))
      b_up <- max(sapply(seq_len(num_tags_up), function(j) {
        up_tags_position[j]/num_genes - (j - 1)/num_tags_up
      }))
      if (a_up > b_up) {
        ks_up <- a_up
      }
      else {
        ks_up <- -b_up
      }
    }
    else {
      ks_up <- 0
    }
    if (num_tags_down > 1) {
      a_down <- 0
      b_down <- 0
      a_down <- max(sapply(seq_len(num_tags_down), function(j) {
        j/num_tags_down - down_tags_position[j]/num_genes
      }))
      b_down <- max(sapply(seq_len(num_tags_down), function(j) {
        down_tags_position[j]/num_genes - (j - 1)/num_tags_down
      }))
      if (a_down > b_down) {
        ks_down <- a_down
      }
      else {
        ks_down <- -b_down
      }
    }
    else {
      ks_down <- 0
    }
    if (ks_up == 0 & ks_down != 0) {
      connectivity_score <- -ks_down
    }
    else if (ks_up != 0 & ks_down == 0) {
      connectivity_score <- ks_up
    }
    else if (sum(sign(c(ks_down, ks_up))) == 0) {
      connectivity_score <- ks_up - ks_down
    }
    else {
      connectivity_score <- ks_up - ks_down
    }
    return(connectivity_score)
  }
  
  up <- signature$GeneID[signature[["log2FoldChange"]] > 0]
  down <- signature$GeneID[signature[["log2FoldChange"]] < 0]
  dz_cmap_scores <- pbapply::pbapply(cmap_exp_sig,
                                     2, FUN = function(x) cmap_score_ultimate(up,
                                                                              down, drug_signature = x))
  
  results <- data.frame(id = colnames(cmap_exp_sig),
                        cmap_score = dz_cmap_scores)
  results <- left_join(results, info_sig, by = c("id" = "sig_id"))
  results <- results[order(results$cmap_score), ]
  
  return(results)
  
}

get_cmap = function(signature, lincs_signatures){
  
  id <- list()
  
  # 计算疾病特征和敲除扰动特征的连通性分值
  sc <- calc_reverscore(signature, lincs_signatures) 
  
  # 对于每个敲除类型，选择接近于中值连通性分值的谱，每个敲除基因一个谱
  sc_median <- aggregate(cmap_score ~ pert_iname, sc, median)
  instance_median <- merge(sc_median, sc, by = "pert_iname")
  instance_median$diff <- abs(instance_median$cmap_score.x - instance_median$cmap_score.y)
  instance_min_diff <- aggregate(diff ~ pert_iname, instance_median, min) # 每个treatment距离median的最小diff
  instance_select <- merge(instance_median, instance_min_diff, by=c("pert_iname", "diff"))
  instance_select <- instance_select[!duplicated(instance_select$pert_iname), ]
  id <- instance_select$id
  
  return(id)
  
}

# ESCC疾病特征
ncbi <- fread('./data/seq/Homo_sapiens.gene_info', check.names = F) # 基因文件
load('./data/seq/preprocess_result/sig_tvsn.RData')
sig_tvsn <- left_join(sig_tvsn, ncbi[,2:3], by = c('gene' = 'Symbol'))

# 基因敲除扰动谱
# info_sig <- fread('./data/seq/resource/LINCS/GSE106127_sig_info.txt')
# info_gene <- fread('./data/seq/resource/LINCS/GSE106127_gene_info.txt')
info_sig <- fread('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_sig_info.txt') %>% filter(pert_type %in% c('trt_sh', 'trt_xpr'))
info_gene <- fread('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_gene_info.txt')

# unique(info_sig$pert_iname) # 4320个敲除基因
# table(info_sig$cell_id) # 15个细胞系
# table(info_sig$pert_type) # 3种处理类型
unique(info_sig$pert_iname) # 4371个敲除基因
table(info_sig$cell_id) # 20个细胞系
table(info_sig$pert_type) # 1种处理类型

# 排除掉扰动谱中的对照样本
cid <- filter(info_sig, pert_type != 'ctl_vector')
# gtx <- parse_gctx('./data/seq/resource/LINCS/GSE106127_level_5_modz_n119013x978.gctx', cid = cid$sig_id)
gtx <- parse_gctx('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx', cid = cid$sig_id)
gc()

# 使用ESCC疾病特征筛选敲除扰动谱
id <- get_cmap(signature = sig_tvsn, lincs_signatures = gtx@mat)
gc()

cid <- filter(cid, sig_id %in% id) 
# gtx_f <- parse_gctx('./data/seq/resource/LINCS/level_5_modz_n119013x978.gctx', cid = cid$sig_id)
gtx_f <- parse_gctx('./data/seq/resource/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx', cid = cid$sig_id)

# 整理扰动谱
lincs <- gtx_f@mat

info_sig <- filter(info_sig, sig_id %in% colnames(lincs))
info_gene <- info_gene[match(rownames(lincs), info_gene$pr_gene_id),]

colnames(lincs) <- info_sig$pert_iname # 列为敲除基因
rownames(lincs) <- info_gene$pr_gene_symbol # 行为扰动基因

lincs_rk <- Rfast::colRanks(-1 * lincs,
                            method = "max") # 由大到小排序，最大值排序1
dimnames(lincs_rk) <- dimnames(lincs)

# 保存结果
save(lincs, lincs_rk, file = './data/seq/preprocess_result/lincs_knockdown_escc_pr.RData')

rm(list = ls())
gc()

