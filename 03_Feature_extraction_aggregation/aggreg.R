## WSI病理特征整理和聚合
## 20240402


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

path_feat = '/home/ljc/0_project/0_ESCC/data/patch/xiaoshu_ren_2048/feature/'
path_feat_2 = '/home/ljc/0_project/0_ESCC/data/patch/xiaoshu_ren_2048/feature/'
path_output = '/home/ljc/0_project/0_ESCC/data/patch/xiaoshu_ren_2048/'


# 特征整理和聚合 ---------------------------------------------------------------------
agg_type = c("mean", "sd", "skewness", "kurtosis") # 聚合方法，4个反应分布的统计量
res <- list()
res_agg <- list()

# todo_list = gsub('_2.csv', '', list.files(path_feat, pattern = '_2'))
todo_list = gsub('.csv', '', list.files(path_feat))

for (file in todo_list) {
  
  print(file)
  
  # 读入
  feat_graph <- fread(paste0(path_feat, file, '.csv'), header = T) # 图特征
  feat_cluster <- fread(paste0(path_feat_2, file, '_2.csv'), header = T) # 聚类特征
  colnames(feat_cluster) <- colnames(feat_graph)
  
  # 填充na
  feat <- rbind(feat_graph, feat_cluster)
  feat <- feat[,-1]
  feat <- apply(feat, 1, function(x){
    if (all(is.na(x))) {
      x[is.na(x)] <- 0 # 如果特征均为na，则填充为0
    } else {
      x[is.na(x)] <- mean(x, na.rm = TRUE) # 如果特征不均为na，则用该特征的平均值填充na
    }
    return(x)
  })
  
  # 整理特征
  feat_df <- as.data.frame(t(feat))
  rownames(feat_df) <- c(feat_graph$V1, sub('_', '-', feat_cluster$V1))
  res[[file]] <- as.data.frame(t(feat_df))
  
  # 聚合特征
  feat_agg <- apply(feat, 2, function(x){
    
    res <- c(mean(x, na.rm = T), # 均值
             sd(x, na.rm = T), # 标准差
             skewness(x, na.rm = T), # 偏度
             kurtosis(x, na.rm = T) # 峰度
    )
    
    return(res)
    
  })
  
  for (i in seq(agg_type)) {
    type = agg_type[i]
    res_agg[[type]][file] <- list(feat_agg[i,]) 
  }
  
}

# 写出整理特征
save(res, file = paste0(path_output, "feat.RData"))

# 写出聚合特征
wb <- createWorkbook() # 创建工作簿
for (type in agg_type) {
  addWorksheet(wb, type) # 加入工作表
  df <- as.data.frame(res_agg[[type]], check.names = F)
  rownames(df) <- c(feat_graph$V1, sub('_', '-', feat_cluster$V1))
  writeData(wb, type, df, rowNames = T, colNames = T) # 写入数据
}
saveWorkbook(wb, paste0(path_output, "feat_agg.xlsx"))









