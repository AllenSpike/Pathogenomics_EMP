dyn.load('/home/ljc/anaconda3/lib/libhdf5_hl.so.200') # 设置hdf5r的配置文件

library(data.table)
library(dplyr)
library(stringr)
library(openxlsx)
library(ComplexUpset)
library(ggplot2)
library(tidyr)


# 单因素cox特征 ----------------------------------------------------------------
hr_xj <- read.xlsx('./fig1.hr_result.xlsx', sheet = 1, rowNames = T) %>% filter(HR_pvalue < 0.05)
hr_tcga <- read.xlsx('./fig1.hr_result.xlsx', sheet = 2, rowNames = T) %>% filter(HR_pvalue < 0.05)
hr_gdph <- read.xlsx('./fig1.hr_result.xlsx', sheet = 3, rowNames = T) %>% filter(HR_pvalue < 0.05)

Reduce(intersect, list(rownames(hr_xj), rownames(hr_tcga), rownames(hr_gdph))) # 0
intersect(rownames(hr_xj), rownames(hr_tcga)) # 0
intersect(rownames(hr_xj), rownames(hr_gdph)) # 3

#### 
da <- data.frame(feature = union(rownames(hr_xj), union(rownames(hr_tcga), rownames(hr_gdph))),
                 GDPH = NA,
                 XJ = NA,
                 TCGA = NA)
da$GDPH <- ifelse(da$feature %in% rownames(hr_gdph), 1, 0)
da$XJ <- ifelse(da$feature %in% rownames(hr_xj), 1, 0)
da$TCGA <- ifelse(da$feature %in% rownames(hr_tcga), 1, 0)

anno <- substr(da$feature, 1, 3)
da$Interaction <- anno
anno_mapping <-  c(
  'T-T' = '#f14166', 
  'I-I' = '#6a75c2', 
  'S-S' = '#f9a411', 
  'T-I' = '#ae5b94', 
  'T-S' = '#f5733c', 
  'I-S' = '#b28d6a'
)

pdf('fig1.hr_features_upset.pdf', width = 6, height = 4)
upset(da, 
      intersect = list('GDPH', 'XJ', 'TCGA'),
      set_sizes=FALSE,
      base_annotations = list('Intersection size'= intersection_size(counts = T, 
                                                                     mapping = aes(fill = Interaction)) + 
                                scale_fill_manual(values = anno_mapping) 
      ))
dev.off()

##
dg <- merge(hr_xj[,c(1,5)], hr_gdph[,c(1,5)], by = 'row.names')
dg$interaction <- substr(dg$Row.names, 1, 3)
dg$sum <- rowSums(dg[,c(2, 4)], na.rm = T)
dg <- pivot_longer(dg, cols = starts_with('HR'), names_to = 'Group', values_to = 'HR')
dg$Group <- ifelse(dg$Group == 'HR.x', 'XJ', 'GDPH')
dg$color <- ifelse(dg$interaction == 'T-T', '#f14166',
                   ifelse(dg$interaction == 'I-I', '#6a75c2',
                          ifelse(dg$interaction == 'S-S', '#f9a411', 
                                 ifelse(dg$interaction == 'T-I', '#ae5b94',
                                        ifelse(dg$interaction == 'T-S', '#f5733c',
                                               ifelse(dg$interaction == 'I-S', '#b28d6a', NA))))))

dg <- dg[order(dg$sum),]

p1 <- ggplot(data = dg, aes(x = reorder(Row.names, HR), 
                            y = HR, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 6,
                                   color = dg$color[!duplicated(dg$Row.names)]
        )
  ) +
  labs(title = "", x = "Feature", y = "HR", fill = 'Group') +
  scale_fill_manual(values = c("XJ" = '#E41A1C', 'GDPH' = '#377EB8', 'TCGA' = '#4DAF4A'))

p2 <- ggplot(data = dg, aes(x = reorder(Row.names, HR), 
                            y = HR, color = as.factor(interaction))) + 
  geom_point() + 
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
  scale_color_manual(values = c("T-T" = "#f14166", "I-I" = "#6a75c2", "S-S" = "#f9a411", "T-I" = "#ae5b94", 'T-S' = '#f5733c','I-S' = '#b28d6a')) +
  labs(title = "", x = "", y = "", color = 'Interaction') +
  theme(legend.position = c(0.5, 0), legend.direction = 'horizontal')

pdf('fig1.hr_features_intensity.pdf', width = 7, height = 5)
p1 + p2 + plot_layout(heights = c(1, 0))
dev.off()

rm(list = ls())
gc()


# 预后分层特征 ------------------------------------------------------------------
xj <- read.xlsx('./fig1.feature_intensity_result.xlsx', sheet = 1)
tcga <- read.xlsx('./fig1.feature_intensity_result.xlsx', sheet = 2)
gdph <- read.xlsx('./fig1.feature_intensity_result.xlsx', sheet = 3)

Reduce(intersect, list(xj$ind, tcga$ind, gdph$ind)) # 0
intersect(xj$ind, tcga$ind) # 4
intersect(xj$ind, gdph$ind) # 12
intersect(tcga$ind, gdph$ind) # 5

df <- merge(gdph, merge(xj, tcga, by = 'ind', all = T), by = "ind", all = T)
colnames(df) <- c('feature', 'GDPH', 'XJ', 'TCGA')
df[,2:4][!is.na(df[,2:4])] <- 1
df[is.na(df)] <- 0

anno <- substr(df$feature, 1, 3)
df$Interaction <- anno
anno_mapping <-  c(
  'T-T' = '#f14166', 
  'I-I' = '#6a75c2', 
  'S-S' = '#f9a411', 
  'T-I' = '#ae5b94', 
  'T-S' = '#f5733c', 
  'I-S' = '#b28d6a'
)

pdf('fig1.top_features_upset.pdf', width = 6, height = 4)
upset(df, 
      intersect = list('GDPH', 'XJ', 'TCGA'),
      set_sizes=FALSE,
      base_annotations = list('Intersection size'= intersection_size(counts = T, 
                                                                     mapping = aes(fill = Interaction)) + 
                                scale_fill_manual(values = anno_mapping) 
                              ))
dev.off()

####
ds <- merge(gdph, merge(xj, tcga, by = 'ind', all = T), by = "ind", all = T)
ds$sum <- rowSums(ds[,2:4], na.rm = T)
ds <- ds[-which(apply(ds, 1, function(x) {sum(is.na(x)) == 2})),] # 去掉只在单个队列中筛到的特征
ds$interaction <- substr(ds$ind, 1, 3)

ds <- pivot_longer(ds, cols = starts_with('values'), names_to = 'Group', values_to = 'Intensity')
ds$Group <- ifelse(ds$Group == 'values', 'XJ', 
                   ifelse(ds$Group == 'values.x', 'TCGA', 'GDPH'))
ds$color <- ifelse(ds$interaction == 'T-T', '#f14166',
                        ifelse(ds$interaction == 'I-I', '#6a75c2',
                               ifelse(ds$interaction == 'S-S', '#f9a411', 
                                      ifelse(ds$interaction == 'T-I', '#ae5b94',
                                             ifelse(ds$interaction == 'T-S', '#f5733c',
                                                    ifelse(ds$interaction == 'I-S', '#b28d6a', NA))))))

ds$Intensity[is.na(ds$Intensity)] <- 0
ds <- ds[order(ds$sum),]

p1 <- ggplot(data = ds, aes(x = reorder(ind, Intensity), y = Intensity, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(legend.position = "right", 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 6,
                                   color = ds$color[!duplicated(ds$ind)]
        )
  ) +
  labs(title = "", x = "Feature", y = "Intensity", fill = 'Group') +
  scale_fill_manual(values = c("XJ" = '#E41A1C', 'GDPH' = '#377EB8', 'TCGA' = '#4DAF4A'))

p2 <- ggplot(ds, aes(x = reorder(ind, Intensity), 
                     y = Intensity, color = interaction)) + 
  geom_point() + 
  coord_flip() +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) + 
  labs(title = "", x = "", y = "", color = 'Interaction') + 
  scale_color_manual(values = c("T-T" = "#f14166", "I-I" = "#6a75c2", "S-S" = "#f9a411", "T-I" = "#ae5b94", 'T-S' = '#f5733c','I-S' = '#b28d6a')) + 
  theme(legend.position = c(0.5, 0), legend.direction = 'horizontal')

pdf('fig1.top_features_intensity.pdf', width = 7, height = 5)
p1 + p2 + plot_layout(heights = c(1, 0))
dev.off()
