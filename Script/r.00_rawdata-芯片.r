rm(list = ls())
#setwd("/data/nas1/zhuxuying/46.YQSZ-30111-9")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")

library(GEOquery)
library(tidyverse)
#library(lance)

GEO_data <- 'GSE140831'
gene_annotation <- 'GPL15988'

# gset <- getGEO(GEO_data , # 前面创建GEO_data对象
#                destdir = '.',
#                GSEMatrix = T,
#                getGPL = F) # getGEO函数自动从官网中获取对应的数据集


gset <- getGEO(filename = "./GSE140831_series_matrix.txt.gz", getGPL = F)

expr <- exprs(gset[[1]]) # 将获取的数据集的表达矩阵提取出来
#expr_1 <- exprs(gset_1) # 将获取的数据集的表达矩阵提取出来

##判断是否需要标准化
qx <- as.numeric(quantile(expr, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) {
  expr[which(expr <= 0)] <- 0
  expr <- log2(expr + 1)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}

expr <- as.data.frame(expr)
expr <- na.omit(expr)

##下载平台文件，提取注释信息
gpl <- getGEO(gene_annotation, destdir = '.')
gpl <- Table(gpl)
write.csv(gpl, 'gpl.csv')
colnames(gpl)
probe2symbol <- dplyr::select(gpl, 'ID', 'Gene Symbol')
probe2symbol <- filter(probe2symbol, `Gene Symbol` != '')
probe2symbol <- separate(probe2symbol, `Gene Symbol`, into = c('symbol', 'drop'), sep = '///')
probe2symbol <- dplyr::select(probe2symbol, -drop)
names(probe2symbol) <- c('ID', 'symbol')
probe2symbol <- probe2symbol[probe2symbol$symbol != ' --- ', ]
probe2symbol$symbol <- gsub(' ', '', probe2symbol$symbol)
probe2symbol <- na.omit(probe2symbol) 

dat <- expr
dat$ID <- rownames(dat)
dat$ID <- as.character(dat$ID)
probe2symbol$ID <- as.character(probe2symbol$ID)

#处理一个基因对应多行表达值的情况
dat <- dat %>%
  merge(probe2symbol, by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol, everything())%>%     ## 重新排列
  mutate(rowMean = rowMeans(.[grep('GSM', names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol, .keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

#提取信息文件，准备分组
a <- gset[[3]]
pd <- pData(a)
write.csv(pd, 'pd.csv')

colnames(pd)
table(pd$characteristics_ch1)


# pd <- pd[which(pd$characteristics_ch1 %in% c("subject status/group: Control", "subject status/group: patient with RPL")),]
group <- data.frame(sample = pd$geo_accession, group = pd$characteristics_ch1)
table(group$group)
group$group <- ifelse(group$group == 'tissue: adjacent normal tissue', 'control', 'disease')
table(group$group)
# control disease 
# 73      73

group <- group[order(group$group),]
dat <- dat[, group$sample]

write.csv(dat, file = paste0('dat.', GEO_data, '.csv'))
write.csv(group, file = paste0('group.', GEO_data, '.csv'), row.names = F)