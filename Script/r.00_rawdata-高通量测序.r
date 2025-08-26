rm(list = ls())
setwd('/data/nas1/zhuxuying/48.CDZK-20103-10/')
if(!dir.exists('./00_rawdata')){
  dir.create('./00_rawdata')
}
setwd('./00_rawdata/')


#训练集
GEO_data <- 'GSE77938'
AD_1 <- 'GSE140831'

# 构建下载链接
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE140831", "file=GSE140831_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");

tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames="GeneID")
ncol(tbl)
expr <- as.data.frame(tbl)
write.csv(expr, paste0(GEO_data, '__raw_count.csv'))

# load gene annotations 
apath <- paste(urld, "type=rnaseq_counts", "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")

annot <- data.table::fread(apath, header=T, quote="", stringsAsFactors=F, data.table=F)
rownames(annot) <- annot$GeneID
write.csv(annot, 'gene_annoation.csv')

table(annot$GeneType)
annot <- annot[annot$GeneType == 'protein-coding', ]
colnames(annot)
probe2symbol <- dplyr::select(annot, 'GeneID', 'Symbol')
probe2symbol <- filter(probe2symbol, `Symbol` != '')
names(probe2symbol) <- c('ID', 'symbol')
probe2symbol <- na.omit(probe2symbol)

expr$ID <- rownames(expr)
dat <- expr
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
gset <- getGEO(GEO_data ,
               destdir = '.',
               GSEMatrix = T,
               getGPL = F) 
a <- gset[[1]]
pd <- pData(a)
write.csv(pd, 'pd.csv')
colnames(pd)
table(pd$characteristics_ch1)

# pd <- pd[which(pd$characteristics_ch1 %in% c("subject status/group: Control", "subject status/group: patient with RPL")),]
group <- data.frame(sample = pd$geo_accession, group = pd$characteristics_ch1)
table(group$group)
group_1 <- group %>% filter(group %in% c('disease state: non-KTCN', 'disease state: KTCN'))
group <- group_1
group$group <- ifelse(group$group == 'disease state: non-KTCN', 'Control', 'KC')
table(group$group) 
# Control      KC 
# 25      25 

group <- group[order(group$group),]
dat <- dat[, group$sample]

write.csv(dat, file = paste0('01.dat.', GEO_data, '_count.csv'))
write.csv(group, file = paste0('01.group.', GEO_data, '.csv'), row.names = F)