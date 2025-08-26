rm(list = ls())
setwd('/data/nas1/zhuxuying/44.YQXA-10152-9/')
if(!dir.exists('./00_rawdata')){
  dir.create('./00_rawdata')
}
setwd('./00_rawdata/')

library(tidyverse)
library(BiocManager)
#BiocManager::install('TCGAbiolinks')
library(TCGAbiolinks)
cancer_type = "TCGA-PRAD"   #肿瘤类型，这里可修改癌症类型

expquery <- GDCquery(project = cancer_type,
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification",
                     workflow.type = "STAR - Counts"
)
GDCdownload(expquery,directory = "GDCdata")
# 准备数据，转换为SummarizedExperiment对象
expquery2 <- GDCprepare(expquery,directory = "GDCdata",summarizedExperiment = T)
save(expquery2,file = "gbm.gdc_2024.rda") # 保存 rda格式

load("gbm.gdc_2024.rda")

gene_annotation <- read.table("/data/nas2/pipeline/Gencode/gencode.v36/gencode.v36.annotation.gtf.gene.probemap",header = T)
colnames(gene_annotation)[1] <- "ENSEMBL"
colnames(gene_annotation) <- c("ENSEMBL","symbol")
gene_annotation <- gene_annotation[,-c(3,4,5,6)]
#基因名称symbol ENSEMBL
#提取counts 
counts <- expquery2@assays@data@listData[["unstranded"]]
colnames(counts) <- expquery2@colData@rownames
rownames(counts) <- expquery2@rowRanges@ranges@NAMES

#基因ID转换
counts <- counts %>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation,"ENSEMBL") %>% 
  dplyr::select(-ENSEMBL) %>% #ENSEMBL一列多余删除
  dplyr::select(symbol, everything()) %>% #gene列排在前排序
  mutate(rowMean = rowMeans(.[grep('TCGA',names(.))])) %>%
  arrange(desc(rowMean)) %>%  #按照平均值从大到小排序
  distinct(symbol,.keep_all = T) %>% #留下第一个gene
  dplyr::select(-rowMean) %>% #反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1]) #把第一列gene变成行名并删除

#group
sample <- data.frame(sample=colnames(counts),num = "",group = "")
sample$num <- substr(sample$sample,14,15)
sample$num = lapply(sample$num,as.numeric)
sample$group <-  ifelse(sample$num < 10,"PCa","Control")
sample <- sample[,-2]
sample$sample <- substring(sample$sample,1,16)
sample <- sample[order(sample[,2],decreasing = T),] 
table(sample$group)  

write.table(sample,"00.group.csv",sep = "\t",row.names = T,col.names = NA,quote = F)

colnames(counts) <- substring(colnames(counts),1,16)
counts <- counts[,sample$sample]
write.table(counts,"00.counts.csv",sep = "\t",row.names = T,col.names = NA,quote = F)

#提取fpkm
fpkm <- expquery2@assays@data@listData[["fpkm_unstrand"]]
colnames(fpkm) <- expquery2@colData@rownames
rownames(fpkm) <- expquery2@rowRanges@ranges@NAMES
fpkm <- fpkm%>% 
  as.data.frame() %>% 
  rownames_to_column("ENSEMBL") %>% 
  inner_join(gene_annotation,"ENSEMBL") %>% 
  dplyr::select(-ENSEMBL) %>% #ENSEMBL一列多余删除
  dplyr::select(symbol, everything()) %>% #gene列排在前排序
  mutate(rowMean = rowMeans(.[grep('TCGA',names(.))])) %>%
  arrange(desc(rowMean)) %>%  #按照平均值从大到小排序
  distinct(symbol,.keep_all = T) %>% #留下第一个gene
  dplyr::select(-rowMean) %>% #反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1]) #把第一列gene变成行名并删除

# 把TCGA barcode切割为16位字符,并去除重复样本
colnames(fpkm) <- substring(colnames(fpkm),1,16)
fpkm <- fpkm[rownames(counts),sample$sample]
write.table(fpkm,"00.fpkms.csv",sep = "\t",row.names = T,col.names = NA,quote = F)

group<-sample

#phe
phe1 <- read_tsv("TCGA-PRAD.GDC_phenotype.tsv")
cli_dat1 <- cbind(sample = phe1$sample,
                  gender = phe1$gender.demographic,
                  age = phe1$age_at_index.demographic,
                  Gleason=phe1$primary_gleason_grade.diagnoses,
                  T_stage = phe1$ajcc_pathologic_t.diagnoses,
                  N_stage = phe1$ajcc_pathologic_n.diagnoses,
                  M_stage = phe1$ajcc_clinical_m.diagnoses
)
cli_dat <- as.data.frame(cli_dat1)
cli_dat[cli_dat == ""] <- NA
cli_dat <- na.omit(cli_dat)
rownames(cli_dat)<-cli_dat$sample
cli_dat2 <- cli_dat[group$sample[group$group == "PCa"],]
cli_dat2 <- na.omit(cli_dat2)
cli_dat2 <- cli_dat[,-1]

# #surv
phe2 <- read_tsv("TCGA-PRAD.GDC_phenotype.tsv")
phenotype <-phe2[,c('submitter_id.samples',
                    "biochemical_recurrence",
                    "days_to_first_biochemical_recurrence",
                    "days_to_second_biochemical_recurrence",
                    'days_to_third_biochemical_recurrence',
                    "days_to_last_follow_up.diagnoses")]
survival.yes <- phenotype[which(phenotype$biochemical_recurrence=='YES'),]
survival.no <- phenotype[which(phenotype$biochemical_recurrence=='NO'),]

survival.no  <- survival.no[order(survival.no$days_to_last_follow_up.diagnoses),]
survival.no <- survival.no[,c(1,2,6)]%>% na.omit()
colnames(survival.no)<-c('sample',"BCR","BCR.time")

survival.yes  <- survival.yes[order(survival.yes$days_to_first_biochemical_recurrence),]
survival.yes <- survival.yes[,c(1,2,3)]%>% na.omit()
colnames(survival.yes)<-c('sample',"BCR","BCR.time")

survival<-rbind(survival.no,survival.yes)

ss<- Reduce(intersect,list(survival$sample,group$sample[group$group == "PCa"]))
rownames(survival)<-survival$sample
survival<-survival[ss,]
survival<- na.omit(survival)#397个有生存BCR信息的疾病样本

table(survival$BCR)
# NO YES 
# 370  27 
survival$BCR<-ifelse(survival$BCR=='YES',1,0)
write.csv(survival,file = "01.survival_BCR.csv")

count1 <- counts[,survival$sample]
data1 <- fpkm[,survival$sample]
cli_dat3 <- cli_dat2[survival$sample,]

write.csv(count1,file = '01.count_PCa.csv') #397  #count为原始数据，无log
write.csv(data1,file = '01.fpkm_PCa.csv') #397  #fpkm已log化
write.csv(cli_dat3,file = "01.phenotype_PCa.csv") #397


# #验证集
# #GSE46602  预后模型验证
# gset<-getGEO("GSE46602",
#              destdir = '.',
#              GSEMatrix = T,
#              getGPL = F)
# gset$GSE46602_series_matrix.txt.gz@annotation
# expr<-as.data.frame(exprs(gset[[1]]))
# gpl<-getGEO("GPL570",destdir = '.')
# a=gset[[1]]
# #gpl<-Table(gpl)
# gpl <- idmap(gpl = 'GPL570')
# colnames(gpl)
# probe2symobl<-gpl
# # probe2symobl<-gpl %>%
# #   dplyr::select('ID','Gene Symbol')%>%
# #   filter('Gene Symbol'!='')%>%
# #   separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
# #   dplyr::select(-drop)
# probe2symobl=probe2symobl[probe2symobl$symbol!='',]
# colnames(probe2symobl)<-c('ID','symbol')
# dat.va<-expr
# dat.va$ID<-rownames(dat.va)
# dat.va$ID<-as.character(dat.va$ID)
# probe2symobl$ID<-as.character(probe2symobl$ID)
# dat.va<-dat.va %>%
#   inner_join(probe2symobl,by='ID')%>%
#   dplyr::select(-ID)%>%     ## 去除多余信息
#   dplyr::select(symbol,everything())%>%     ## 重新排列
#   mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
#   arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
#   distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
#   dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
#   tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# pd<-pData(a)
# table(pd$characteristics_ch1)
# pd2 <- subset(pd,characteristics_ch1=='tissue: prostate tumor')
# 
# survival.va <- pd2[,c("characteristics_ch1.5","characteristics_ch1.6")] #1=yes, 0=no
# survival.va$characteristics_ch1.5<-gsub('bcr: ','',survival.va$characteristics_ch1.5,fixed = T)
# survival.va$characteristics_ch1.6<-gsub('bcr_free_time: ','',survival.va$characteristics_ch1.6,fixed = T)
# survival.va2<-na.omit(survival.va)
# survival.va2$characteristics_ch1.5 <- ifelse(survival.va2$characteristics_ch1.5=='YES','1','0')
# colnames(survival.va2) <- c('bcr','bcr.time')
# survival.va2$bcr.time<-as.numeric(survival.va2$bcr.time)*30
# survival.va3<-data.frame(sample=rownames(survival.va2),bcr=survival.va2$bcr,bcr.time=survival.va2$bcr.time)
# dat.va2 <- dat.va[,survival.va3$sample]
# write.csv(dat.va,file = '02.GSE46602_expr.csv',quote = F,row.names = T)
# write.csv(survival.va3,"02.GSE46602_BCR.csv",row.names = F)   #36
