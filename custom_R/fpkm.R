#!/usr/bin/env Rscript: 这是一个用于告诉操作系统使用哪种解释器来执行脚本的声明。在这里，它指定使用R来执行脚本。

#  Extract Command Line Arguments
# args=commandArgs(TRUE): 这一行将命令行参数存储在变量args中。commandArgs(TRUE)用于获取R脚本的命令行参数。
args=commandArgs(TRUE)

# if (length(args) != 2) { ... }: 这个条件语句检查命令行参数的数量是否为2。如果不是，它会打印出用法信息并退出脚本。
if (length(args) != 2) {
  print ("usage: [input gtf] [output name]")
  q()
}

#  input <- as.character(args[1]) 和 prefix <- as.character(args[2]): 这两行将第一个和第二个命令行参数分别存储在input和prefix变量中。
# input是GTF文件的路径，而prefix是输出文件的前缀。
input <- as.character(args[1])
prefix <- as.character(args[2])
#max_y <- as.numeric(args[3])

library(tidyverse)
library(dplyr)
#自定义函数:各种转换的方法
# 这些函数对于在RNA测序数据分析中进行表达量标准化和转换非常有用，可以根据需求选择其中之一来进行数据处理。例如，countToFpkm函数用于将计数数据转换为FPKM表达量，而countToCPM函数用于将计数数据转换为CPM表达量。
# 这些标准化方法帮助比较不同基因或样本的表达量，同时考虑了数据的大小和长度等因素。

# countToTpm 函数:
# 计算每个基因或转录本的表达率（rate），其中rate = log(counts) - log(effLen)。
# 计算所有基因或转录本的表达率之和的对数（denom），并将其用于标准化。
# 将每个基因或转录本的表达率标准化为TPM，并返回结果。
# 输入参数：counts（基因或转录本的计数值），effLen（有效长度，通常是基因或转录本的长度）。
# 功能：将基因或转录本的计数值（Counts）转换为TPM（Transcripts Per Million，每百万转录本数）表达量。
countToTpm <- function(counts, effLen)
{
  #  Logarithms and Exponentials
  rate <- log(counts) - log(effLen)
  #  computes the exponential function.
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

# countToFpkm 函数:
# 计算所有基因或转录本的总计数（N），即N = sum(counts)。
# 使用公式计算FPKM，其中包括计算对数、除法和乘法运算，最终返回FPKM表达量。
# 输入参数：counts（基因或转录本的计数值），effLen（有效长度，通常是基因或转录本的长度）。
# 功能：将基因或转录本的计数值（Counts）转换为FPKM（Fragments Per Kilobase of transcript per Million mapped reads，每百万映射的读片段数）表达量。
countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

# fpkmToTpm 函数:
# 计算所有基因或转录本的FPKM的总和（sum(fpkm)）。
# 使用公式计算TPM，其中包括计算对数、减法和乘法运算，最终返回TPM表达量。
# 输入参数：fpkm（基因或转录本的FPKM表达量）。
# 功能：将基因或转录本的FPKM表达量转换为TPM（Transcripts Per Million，每百万转录本数）表达量。
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

# countToCPM 函数:
# 计算所有基因或转录本的总计数（N），即N = sum(counts)。
# 使用公式计算CPM，其中包括计算对数、除法和乘法运算，最终返回CPM表达量。
# 输入参数：counts（基因或转录本的计数值）。
# 功能：将基因或转录本的计数值（Counts）转换为CPM（Counts Per Million，每百万计数数）表达量。
# 实现步骤：
countToCPM <- function(counts)
{
  N <- sum(counts)
  exp( log(counts) + log(1e6) - log(N) )
}

#读gtf文件，计算所有外显子的长度
# 读取GTF文件：这一部分代码使用read_tsv函数从GTF文件中读取数据，筛选出外显子类型的记录，并计算每个外显子的长度。
# 然后，它从GTF属性中提取基因ID，计算每个基因的有效长度（非冗余外显子的长度），最后得到一个包含基因ID和有效长度的数据框。

gtf <- read_tsv(input, 
                comment="#", 
                col_names=c('chr','source','type','start','end',
                            'score','strand','phase','attributes')) %>%
  filter(type=='exon') %>% 
  mutate(len = end - start + 1) %>% 
  select(start, end, attributes,len)

# 计算基因的非冗余外显子的长度，即获得有效基因长度
#  读取基因计数文件：这一部分代码使用read.table函数读取基因计数文件，并将列重命名为gene_id和count。
# 然后，它将这个表格与之前计算的有效基因长度合并，以得到一个包含基因ID、计数和有效长度的数据框。
gtf$attributes %>% str_extract(., "gene_id \"[\\w|\\.]+") %>%
  str_remove(., "gene_id \"") -> gtf$gene_id
gtf %>% select(start, end, gene_id, len) %>%
  distinct(start,end,gene_id, .keep_all = T) %>%
  select(gene_id,len) %>% group_by(gene_id) %>% 
  summarise(est_len=sum(len)) -> gtf
# write.table(gtf,"gene_length.txt",row.names = F,quote = F)

#读取count的表达量矩阵，其中gene_id是和gtf文件一致的，然后和刚才计算得到的有效基因长度合并
expmat <- read.table(paste(prefix),header = F) %>%
  rename(gene_id = V1, count = V2) %>%
  inner_join(gtf , by = 'gene_id' ) %>%
  drop_na()

#计算fpkm
#  MARGIN	是一个向量，给出了函数要应用的下标。例如，对于矩阵，1 表示行，2 表示列，c(1, 2) 表示行和列。当 X 已命名为 dimnames 时，它可以是一个选择维度名称的字符向量。
#  计算FPKM：这一部分代码使用apply函数，对每一列的计数数据进行FPKM的计算，最终得到一个包含基因ID和FPKM的结果数据框。
fpkm=apply(expmat[,2:3],MARGIN = 2,function(x) countToFpkm(x,expmat$est_len ))
result <- cbind(expmat[,1], fpkm[,1])
colnames(result) <- c("gene_id", "fpkm")

#  写入结果：最后，结果数据框被写入一个文本文件，文件名为prefix加上"fpkm.txt"。
write.table(result,file = paste(prefix, "fpkm.txt", sep="."),quote = F,row.names = F)