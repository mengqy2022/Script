library(argparse)
library(tidyverse)
library(data.table)

# 创建解析器
parser <- ArgumentParser(description = "计算TPM值的脚本")

# 添加参数
parser$add_argument("--counts_file", type = "character", default = "counts.txt",
                    help = "输入文件路径，包含基因计数数据")
parser$add_argument("--output_file", type = "character", default = "tpm_results.csv",
                    help = "输出文件路径，保存TPM结果")

# 解析参数
args <- parser$parse_args()

# 读取输入文件
a1 <- fread(args$counts_file, header = TRUE, data.table = FALSE)

# counts矩阵的构建
counts <- a1[, 7:ncol(a1)]
rownames(counts) <- a1$Geneid

# 从featurecounts原始输出文件中提取Geneid、Length(转录本长度)
geneid_efflen <- subset(a1, select = c("Geneid", "Length"))
colnames(geneid_efflen) <- c("geneid", "efflen")
geneid_efflen_fc <- geneid_efflen

# 取出counts中geneid对应的efflen
efflen <- geneid_efflen[match(rownames(counts), geneid_efflen$geneid), "efflen"]

# 计算TPM
counts2TPM <- function(count = count, efflength = efflen) {
  RPK <- count / (efflength / 1000)
  PMSC_rpk <- sum(RPK) / 1e6
  RPK / PMSC_rpk
}

tpm <- as.data.frame(apply(counts, 2, counts2TPM))

# 保存结果到输出文件
write.csv(tpm, file = args$output_file)

# 打印TPM的列和
print(colSums(tpm))