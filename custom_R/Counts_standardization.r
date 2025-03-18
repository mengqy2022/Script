# 安装并加载必要的包
if (!requireNamespace("R6", quietly = TRUE))
    install.packages("R6")
if (!requireNamespace("argparse", quietly = TRUE))
    install.packages("argparse")

library(R6)
library(argparse)

# 创建一个R6类来包含这些函数
CountStandardizer <- R6::R6Class(
    "CountStandardizer",
    public = list(
        counts = NULL,
        effLen = NULL,
        len = NULL,
        
        initialize = function(counts = NULL, effLen = NULL, len = NULL) {
            self$counts <- counts
            self$effLen <- effLen
            self$len <- len
        },
        
        countToTpm = function() {
            if (is.null(self$counts) || is.null(self$effLen))
                stop("counts和effLen参数必须提供")
            rate <- log(self$counts) - log(self$effLen)
            denom <- log(sum(exp(rate)))
            return(exp(rate - denom + log(1e6)))
        },
        
        countToFpkm = function() {
            if (is.null(self$counts) || is.null(self$effLen))
                stop("counts和effLen参数必须提供")
            N <- sum(self$counts)
            return(exp(log(self$counts) + log(1e9) - log(self$effLen) - log(N)))
        },
        
        fpkmToTpm = function() {
            fpkm <- self$countToFpkm()
            return(exp(log(fpkm) - log(sum(fpkm)) + log(1e6)))
        },
        
        countToEffCounts = function() {
            if (is.null(self$counts) || is.null(self$effLen) || is.null(self$len))
                stop("counts, effLen和len参数必须提供")
            return(self$counts * (self$len / self$effLen))
        }
    )
)

# 创建argparse对象
parser <- ArgumentParser(description = "RNA-seq数据标准化",
                        epilog= "\n需要分别准备不同的文件：counts矩阵，有效长度矩阵，原始长度矩阵（如果执行countToEffCounts需要）。第一列是gemID，第二列是对于的计数值，要有列名。")
parser$add_argument(
    "--counts", "-c", type = "character", required = TRUE, 
    help = "counts矩阵的文件路径"
)
parser$add_argument(
    "--effLen", "-e", type = "character", required = TRUE, 
    help = "有效长度矩阵的文件路径"
)
parser$add_argument(
    "--len", "-l", type = "character", required = FALSE, 
    help = "原始长度矩阵的文件路径，如果执行countToEffCounts需要"
)
parser$add_argument(
    "--output", "-o", type = "character", required = TRUE, 
    help = "输出结果的文件路径"
)
parser$add_argument(
    "--method", "-m", type = "character", required = TRUE, choices = c("countToTpm", "countToFpkm", "fpkmToTpm", "countToEffCounts"), 
    help = "选择标准化方法"
)

# 解析参数
args <- parser$parse_args()

# 读取输入文件
counts <- read.table(args$counts, header = TRUE, row.names = 1, check.names = FALSE)
effLen <- read.table(args$effLen, header = TRUE, row.names = 1, check.names = FALSE)
len <- if (!is.null(args$len)) read.table(args$len, header = TRUE, row.names = 1, check.names = FALSE) else NULL

# 创建CountStandardizer实例
standardizer <- CountStandardizer$new(counts = counts, effLen = effLen, len = len)

# 根据方法调用相应的函数
result <- switch(
    args$method,
    "countToTpm" = standardizer$countToTpm(),
    "countToFpkm" = standardizer$countToFpkm(),
    "fpkmToTpm" = standardizer$fpkmToTpm(),
    "countToEffCounts" = standardizer$countToEffCounts(),
    stop("无效的方法, 请检查输入")
)

# 写入输出文件
write.table(result, file = args$output, quote = FALSE, sep = "\t", col.names = NA)