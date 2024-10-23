#!/bin/bash

input="$1" # ./*_FPKM.xls  # 结果与其相同文件夹

#  在shell编程中，“EOF”通常与“<<”结合使用，“<<EOF”表示后续的输入作为子命令或子shell的输入，直到遇到“EOF”，再次返回到主调用shell，可将其理解为分界符（delimiter）。
# EOF是End of file的缩写，自定义终止符。
# 既然是分界符，那么形式自然不是固定的，这里可以将”EOF“可以进行自定义，但是前后的”EOF“必须成对出现且不能和shell命令冲突。其使用形式如下：

Rscript - <<EOF
    #  允许用户设置和检查各种全局选项，这些选项会影响 R 计算和显示结果的方式。
    options(stringsAsFactors = F)

    # input <- "./ALL_FPKM.xls"  
    input <- "$input" # 必须加引号。原因？ 字符串

    #  模式匹配和替换
    #  停止执行当前表达式并执行出错操作。
    #  sub and gsub perform replacement of the first and all matches respectively. 
    if(!grepl("_FPKM.xls", input))stop("The input must be '{Sample}_FPKM.xls'")
    output <- gsub(pattern = "_FPKM.xls",replacement = "_TPM.xls", input)

    # chcek input
    #  根据参数生成诊断信息。
    if(!file.exists(input)) stop("FileNotExist: '",input,"'")
    message("==> Start to convert '",input,"' <==")

    # read permission of input # 可读时为0
    if(file.access(input, mode = 4) != 0) stop("File '",input,"' does not have read permission")

    #  读取连接中的部分或全部文本行。
    #  根据字符向量 x 中的匹配子串拆分为子串。
    message("Trying to read ",input,"...")
    fpkm <- read.table(input, header = TRUE, sep = "\t",comment.char = "")
    header <- strsplit(readLines(input, n = 1L), "\t")[[1]]

    colnames(fpkm) <- header
    
    geneid <- colnames(fpkm)[1] # "GeneID"
    #if(!"GeneID" %in% colnames(fpkm[,2:ncol)) stop("Not find 'GeneID' column in '",input,"'")
    # if(!geneid %in% colnames(fpkm[,2:ncol)) stop("Not find '"geneid"' column in '",input,"'")
    if(sum(duplicated(fpkm[,geneid])) > 0 ) stop( "There are some genes is repeated.\n", paste0(unique(fpkm[,geneid][duplicated(fpkm[,geneid])]), collapse = ", ") )


    message("Convertting FPKM to TPM...")

    tpm <- apply(fpkm[,2:ncol(fpkm)], 2, function(t){
      exp(log(t)-log(sum(t))+log(1e6))
    })

    result_tpm <- data.frame(GeneID = fpkm[,geneid], # GeneID仅作为临时基因ID名称
                             tpm
                             )
    colnames(result_tpm)[1] <- geneid # 基因ID名称转为原文件中的定义

    write.table(result_tpm, file = output, sep = "\t", quote = F, na = "", col.names = TRUE, row.names = FALSE)
    message("Result: '",output,"'")
    message("==> Succeed to convert '",input,"' <==\n")
EOF