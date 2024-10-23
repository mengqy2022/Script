gene_correlation <- function(data, gene) {
  # data: 数据框，每行为一个样本的基因表达值，每列为一个基因
  # gene: 感兴趣的基因名称
  
  # 提取感兴趣基因和其他基因的表达值
  gene_data <- data[,gene]
  #gene_list <- setdiff(colnames(data),gene)
  #other_data <- data[,gene_list]
  other_data <- data[,!colnames(data) == gene, drop=FALSE]
  
  # 计算感兴趣基因和其他基因的相关性
  cor_results <- apply(other_data, 2, function(x) cor.test(gene_data, x, method = "pearson",adjust = "fdr"))
  
  # 提取相关性和P值
  cor_coef <- sapply(cor_results, function(x) x$estimate)
  p_value <- sapply(cor_results, function(x) x$p.value)
  
  # 将结果存储在数据框中
  cor_result_df <- data.frame(
    target_gene = gene,
    gene_all = colnames(other_data),
    Correlation = cor_coef,
    FDR = p_value,
    #adjusted_pvalue = adjusted_pvalue,
    stringsAsFactors = FALSE,
    row.names = NULL
  ) %>%
    na.omit() %>%
    # 过滤
    #filter(p_value < 0.05, abs(cor_coef) > 0.5)
    filter(FDR < 0.05) %>%
    return(cor_result_df)
}