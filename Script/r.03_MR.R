gc()
rm(list = ls())
options(stringsAsFactors = FALSE)

# 设置工作目录
project_dir <- "/data/nas1/mengqingyao_OD/project/Project_0899/"

setwd(project_dir)
if (!dir.exists("./03_MR")) 
  dir.create("./03_MR")

setwd("./03_MR")

# 加载所需包 ----------------------------------------------------------------
required_packages <- c(
  "TwoSampleMR", "ieugwasr", "plinkbinr", "dplyr", "data.table",
  "gwasglue", "VariantAnnotation", "MungeSumstats", "org.Hs.eg.db",
  "readr", "ggplot2", "forestplot", "patchwork", "metafor"
)

# 安装缺失的包
new_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]
if(length(new_packages)) install.packages(new_packages)

# 加载所有包
invisible(lapply(required_packages, library, character.only = TRUE))

# 1. 数据准备 ----------------------------------------------------------------
# 加载结局数据 
outcome_vcf <- "ebi-a-GCST90038656.vcf.gz"
vcf <- readVcf(outcome_vcf)
outcome <- gwasvcf_to_TwoSampleMR(vcf, type = "outcome")
outcome_all <- outcome

# 设置结局ID和名称
outcome_all$id.outcome <- "Compression_Fracture"
outcome_all$outcome <- "Compression Fracture || id:ebi-a-GCST90018828"
save(outcome_all, file = "outcome_dat.Rdata")

# 加载候选基因
OP_DEGs <- read_csv("../02_DEGs/DEG_sig.csv")
symbols <- OP_DEGs$symbol

# 基因Symbol转ENSEMBL ID
ensembl_ids <- mapIds(org.Hs.eg.db,
                      keys = symbols,
                      keytype = "SYMBOL",
                      column = "ENSEMBL")

# 2. 工具变量筛选 ------------------------------------------------------------
# 获取eQTL数据列表
ids <- read.table("/data/nas2/database/MR/eQTL_vcf/eQTL_list") %>%
  dplyr::select(V1) %>%
  pull()
cleaned_ids <- sub("eqtl-a-", "", ids)

# 初始化结果存储
results_list <- list(
  expo.res = list(),
  dat.res = list(),
  steiger.res = list(),
  pleiotropy = data.frame(),
  heterogeneity = data.frame(),
  significant_genes = data.frame()  # 修复：初始化significant_genes
)

# 存储所有协调数据
all_harmonised_data <- list()

# 主分析循环
for (i in seq_along(ensembl_ids)) {
  current_id <- ensembl_ids[i]
  current_symbol <- symbols[i]
  
  if (is.na(current_id) || !current_id %in% cleaned_ids) next
  
  cat("\nProcessing gene:", current_symbol, "(", current_id, ")\n")
  
  # 2.1 处理暴露数据
  expos_vcf <- paste0('/data/nas2/database/MR/eQTL_vcf/eqtl-a-', current_id, '.vcf.gz')
  
  if (!file.exists(expos_vcf)) {
    cat("VCF file not found for", current_id, "\n")
    next
  }
  
  vcf1 <- tryCatch(
    VariantAnnotation::readVcf(expos_vcf, "hg19"),
    error = function(e) {
      cat("Error reading VCF for", current_id, ":", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(vcf1)) next
  
  exposure_dat <- gwasvcf_to_TwoSampleMR(vcf1, type = "exposure")
  exposure_dat$id.exposure <- current_id
  exposure_dat$exposure <- paste0(current_symbol, ' || id:', current_id)
  
  # 2.2 筛选显著SNP (p < 5e-08)
  temp_dat <- exposure_dat[exposure_dat$pval.exposure < 5e-06, ]
  
  if (nrow(temp_dat) < 3) {
    cat("Less than 3 significant SNPs for", current_id, "\n")
    next
  }
  
  # 准备LD clump数据
  temp_dat$id <- temp_dat$id.exposure
  temp_dat$rsid <- temp_dat$SNP
  temp_dat$pval <- temp_dat$pval.exposure
  
  # 2.3 连锁不平衡去除
  exp_dat <- tryCatch(
    ld_clump(
      temp_dat,
      plink_bin = get_plink_exe(),
      bfile = '/data/nas2/database/MR/g1000_eur/g1000_eur',
      clump_kb = 10,
      clump_r2 = 0.001
    ),
    error = function(e) {
      cat("LD clumping failed for", current_id, ":", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(exp_dat)) next
  
  # 2.4 匹配结局数据
  outcome_dat <- subset(outcome_all, SNP %in% exp_dat$SNP)
  
  if (nrow(outcome_dat) == 0) {
    cat("No matching SNPs in outcome for", current_id, "\n")
    next
  }
  
  # 2.5 数据协调
  dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat)
  dat <- subset(dat, mr_keep == TRUE)
  
  if (nrow(dat) == 0) {
    cat("No valid harmonised data for", current_id, "\n")
    next
  }
  
  # 2.6 计算F统计量 (工具变量强度)
  dat <- dat %>%
    mutate(
      R = get_r_from_bsen(beta.exposure, se.exposure, samplesize.exposure),
      F = (samplesize.exposure - 2) * ((R^2) / (1 - R^2))
    ) %>%
    filter(F > 10)
  
  if (nrow(dat) < 3) {
    cat("Less than 3 SNPs with F > 10 for", current_id, "\n")
    next
  }
  
  # 存储协调数据
  all_harmonised_data[[current_id]] <- dat
  
  # 2.7 异质性检验
  het <- mr_heterogeneity(dat)
  
  # 2.8 多效性检验
  ple <- mr_pleiotropy_test(dat)
  
  # 2.9 Steiger方向性检验
  steiger_sl <- directionality_test(dat)
  
  # 存储检验结果
  results_list$expo.res[[current_id]] <- exp_dat
  results_list$dat.res[[current_id]] <- dat
  results_list$steiger.res[[current_id]] <- steiger_sl
  results_list$pleiotropy <- rbind(results_list$pleiotropy, ple)
  results_list$heterogeneity <- rbind(results_list$heterogeneity, het)
  
  # 检查假设条件，只标记显著基因
  if (is.na(steiger_sl$steiger_pval) | is.na(ple$pval) | is.na(het$Q_pval[2])) next
  if (ple$pval < 0.05 | !steiger_sl$correct_causal_direction | steiger_sl$steiger_pval > 0.05) next
  
  # 存储显著基因信息
  sig_gene <- data.frame(
    gene_symbol = current_symbol,
    ensembl_id = current_id,
    nsnp = nrow(dat),
    het_qvalue = het$Q_pval[2],
    ple_pvalue = ple$pval,
    steiger_direction = steiger_sl$correct_causal_direction,
    steiger_pvalue = steiger_sl$steiger_pval,
    stringsAsFactors = FALSE
  )
  results_list$significant_genes <- rbind(results_list$significant_genes, sig_gene)
}

# 保存初步筛选结果
saveRDS(results_list, file = "IVs_analysis_results.rds")
saveRDS(all_harmonised_data, file = "all_harmonised_data.rds")

results_list <- readRDS(file = "IVs_analysis_results.rds")
all_harmonised_data <- readRDS(file = "all_harmonised_data.rds")

# 3.对显著基因进行多方法MR分析 -----------------------------------------------------------
# 读取显著基因
sig_genes <- results_list$significant_genes
sig_ids <- unique(sig_genes$ensembl_id)

if (length(sig_ids) > 0) {
  # 多方法MR分析
  mr_methods <- c(
    "mr_egger_regression",
    "mr_ivw",
    "mr_weighted_median",
    "mr_simple_mode",
    "mr_weighted_mode"
  )
  
  # 执行所有显著基因的多方法分析
  all_mr_results <- lapply(sig_ids, function(id) {
    current_dat <- all_harmonised_data[[id]]
    gene_symbol <- sig_genes[sig_genes$ensembl_id == id, "gene_symbol"][1]
    
    # 根据异质性选择主要方法（用于significant_results）
    het <- mr_heterogeneity(current_dat)
    ivw_het_pval <- het[het$method == "Inverse variance weighted", "Q_pval"]
    
    if (length(ivw_het_pval) > 0 && ivw_het_pval < 0.05) {
      primary_result <- mr(current_dat, method_list = "mr_ivw_mre") %>% generate_odds_ratios()
    } else {
      primary_result <- mr(current_dat, method_list = "mr_ivw") %>% generate_odds_ratios()
    }
    
    # 所有方法的完整分析（用于final_results_table）
    all_methods_result <- mr(current_dat, method_list = mr_methods) %>% generate_odds_ratios()
    
    # 异质性检验
    het_results <- mr_heterogeneity(current_dat) %>%
      filter(method == "Inverse variance weighted") %>%
      dplyr::select(Q, Q_df, Q_pval)
    
    # 合并结果
    list(
      primary = cbind(primary_result, het_results) %>% 
        mutate(gene_symbol = gene_symbol, ensembl_id = id),
      all_methods = cbind(all_methods_result, het_results) %>% 
        mutate(gene_symbol = gene_symbol, ensembl_id = id)
    )
  })
  
  # 分离主要结果和所有方法结果
  primary_results <- lapply(all_mr_results, function(x) x$primary) %>% bind_rows()
  all_methods_results <- lapply(all_mr_results, function(x) x$all_methods) %>% bind_rows()
  
  # 保存显著结果 - 只保留Inverse variance weighted方法
  significant_results <- primary_results %>%
    mutate(
      OR_CI = sprintf("%.2f (%.2f-%.2f)", or, or_lci95, or_uci95),
      P_value = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval))
    ) %>%
    dplyr::select(
      gene_symbol, ensembl_id, method, nsnp, b, se, pval,
      or, or_lci95, or_uci95, OR_CI, P_value, Q, Q_df, Q_pval
    )
  
  write.csv(significant_results, "IVs_significant_MR_results.csv", row.names = FALSE)
  
  # 保存所有方法结果
  final_results_table <- all_methods_results %>%
    mutate(
      OR_CI = sprintf("%.2f (%.2f-%.2f)", or, or_lci95, or_uci95),
      P_value = ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval)),
      Heterogeneity = sprintf("Q=%.2f, p=%.3f", Q, Q_pval)
    ) %>%
    dplyr::select(
      gene_symbol, ensembl_id, method, nsnp, b, se, pval,
      or, or_lci95, or_uci95, OR_CI, P_value, Heterogeneity
    ) %>%
    arrange(gene_symbol, method)
  
  write.csv(final_results_table, "Algorithm_5_results_table.csv", row.names = FALSE)
  
  # 回答你的问题：两个文件中"Inverse variance weighted"的OR值应该相同
  ivw_sig <- significant_results %>% filter(method == "Inverse variance weighted")
  ivw_alg <- final_results_table %>% filter(method == "Inverse variance weighted")
  
  if (nrow(ivw_sig) == nrow(ivw_alg) && all(ivw_sig$or == ivw_alg$or)) {
    cat("✓ IVs_significant_MR_results.csv和Algorithm_5_results_table.csv中Inverse variance weighted的OR值相同\n")
  } else {
    cat("⚠ 警告：两个文件中Inverse variance weighted的OR值不一致\n")
  }
  
} else {
  cat("No significant genes found.\n")
  # 创建空的输出文件
  write.csv(data.frame(), "IVs_significant_MR_results.csv", row.names = FALSE)
  write.csv(data.frame(), "Algorithm_5_results_table.csv", row.names = FALSE)
}

# 4 多方法结果筛选 --------------------------------------------------------
# 统计每个基因的显著方法数量
method_counts <- final_results_table %>%
  group_by(ensembl_id, gene_symbol) %>%
  summarise(
    total_methods = n(),
    sig_methods = sum(pval < 0.05),
    ivw_pval = pval[which(method == "Inverse variance weighted")][1],  # 确保只取第一个匹配项
    .groups = "drop"
  ) %>%
  filter(!is.na(ivw_pval))

# 应用筛选标准
filtered_results <- method_counts %>%
  mutate(
    result_status = case_when(
      sig_methods >= 3 ~ "Robust (≥3 methods significant)",
      sig_methods >= 1 & ivw_pval < 0.05 ~ "Acceptable (IVW significant)",
      TRUE ~ "Excluded"
    )
  ) %>%
  filter(result_status != "Excluded")

write.csv(filtered_results, "Method_5_screening_report.csv", row.names = FALSE)

filtered_results <- read_csv("Method_5_screening_report.csv")
final_results_table <- read_csv(file = "Algorithm_5_results_table.csv")

final_results_table <- final_results_table %>% 
  filter(gene_symbol %in% filtered_results$gene_symbol)

# 森林图 - 仅保留IVW方法
forest_data <- final_results_table %>%
  # 仅保留IVW方法的结果
  filter(method == "Inverse variance weighted") %>%
  # 按基因符号排序
  arrange(gene_symbol) %>%
  # 添加分组标识
  mutate(group = gene_symbol)

# 重新生成标签矩阵（确保与过滤后的数据一致）
label_matrix <- cbind(
  forest_data$gene_symbol,
  as.character(forest_data$method),
  forest_data$nsnp,
  forest_data$OR_CI,
  forest_data$P_value,
  forest_data$Heterogeneity
)

# 添加表头
label_matrix <- rbind(
  c("Gene", "Method", "nSNP", "OR (95% CI)", "P-value", "Heterogeneity"),
  label_matrix
)

pdf("Method_forest_plot_grouped.pdf", width = 16, height = nrow(forest_data)*0.3 + 2)
forestplot(
  labeltext = label_matrix,
  mean = c(NA, forest_data$or),
  lower = c(NA, forest_data$or_lci95),
  upper = c(NA, forest_data$or_uci95),
  is.summary = c(TRUE, rep(FALSE, nrow(forest_data))),
  xlog = TRUE,
  col = fpColors(box = "#1c61b6", line = "#1c61b6", summary = "red3"),
  graphwidth = unit(80, "mm"),
  xticks = c(0.5, 1, 2, 5),
  hrzl_lines = list("2" = gpar(lty = 2, col = "#999999")),
  fn.ci_norm = fpDrawNormalCI,
  boxsize = 0.2,
  line.margin = unit(0.5, "cm")
)
dev.off()

# 设置PNG输出参数（调整高度和分辨率）
png("Method_forest_plot_grouped.png", 
    width = 16,              # 宽度16英寸
    height = nrow(forest_data)*0.3 + 2,  # 动态高度
    units = "in",            # 尺寸单位为英寸
    res = 300)               # 分辨率300dpi（高清输出）

forestplot(
  labeltext = label_matrix,
  mean = c(NA, forest_data$or),
  lower = c(NA, forest_data$or_lci95),
  upper = c(NA, forest_data$or_uci95),
  is.summary = c(TRUE, rep(FALSE, nrow(forest_data))),
  xlog = TRUE,
  col = fpColors(
    box = "#1c61b6", 
    line = "#1c61b6", 
    summary = "red3"
  ),
  graphwidth = unit(80, "mm"),
  xticks = c(0.5, 1, 2, 5),
  hrzl_lines = list(
    "2" = gpar(lty = 2, col = "#999999")  # 灰色虚线分隔
  ),
  fn.ci_norm = fpDrawNormalCI,  # 统一图形样式
  boxsize = 0.2,
  line.margin = unit(0.5, "cm"),
  # 添加标题（可选）
  title = "Mendelian Randomization Results by Gene (IVW Method Only)"
)
dev.off()

sig_results <- read.csv("Method_5_screening_report.csv")
dat_all <- bind_rows(results_list$dat.res)
sig_ids <- unique(sig_results$ensembl_id)
dat_sig <- subset(dat_all, id.exposure %in% sig_ids)

# 创建基因ID到Symbol的映射
gene_mapping <- data.frame(
  ensembl_id = ensembl_ids,
  gene_symbol = symbols,
  stringsAsFactors = FALSE
)

# 5.2 散点图
scatter_plots <- lapply(sig_ids, function(id) {
  current_dat <- dat_sig[dat_sig[,"id.exposure"] == id, ]
  gene_symbol <- gene_mapping %>%
    filter(ensembl_id == id, !is.na(gene_symbol)) %>%
    pull(gene_symbol)
  
  p <- mr_scatter_plot(mr(current_dat, method_list = mr_methods), current_dat)[[1]] +
    ggtitle(gene_symbol) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
  return(p)
})

# 合并散点图
n_plots <- length(scatter_plots)
ncol <- ceiling(sqrt(n_plots))
nrow <- ceiling(n_plots / ncol)

combined_scatter <- wrap_plots(scatter_plots, ncol = ncol, nrow = nrow) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
  )

ggsave("Correlation_scatter_plots.pdf", combined_scatter,
       width = 6.5 * ncol, height = 5 * nrow, limitsize = FALSE)
ggsave("Correlation_scatter_plots.png", combined_scatter,
       width = 6.5 * ncol, height = 5 * nrow, dpi = 300, limitsize = FALSE)

# 诊断效能分析 - SNP层面森林图
# 读取显著结果
sig_results <- read.csv("Method_5_screening_report.csv")
dat_all <- bind_rows(results_list$dat.res)
sig_ids <- unique(sig_results$ensembl_id)
dat_sig <- subset(dat_all, id.exposure %in% sig_ids)

# 为每个显著基因创建单独的SNP层面森林图
top_genes <- sig_results %>%
  arrange(ivw_pval) %>%
  pull(ensembl_id)

# 创建存储单独森林图的列表
individual_plots <- list()

# 创建输出文件夹
output_dir <- "Diagnostic_SNPs_forest_single"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 为每个基因创建单独的森林图（根据SNP数量调整图片大小）
for (i in seq_along(top_genes)) {
  id <- top_genes[i]
  current_dat <- dat_sig[dat_sig$id.exposure == id, ]
  gene_symbol <- gene_mapping %>%
    filter(ensembl_id == id, !is.na(gene_symbol)) %>%
    pull(gene_symbol)
  
  if (nrow(current_dat) < 1) next
  
  # 计算SNP数量
  snp_count <- nrow(current_dat)
  
  # 根据SNP数量动态调整图片高度
  # 基础高度 + 每个SNP增加的高度
  base_height <- 120  # 基础高度
  per_snp_height <- 40  # 每个SNP增加的高度
  plot_height <- base_height + (snp_count * per_snp_height)
  
  # 确保高度在合理范围内
  plot_height <- max(400, min(plot_height, 8000))  # 最小400px，最大2000px
  
  # 计算IVW结果
  ivw_result <- mr(current_dat, method_list = "mr_ivw") %>% generate_odds_ratios()
  
  # 准备SNP层面数据
  snp_data <- current_dat %>%
    mutate(
      or = exp(beta.outcome),
      or_lci95 = exp(beta.outcome - 1.96 * se.outcome),
      or_uci95 = exp(beta.outcome + 1.96 * se.outcome),
      p_value = 2 * pnorm(-abs(beta.outcome / se.outcome))
    )
  
  # 创建合并数据（SNP + IVW汇总）
  forest_data <- list()
  
  # 添加SNP数据
  for (j in 1:nrow(snp_data)) {
    forest_data[[paste0("snp_", j)]] <- data.frame(
      snp_name = snp_data$SNP[j],
      or = snp_data$or[j],
      or_lci95 = snp_data$or_lci95[j],
      or_uci95 = snp_data$or_uci95[j],
      p_value = snp_data$p_value[j],
      is_summary = FALSE
    )
  }
  
  # 添加IVW汇总数据
  forest_data[["ivw"]] <- data.frame(
    snp_name = "ALL",
    or = ivw_result$or,
    or_lci95 = ivw_result$or_lci95,
    or_uci95 = ivw_result$or_uci95,
    p_value = ivw_result$pval,
    is_summary = TRUE
  )
  
  forest_df <- bind_rows(forest_data)
  
  # 创建标签文本
  labeltext <- cbind(
    c("SNP", forest_df$snp_name),
    c("OR (95% CI)", 
      sprintf("%.2f (%.2f-%.2f)", forest_df$or, forest_df$or_lci95, forest_df$or_uci95)),
    c("P-value", 
      ifelse(forest_df$p_value < 0.001, "< 0.001", sprintf("%.3f", forest_df$p_value)))
  )
  
  # 根据SNP数量调整字体大小
  text_size <- ifelse(snp_count > 10, 0.7, 0.8)
  
  # 保存单个森林图到指定文件夹
  png_filename <- file.path(output_dir, paste0("Diagnostic_SNP_forest_", gene_symbol, "_", id, ".png"))
  png(png_filename, width = 1600, height = plot_height, res = 150)
  
  forest_plot <- forestplot(
    labeltext = labeltext,
    mean = c(NA, forest_df$or),
    lower = c(NA, forest_df$or_lci95),
    upper = c(NA, forest_df$or_uci95),
    is.summary = c(TRUE, forest_df$is_summary),
    xlog = TRUE,
    col = fpColors(box = "#1c61b6", line = "#1c61b6", summary = "red3"),
    graphwidth = unit(80, "mm"),
    xticks = c(0.2, 0.5, 1, 2, 5),
    hrzl_lines = list("2" = gpar(lty = 2, col = "#999999")),
    txt_gp = fpTxtGp(
      label = gpar(cex = text_size),
      ticks = gpar(cex = text_size * 0.9),
      xlab = gpar(cex = text_size * 1.1)),
    title = paste0("Diagnostic Performance: ", gene_symbol, " (", id, ")"),
    boxsize = 0.2,
    clip = c(0.1, 10)
  )
  
  print(forest_plot)
  dev.off()
  
  # 存储图形对象用于后续合并
  individual_plots[[i]] <- forest_plot
  
  # 同时保存PDF版本到文件夹（同样调整高度）
  pdf_height <- plot_height / 100 * 0.35  # 转换为英寸
  pdf_filename <- file.path(output_dir, paste0("Diagnostic_SNP_forest_", gene_symbol, "_", id, ".pdf"))
  pdf(pdf_filename, width = 8, height = pdf_height)
  print(forest_plot)
  dev.off()
  
  cat("已保存", gene_symbol, "的森林图，包含", snp_count, "个SNP，图片高度:", plot_height, "px\n")
}

# 合并所有小图为一张大图
library(patchwork)
library(ggplot2)

# 首先检查 individual_plots 的结构
if (is.list(individual_plots)) {
  
  # 转换每个 gforge_forestplot 对象为 ggplot 或 grob
  converted_plots <- lapply(individual_plots, function(plot_obj) {
    
    if (inherits(plot_obj, "gforge_forestplot")) {
      # 方法1: 尝试提取 ggplot 组件
      if (!is.null(plot_obj$plot) && inherits(plot_obj$plot, "ggplot")) {
        return(plot_obj$plot)
      }
      # 方法2: 尝试访问第一个元素
      else if (length(plot_obj) > 0 && inherits(plot_obj[[1]], "ggplot")) {
        return(plot_obj[[1]])
      }
      # 方法3: 转换为 grob 对象
      else {
        return(grid::grid.grabExpr(print(plot_obj)))
      }
    } else if (inherits(plot_obj, "ggplot")) {
      return(plot_obj)  # 已经是 ggplot，直接返回
    } else {
      # 对于其他类型，尝试转换为 grob
      tryCatch({
        grid::grid.grabExpr(print(plot_obj))
      }, error = function(e) {
        warning("无法转换对象: ", class(plot_obj))
        return(NULL)
      })
    }
  })
  
  # 移除转换失败的对象（NULL）
  converted_plots <- converted_plots[!sapply(converted_plots, is.null)]
  
  # 现在使用 wrap_plots
  combined_plot <- wrap_plots(converted_plots, ncol = ) + 
    plot_annotation(
      title = "Diagnostic Performance of Significant Genes on Outcome (IVW Method)",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    )
  
  print(combined_plot)
  
} else if (inherits(individual_plots, "gforge_forestplot")) {
  # 如果 individual_plots 是单个 gforge_forestplot 对象
  # 转换为 ggplot 或直接显示
  if (!is.null(individual_plots$plot)) {
    print(individual_plots$plot)
  } else {
    print(individual_plots)
  }
}

# 保存合并的大图到当前目录
png("Diagnostic_Combined_SNPs_forest.png", width = 8000, height = 6000, res = 150)
print(combined_plot)
dev.off()

pdf("Diagnostic_Combined_SNPs_forest.pdf", width = 60, height = 50)
print(combined_plot)
dev.off()

# 5.3 漏斗图 (随机性判断)
funnel_plots <- lapply(sig_ids, function(id) {
  current_dat <- dat_sig[dat_sig[,"id.exposure"] == id, ]
  gene_symbol <- gene_mapping %>%
    filter(ensembl_id == id, !is.na(gene_symbol)) %>%
    pull(gene_symbol)
  
  # 使用单个SNP的分析结果
  res <- mr_singlesnp(current_dat)
  
  # 跳过无效结果
  if (nrow(res) == 0 || all(is.na(res$b))) return(NULL)
  
  # 准备漏斗图数据
  funnel_data <- data.frame(
    beta = res$b,
    se = res$se,
    snp = res$SNP
  )
  
  # 安全获取IVW估计值（处理可能的多行结果）
  ivw_res <- tryCatch({
    mr(current_dat, method_list = "mr_ivw")
  }, error = function(e) NULL)
  
  ivw_beta <- if (!is.null(ivw_res) && nrow(ivw_res) > 0) {
    ivw_res$b[which(ivw_res$method == "Inverse variance weighted")]  # 精确匹配IVW行
  } else NA_real_
  
  # 绘制漏斗图
  p <- ggplot(funnel_data, aes(x = beta, y = 1/se)) +
    geom_point(shape = 1, size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    {
      if (length(ivw_beta) == 1 && !is.na(ivw_beta)) 
        geom_vline(xintercept = ivw_beta, color = "red")
    } +
    labs(title = paste0(gene_symbol, "\n", id),
         x = "SNP Effect Size", 
         y = "Precision (1/SE)") +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 8))
  
  return(p)
}) %>% 
  purrr::compact()

# 合并图形
n_plots <- length(funnel_plots)
ncol <- ceiling(sqrt(n_plots))
nrow <- ceiling(n_plots / ncol)

combined_funnel <- wrap_plots(funnel_plots, ncol = ncol, nrow = nrow) +
  plot_annotation(title = "MR Funnel Plots (Individual SNP Effects)",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

# 保存图形
ggsave("Randomness_combined_funnel_plots.pdf", combined_funnel,
       width = min(6 * ncol, 49), 
       height = min(4 * nrow, 49), 
       limitsize = FALSE)

ggsave("Randomness_combined_funnel_plots.png", combined_funnel,
       width = min(6 * ncol, 49), 
       height = min(4 * nrow, 49), 
       limitsize = FALSE)

# 5. 敏感性分析 ------------------------------------------------------------
# 异质性检验结果汇总表
het_summary <- results_list$heterogeneity %>%
  filter(id.exposure %in% sig_ids) %>% 
  filter(method == "Inverse variance weighted") %>%
  dplyr::mutate(
    gene_symbol = gene_mapping$gene_symbol[match(id.exposure, gene_mapping$ensembl_id)],
    Q_pvalue_label = ifelse(Q_pval < 0.001, "< 0.001", sprintf("%.3f", Q_pval)),
    Model_choice = ifelse(Q_pval < 0.05, "Random effects", "Fixed effects")
  ) %>%
  dplyr::select(
    Gene = gene_symbol,
    ENSEMBL_ID = id.exposure,
    Q_value = Q,
    Q_df,
    Q_pvalue = Q_pval,
    Q_pvalue_label,
    Model_choice
  ) %>%
  arrange(Q_pvalue)

# 保存异质性检验结果
write.csv(het_summary, "Heterogeneity_test_results.csv", row.names = FALSE)

# 异质性检验漏斗图
het_funnel_plots <- lapply(sig_ids, function(id) {
  current_dat <- dat_sig[dat_sig$id.exposure == id, ]
  gene_symbol <- gene_mapping %>%
    filter(ensembl_id == id, !is.na(gene_symbol)) %>%
    pull(gene_symbol)
  
  # 计算IVW估计值
  ivw_res <- mr(current_dat, method_list = "mr_ivw")
  ivw_beta <- ivw_res$b[1]
  
  # 准备漏斗图数据
  funnel_data <- data.frame(
    beta = current_dat$beta.exposure / current_dat$beta.outcome,
    se = abs(current_dat$se.exposure / current_dat$beta.outcome),
    snp = current_dat$SNP
  )
  
  # 绘制漏斗图
  p <- ggplot(funnel_data, aes(x = beta, y = 1/se)) +
    geom_point(shape = 1, size = 2, color = "#1c61b6") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = ivw_beta, color = "red", linetype = "solid") +
    labs(title = paste0(gene_symbol, " (", id, ")"),
         x = "MR Estimate", 
         y = "Precision (1/SE)",
         caption = paste0("IVW estimate: ", round(ivw_beta, 3))) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 8, hjust = 0.5),
          plot.caption = element_text(size = 6))
  
  return(p)
})

# 合并异质性漏斗图
n_plots <- length(het_funnel_plots)
ncol <- ceiling(sqrt(n_plots))
nrow <- ceiling(n_plots / ncol)

combined_het_funnel <- wrap_plots(het_funnel_plots, ncol = ncol, nrow = nrow) +
  plot_annotation(title = "Heterogeneity Assessment Funnel Plots",
                  subtitle = "Red line represents IVW estimate",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

# 保存图形
ggsave("Heterogeneity_funnel_plots.pdf", combined_het_funnel,
       width = min(6 * ncol, 49), 
       height = min(4 * nrow, 49), 
       limitsize = FALSE)

ggsave("Heterogeneity_funnel_plots.png", combined_het_funnel,
       width = min(6 * ncol, 49), 
       height = min(4 * nrow, 49), 
       dpi = 300, limitsize = FALSE)

# 水平多效性检验 (MR-Egger和MR-PRESSO)
if (all(c("TwoSampleMR", "MRPRESSO") %in% installed.packages())) {
  library(TwoSampleMR)
  library(MRPRESSO)
  
  # 创建结果数据框
  pleiotropy_results <- data.frame(
    gene_symbol = character(),
    ensembl_id = character(),
    egger_intercept = numeric(),
    egger_pval = numeric(),
    presso_global_pval = numeric(),
    presso_outliers = character(),
    stringsAsFactors = FALSE
  )
  
  for (id in sig_ids) {
    current_dat <- dat_sig[dat_sig$id.exposure == id, ]
    
    gene_symbol <- gene_mapping %>%
      dplyr::filter(ensembl_id == id, !is.na(gene_symbol)) %>%
      dplyr::pull(gene_symbol) %>%
      dplyr::first()
    
    # 1. MR-Egger检验
    egger_test <- mr_pleiotropy_test(current_dat)
    
    # 2. MR-PRESSO检验
    presso <- tryCatch({
      mr_presso(BetaOutcome = "beta.outcome",
                BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome",
                SdExposure = "se.exposure",
                data = current_dat,
                OUTLIERtest = TRUE,
                DISTORTIONtest = TRUE,
                NbDistribution = 1000,
                SignifThreshold = 0.05)
    }, error = function(e) {
      return(NULL)
    })
    
    # 处理PRESSO结果
    presso_global_pval <- NA
    presso_outliers <- "None"
    
    if (!is.null(presso)) {
      presso_global_pval <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
      if (!is.null(presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)) {
        presso_outliers <- paste(presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`, collapse = ", ")
      }
    }
    
    # 添加到结果表
    pleiotropy_results <- rbind(pleiotropy_results, data.frame(
      gene_symbol = gene_symbol,
      ensembl_id = id,
      egger_intercept = egger_test$egger_intercept,
      egger_pval = egger_test$pval,
      presso_global_pval = presso_global_pval,
      presso_outliers = presso_outliers,
      stringsAsFactors = FALSE
    ))
  }
  
  # 添加结果解释列
  pleiotropy_results$egger_conclusion <- ifelse(pleiotropy_results$egger_pval > 0.05, 
                                                "No significant pleiotropy", 
                                                "Significant pleiotropy detected")
  
  pleiotropy_results$presso_conclusion <- ifelse(pleiotropy_results$presso_global_pval > 0.05 | is.na(pleiotropy_results$presso_global_pval),
                                                 "No significant pleiotropy",
                                                 "Significant pleiotropy detected")
  
  # 保存结果为CSV文件
  write.csv(pleiotropy_results, file = "pleiotropy_test_results.csv", row.names = FALSE)
}

# 留一法检验
loo_plots <- lapply(sig_ids, function(id) {
  # 获取当前基因的数据
  current_dat <- dat_sig[dat_sig$id.exposure == id, ]
  
  # 安全获取基因符号
  gene_symbol <- gene_mapping %>%
    dplyr::filter(ensembl_id == id, !is.na(gene_symbol)) %>%
    dplyr::pull(gene_symbol) %>%
    dplyr::first()
  
  # 执行留一法分析
  loo_res <- TwoSampleMR::mr_leaveoneout(current_dat)
  
  # 添加颜色列，标记是否为"ALL"结果
  loo_res$color <- ifelse(loo_res$SNP == "All", "All", "Individual")
  
  # 手动创建ggplot对象 - 更可靠的方法
  p <- ggplot2::ggplot(loo_res, ggplot2::aes(y = SNP, x = b, color = color)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = b - 1.96*se, 
                                         xmax = b + 1.96*se, 
                                         color = color), 
                            height = 0) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c("All" = "red", "Individual" = "black"), 
                                guide = "none") +  # 不显示图例
    ggplot2::labs(
      title = gene_symbol,
      x = "Leave-one-out MR estimate",
      y = "SNP removed"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 6)
    )
  
  return(p)
})

# 合并留一法检验图
combined_loo <- wrap_plots(loo_plots, ncol = ncol, nrow = nrow)
ggsave("Leave_one_combined_out_plots.pdf", combined_loo,
       width = 6 * ncol, height = 5 * nrow, limitsize = FALSE)
ggsave("Leave_one_combined_out_plots.png", combined_loo,
       width = 6 * ncol, height = 5 * nrow, limitsize = FALSE)

# 6. Steiger方向性检验 ------------------------------------------------------
# 提取并处理Steiger检验结果
steiger_results <- bind_rows(results_list$steiger.res) %>%
  # 添加基因符号信息
  mutate(
    gene_symbol = gene_mapping$gene_symbol[match(id.exposure, gene_mapping$ensembl_id)],
    # 将逻辑值转换为更易读的标签
    direction_correct = ifelse(correct_causal_direction, "TRUE", "FALSE"),
    steiger_significant = ifelse(steiger_pval < 0.05, "Significant", "Non-significant"),
    # 计算R2比例
    r2_ratio = snp_r2.exposure/snp_r2.outcome,
    # 添加假设检验结果
    assumption_met = ifelse(correct_causal_direction & steiger_pval < 0.05, 
                            "Valid", "Invalid")
  ) %>%
  # 选择并重命名关键列
  dplyr::select(
    Gene = gene_symbol,
    ENSEMBL_ID = id.exposure,
    Exposure = exposure,
    Outcome = outcome,
    SNP_R2_Exposure = snp_r2.exposure,
    SNP_R2_Outcome = snp_r2.outcome,
    R2_Ratio = r2_ratio,
    Causal_Direction = direction_correct,
    Steiger_Test = steiger_significant,
    Steiger_Pvalue = steiger_pval,
    Assumption_Status = assumption_met
  ) %>%
  # 按P值排序
  arrange(Steiger_Pvalue)

# 保存完整Steiger结果
write.csv(steiger_results, "Steiger_directionality_test_results.csv", row.names = FALSE)

# 筛选通过Steiger检验的结果 (方向正确且P<0.05)
valid_steiger <- steiger_results %>%
  filter(Causal_Direction == "TRUE" & Steiger_Test == "Significant")

final_results_table_genes <- unique(sig_results$gene_symbol)

valid_steiger <- valid_steiger %>% filter(Gene %in% final_results_table_genes)

# 保存有效结果
write.csv(valid_steiger, "Steiger_Valid_results.csv", row.names = FALSE)

