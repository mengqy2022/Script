#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Email     : 15877464851@163.com
# @Time      : 2024/09/10 10:00

# 加载必要的库
if (!require(argparse)) {
  install.packages("argparse", repos = "https://mirrors.ustc.edu.cn/CRAN/")
  library(argparse)
}

# 创建解析器
parser <- ArgumentParser(description = 'Metabolic pathway heat map')
parser$add_argument("--k_module_file", required=TRUE, help="Path to the K_Module file")
parser$add_argument("--module_name_class_file", required=TRUE, help="Path to the Module_name_class file")
parser$add_argument("--sepice_name_file", required=TRUE, help="Path to the sepice_name file")
parser$add_argument("--base_dir", required=TRUE, help="Base directory for data files")
parser$add_argument("--output_file", default="heatmap_result.pdf", help="Output file for the Heatmap")
parser$add_argument("--color_scheme", choices=c("Blues", "Greens", "Reds"), default="Blues", help="Color scheme for the Heatmap")
parser$add_argument("--categories", nargs='+', help="List of categories to include in the output, example: [Amino acid metabolism], [Carbohydrate metabolism], [Lipid metabolism], [Metabolism of cofactors and vitamins], [Energy metabolism]")
parser$add_argument("--include_result_1_7", action='store_true', help="Include result_1_7 in the output")
parser$add_argument("--output_all", action='store_true', help="Output all categories and result_1_7")

# 解析参数
args <- parser$parse_args()

# 设置镜像源并加载必要的包
options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

base_packages <- c("BiocManager", "optparse", "readr", "viridisLite", "purrr", "readxl")
Bioconductor_packages <- c("circlize", "tidyfst", "ComplexHeatmap", "ggplot2", "tidyverse")

# 安装和加载CRAN包
for (pkg in base_packages) {
  if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
    install.packages(pkg, ask = FALSE, update = FALSE)
    require(pkg, character.only = TRUE)
  }
}

# 安装和加载Bioconductor包
for (pkg in Bioconductor_packages) {
  if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    require(pkg, character.only = TRUE)
  }
}

# 读取固定文件
K_module <- read_delim(args$k_module_file, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, show_col_types = FALSE)
K_module_sum <- K_module %>% count_dt(X2) %>% rename_dt(m_id = X2, sum_n = n)
Module_name_class <- read_delim(args$module_name_class_file, delim = "\t", escape_double = FALSE, trim_ws = TRUE, show_col_types = FALSE)
Module_name_class <- Module_name_class %>% select_dt(B_class, m_id, m_name)
file_name <- read_csv(args$sepice_name_file, col_names = FALSE)$X1

# 读取所有的gene_anno文件
all_gene_ann <- list()
single_k <- list()
Class_K <- list()
Class_M <- list()
result <- list()

for (i in file_name) {
  filename <- paste0(args$base_dir, "/", i, "/", i, "_gene_anno.txt")
  m <- paste0(i, "_completeness")
  
  if (!file.exists(filename)) {
    message(sprintf("File not found: %s", filename))
    next
  }
  
  all_gene_ann[[i]] <- read_delim(filename, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  single_k[[i]] <- all_gene_ann[[i]]$gene_id[!duplicated(all_gene_ann[[i]]$gene_id)]
  Class_K[[i]] <- as.data.frame(single_k[[i]]) %>% mutate_dt(class = i)
  colnames(Class_K[[i]])[1] <- "k_ids"
  all_protein <- list_rbind(Class_K)
  
  if (!exists("k_module_map")) {
    K_module_map <- K_module %>% left_join_dt(all_protein, by = c("X1" = "k_ids")) %>% na.omit("class")
  }
  
  Class_M[[i]] <- K_module_map %>% filter_dt(class == i) %>% select_dt(m_id = X2, class) %>% count_dt(m_id)
  colnames(Class_M[[i]])[2] <- i
  result[[i]] <- K_module_sum %>% left_join_dt(Class_M[[i]])
  m_id <- result[[i]]$m_id
  result[[i]]$m <- result[[i]][,3] / result[[i]]$sum_n
  result[[i]] <- result[[i]][,-1:-3]
  colnames(result[[i]])[1] <- m
}

result_all <- cbind.data.frame(result) %>% mutate_dt(m_id = m_id)
result_all <- result_all %>% merge(Module_name_class) %>% replace_na_dt(to = 0)
rm(Module_name_class)

# 自定义绘制Heatmap的函数
plot_heatmap <- function(data, title) {
  row_labels <- structure(data$m_name, names = rownames(data))
  heatmap_data <- data %>% select_dt(-B_class, -m_name, -m_id) %>% as.matrix()
  rownames(heatmap_data) <- rownames(data)
  
  p <- Heatmap(heatmap_data, col = RColorBrewer::brewer.pal(9, args$color_scheme), name = "Completedness", column_title = title,
               rect_gp = gpar(col = "white", lty = 1, lwd = 2), row_title_side = "right",
               row_dend_width = unit(1, "cm"), cluster_columns = FALSE, column_names_rot = 90,
               row_labels = row_labels[rownames(heatmap_data)], row_names_max_width = max_text_width(row_labels),
               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 5),
               heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1),
                                           title = "Completedness", legend_height = unit(5, "cm"), title_position = "lefttop-rot"),
               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               column_split = c(rep("a_symbio_gene", 11), rep("b_symbio_gene_draft", 4), rep("c_eup_gene", 5), rep("d_cili_gene", 8)))
  
  pdf(file = args$output_file, height = 12, width = 15)
  draw(p)
  dev.off()
  
  message(sprintf("Output saved to %s", args$output_file))
}

# 处理选择的类别
if (args$output_all) {
  args$categories <- unique(result_all$B_class)
  args$include_result_1_7 <- TRUE
}

if (!is.null(args$categories)) {
  for (cat in args$categories) {
    data <- result_all %>% filter_dt(B_class == cat)
    valid_columns <- names(result)
    conditions <- paste0("(", paste0("data$", valid_columns, "_completeness != 0"), ")", collapse = " | ")
    
    # 这些函数提供了处理异常情况（包括错误和警告）的机制。
    tryCatch({
      data <- data %>% filter_dt(eval(parse(text = conditions)))
      plot_heatmap(data, cat)
    }, error = function(e) {
      message(sprintf("An error occurred while processing category '%s': %s", cat, e$message))
    })
  }
}

# 如果用户选择包含result_1_7
if (args$include_result_1_7) {
  result_1_7 <- result_all %>% filter_dt(B_class %in% c(
    "Biosynthesis of other secondary metabolites", "Biosynthesis of terpenoids and polyketides", 
    "Gene set", "Xenobiotics biodegradation", "Nucleotide metabolism", "Module set", "Glycan metabolism"
  ))
  valid_columns <- names(result)
  conditions <- paste0("(", paste0("result_1_7$", valid_columns, "_completeness != 0"), ")", collapse = " | ")
  
  tryCatch({
    result_1_7 <- result_1_7 %>% filter_dt(eval(parse(text = conditions)))
    plot_heatmap(result_1_7, "Resistance and other")
  }, error = function(e) {
    message(sprintf("An error occurred while processing result_1_7: %s", e$message))
  })
}
