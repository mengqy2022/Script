#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Email     : 15877464851@163.com
# @Time      : 2023/9/27 10:00
# @File      : Symbionts_tools_RStudio.r

#-----------------------------------------------------------------------------------------
#  本脚本支持在RStudio中使用。                                                            |
#  绘制KAAS 5大分类模块代谢热图，通过完整度比对不同物种之间的代谢差异和相同点。            |
#  --------------------------------注意！！！！！---------------------------------------  |
#           只需要更改小部分路径，以及筛选的条件，就可以全脚本运行的到想要的结果          |
#           ###################################################################           |
#           #                            更改的地方                           #           |
#           ###################################################################           |
#  更改完成后，ctrl + shift + 回车,运行全部脚本即可的到结果                               |
#  或者ctrl + 回车,一行一行运行。                                                         |       
#-----------------------------------------------------------------------------------------

#  清空环境变量
rm(list = ls())

#  检查包是否安装  ----
#  设置镜像源
#  安装所需包
#  载入所需包

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

#  运行代码  ----
#############################################################
#  读取固定文件。                                           #
#  在不同的电脑，更改文件绝对目录。                         #
#  详细内容/etc/KEGG_General_Documentation_module/文件夹中。#
#############################################################
#  模块和K号对应信息
K_module <- read_delim("F:/database/KEGG/KEGG_General_Documentation_module/K_Module_2023_11_08.txt", delim = "\t", 
                       escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE, show_col_types = FALSE)
K_module_sum <- K_module %>% count_dt(X2) %>% rename_dt(m_id = X2, sum_n = n)

#  模块分类和名称
Module_name_class <- read_delim("F:/database/KEGG/KEGG_General_Documentation_module/Module_name_class.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE, show_col_types = FALSE)
Module_name_class <- Module_name_class %>% select_dt(B_class, m_id, m_name)

###########################################
#  读取之前包含文件名称的文件。           #
#  不同工作需要更改文件路径               #
###########################################
file_name <- read_csv("./sepice_name.txt", col_names = FALSE)$X1

#  读取所有的gene_anno文件
all_gene_ann <- list()
single_k <- list()
Class_K <- list()
Class_M <- list()
result <- list()

for (i in file_name) {
  ################################################
  #     在paste0中第一个双引号中加入 数据文件加  #
  #             例如：./data_genomes             #
  ################################################
  filename <- paste0("./database_all/",i,"/",i,"_gene_anno.txt")
  m <- paste0(i,"_completeness")
  if (!file.exists(filename)){
    next
  }
  #  读取所有的gene_anno文件
  all_gene_ann[[i]] <- read_delim(filename, 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
  #  挑选不重复的K号
  single_k[[i]] <-  all_gene_ann[[i]]$gene_id[!duplicated(all_gene_ann[[i]]$gene_id)]
  
  #  提取各个物种的K号
  Class_K[[i]] <- as.data.frame(single_k[[i]]) %>% mutate_dt(class = i)
  colnames(Class_K[[i]])[1] <- "k_ids"
  
  #  统计相关信息以及文件转换
  all_protein <- list_rbind(Class_K)
  
  if (!exists("k_module_map")) {
    # 映射模块上的K号（只能执行一次）
    K_module_map <- K_module %>% left_join_dt(all_protein,by = c("X1" = "k_ids")) %>% na.omit("class")
  }

  #  各个物种对应的模块信息
  Class_M[[i]] <- K_module_map %>% filter_dt(class == i) %>% select_dt(m_id = X2,class) %>% count_dt(m_id)
  colnames(Class_M[[i]])[2] <- i
  
  #  结果整合
  result[[i]] <- K_module_sum %>% left_join_dt(Class_M[[i]])
  m_id <- result[[i]]$m_id
  result[[i]]$m <- result[[i]][,3]/result[[i]]$sum_n
  result[[i]] <- result[[i]][,-1:-3]
  colnames(result[[i]])[1] <- m
}


result_all <- cbind.data.frame(result) %>% mutate_dt(m_id = m_id)

#  模块分类
result_all <- result_all %>% merge(Module_name_class) %>% 
  replace_na_dt(to = 0)
rm(Module_name_class)

#  模块大类的总数量
B_class <- result_all %>% count_dt(B_class)

#  热图绘制  ----
if (!dir.exists("./headmap_result_red_0_v_1.0")) {
  dir.create(path = "./headmap_result_red_0_v_1.0")
}

#  自定义颜色
#col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "green", "red"))

#  绘制热图
plot_heatmap <- function(data, title, filename) {
  row_labels <- structure(data$m_name, names = rownames(data))
  heatmap_data <- data %>% select_dt(-B_class, -m_name, -m_id) %>% as.matrix()
  rownames(heatmap_data) <- rownames(data)
  
  p <- Heatmap(heatmap_data, col = RColorBrewer::brewer.pal(9,"Blues"), name = "Completedness", column_title = title,
               rect_gp = gpar(col = "white", lty = 1, lwd = 2), row_title_side = "right",
               row_dend_width = unit(1, "cm"), cluster_columns = FALSE, column_names_rot = 90,
               row_labels = row_labels[rownames(heatmap_data)], row_names_max_width = max_text_width(row_labels),
               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 5),
               heatmap_legend_param = list(at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 0.25, 0.5, 0.75, 1),
                                           title = "Completedness", legend_height = unit(5, "cm"), title_position = "lefttop-rot"),
               column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               column_split = c(rep("a_symbio_gene", 11), rep("b_symbio_gene_draft", 4), rep("c_eup_gene", 5), rep("d_cili_gene", 8)))
  
  pdf(file = paste0("./headmap_result_red_0_v_1.0/", filename, ".pdf"), height = 12, width = 15)
  draw(p)
  dev.off()
}

#  分类筛选和绘图
categories <- c("Amino acid metabolism", "Carbohydrate metabolism", "Lipid metabolism", "Metabolism of cofactors and vitamins", "Energy metabolism")

for (cat in categories) {
  data <- result_all %>% filter_dt(B_class == cat)
  
  # 获取result中的所有列名并构建条件字符串
  valid_columns <- names(result)
  conditions <- paste0("(", paste0("data$", valid_columns, "_completeness != 0"), ")", collapse = " | ")
  
  # 使用条件过滤数据
  tryCatch({
    data <- data %>% filter_dt(eval(parse(text = conditions)))
    plot_heatmap(data, cat, str_replace_all(cat, " ", "_"))
  }, error = function(e) {
    message(sprintf("An error occurred while processing category '%s': %s", cat, e$message))
  })
}

#  合并剩余数据
result_1_7 <- result_all %>% filter_dt(B_class %in% c("Biosynthesis of other secondary metabolites", "Biosynthesis of terpenoids and polyketides", "Gene set", "Xenobiotics biodegradation", "Nucleotide metabolism", "Module set", "Glycan metabolism"))
str_input <- paste0("result_1_7$", names(result), "_completeness != 0 | ", collapse = '')
str_input <- substr(str_input, 0, nchar(str_input) - 3)
result_1_7 <- eval(parse(text = paste0("result_1_7 %>% filter_dt(", str_input, ")")))
plot_heatmap(result_1_7, "Resistance and other", "other")

