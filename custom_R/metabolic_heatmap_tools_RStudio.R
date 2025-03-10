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

base_packages=c("BiocManager",
                "optparse",
                "readr",
                "viridisLite",
                "purrr")

Bioconductor_packages=c("circlize",
                        "tidyfst",
                        "ComplexHeatmap",
                        "ggplot2")

# install packages in CRAN
for (pkg in base_packages){
  if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

# install packages in Bioconductor
for (pkg in Bioconductor_packages){
  if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}

# ignore all warnings, check if the packages could be libraried correctly.
for (pkg in c(Bioconductor_packages,base_packages)){
  require(pkg,character.only=T) 
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
K_module_sum <- K_module %>% count_dt(X2) %>% rename_dt(m_id = X2,sum_n =n )

#  模块分类和名称
Module_name_class <- read_delim("F:/database/KEGG/KEGG_General_Documentation_module/Module_name_class.txt", 
                                delim = "\t", escape_double = FALSE, 
                                trim_ws = TRUE, show_col_types = FALSE)
Module_name_class <- Module_name_class %>% select_dt(B_class,m_id,m_name)

###########################################
#  读取之前包含文件名称的文件。           #
#  不同工作需要更改文件路径               #
###########################################
file_name <- read_csv("./sepice_name.txt", col_names = FALSE)
file_name <- file_name$X1

#  读取所有的gene_anno文件
all_gene_ann <- list()
#  挑选不重复的K号
single_k <- list()
#  提取各个物种的K号
Class_K <- list()
#  各个物种对应的模块信息
Class_M <- list()
#  结果整合
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

#  输出文件
#write.table(K_module_map, file = "tools_RStudio_映射模块上的K号.txt",sep = "\t", quote = F,col.names = F,row.names = F)
#write.table(all_protein, file = "tools_RStudio_唯一K号统计.txt",sep = "\t", quote = F,col.names = F,row.names = F)
rm(all_gene_ann,all_protein,Class_K,Class_M,K_module_map,single_k,i,m)

result_all <- cbind.data.frame(result) %>% mutate_dt(m_id = m_id)

#  模块分类
result_all <- result_all %>% merge(Module_name_class) %>% 
  replace_na_dt(to = 0)
rm(Module_name_class)

#  输出文件
library(xlsx)
#write.xlsx(result_all, "result_all.xlsx",row.names = F)
#write.table(result_all, file = "tools_RStudio_result_all.txt" ,sep = "\t", quote = F,col.names = T,row.names = F)

#  手动更改 详细读README
library(readxl)
result_all <- read_excel("20240731_heatmap_data_all.xlsx")

#  模块大类的总数量
B_class <- result_all %>%  count_dt(B_class)

#  热图绘制  ----
######################################################
#  创建目录，可以根据自己要去更改，但是要全部替换。  #
######################################################
if (!file.exists("./headmap_result_red_0")) {
  dir.create(path = "./headmap_result_red_0")
}

#  自定义颜色
#col_fun = colorRamp2(c(0, 0.5, 1), c( "white","green","red"))
#col_fun = colorRamp2(c(0, 0.2, 0.6, 0.8, 1), c("#def5e5", "#74d1b2","#3591a7","#37699e","#3c3365"))
#col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("#C16D58", "#ECD0B4","#3591a7","#37699e","#3c3365"))

# Amino acid biosynthesis ------
#  输出多个筛选值，由于数量太多了，自动输出，复制粘贴
#  将都为0的代谢模块去除
str_input <- paste0("Amino_acid$",names(result),"_completeness"," != 0 | ",collapse ='')
str_input <- substr(str_input,0,nchar(str_input)-3) 
str_input

Amino_acid <- result_all %>% filter_dt(B_class == "Amino acid metabolism")
######################################################
# 复制str_input内容，到下方filter_dt中，不要带双引号 #
######################################################
Amino_acid <- Amino_acid %>% filter_dt(Amino_acid$Polynucleobacter_necessarius_Eae1_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_Eae2_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_Eae3_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_Eae5_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_Eco1_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_Eda1_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_Ewo1_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_Fsp1.4_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_STIR1_completeness != 0 | Amino_acid$Polynucleobacter_necessarius_amieti_completeness != 0 | Amino_acid$Fl_phosphoraccumulans_VTN8_completeness != 0 | Amino_acid$Ca_Devosia_euplotis_LIV5_completeness != 0 | Amino_acid$Ca_Devosia_symbiotica_Na2_completeness != 0 | Amino_acid$Ca_Protistobacter_heckmanni_Eae6_completeness != 0 | Amino_acid$Ca_Protistobacter_heckmanni_POH1_completeness != 0 | Amino_acid$Euplotes_octiocarinztus_completeness != 0 | Amino_acid$Euplotes_aediculatus_completeness != 0 | Amino_acid$Euplotes_amieti_completeness != 0 | Amino_acid$Euplotes_woodruffi_completeness != 0 | Amino_acid$Euplotes_vannus_completeness != 0 | Amino_acid$Tetrahymena_borealis_completeness != 0 | Amino_acid$Tetrahymena_elliotti_completeness != 0 | Amino_acid$Tetrahymena_malaccensis_completeness != 0 | Amino_acid$Tetrahymena_thermophila_completeness != 0 | Amino_acid$Paramecium_biaurelia_completeness != 0 | Amino_acid$Paramecium_bursaria_completeness != 0 | Amino_acid$Paramecium_caudatum_completeness != 0 | Amino_acid$Paramecium_tetraurelia_completeness != 0)

#  20240731 修改  详细内容看READNE
# library(xlsx)
# write.xlsx(Amino_acid, "Amino_acid.xlsx",row.names = F)
# library(readxl)
# Amino_acid <- read_excel("Amino_acid.xlsx")

#  自定义列名
row_labels_Amino_acid <- structure(Amino_acid$m_name, names = rownames(Amino_acid))
#  绘图数据
metabolize_amino <- Amino_acid %>% select_dt(-B_class,-m_name,-m_id)
metabolize_amino <- as.matrix(metabolize_amino)
#  赋予列因子信息
rownames(metabolize_amino) = rownames(Amino_acid)

#  绘制热图
p1 <- Heatmap(metabolize_amino, 
              #col = col_fun,   
              #col = RColorBrewer::brewer.pal(9,"Blues"),
              #col = RColorBrewer::brewer.pal(9,"BuPu"),
              #viridis(alpha = 0.5, 100 , begin = 1, end = 0,option = "H"),
              mako(100, begin = 1, end = 0),
              #plasma(6, begin = 1, end = 0),
              #cividis(6, begin = 0, end = 1),
              name = "Completedness",
              column_title = "Amino acid biosynthesis",
              rect_gp = gpar(col = "white", lty = 1, lwd = 2),
              row_title_side = "right",
              row_dend_width = unit(1, "cm"),
              # 列不聚类
              cluster_columns = F,
              #column_dend_height = unit(0.5, "cm"),
              #row_names_centered = TRUE, 
              #column_names_centered = TRUE,
              column_names_rot = 90,
              row_labels = row_labels_Amino_acid[rownames(metabolize_amino)],
              row_names_max_width = max_text_width(row_labels_Amino_acid),
              row_names_gp = gpar(fontsize = 8),
              column_names_gp= gpar(fontsize = 5),
              heatmap_legend_param = list(
                at = c(0,0.25,0.5,0.75,1),
                labels = c(0,0.25,0.5, 0.75,1),
                title = "Completedness",
                legend_height = unit(5, "cm"),
                title_position = "lefttop-rot"),
              #  标题字体调整
              column_title_gp = gpar(fontsize = 20, fontface = "bold"),
              #  分割热图
              column_split = c(rep(("a_symbio_gene"), 11),rep(("b_symbio_gene_draft"), 4),rep(("c_eup_gene"), 5),rep(("d_cili_gene"), 8))
)

#  Carbohydrates and lipid metabolism ------
#  输出多个筛选值，由于数量太多了，自动输出，复制粘贴
str_input <- paste0("Carbohydrate_Lipid$",names(result),"_completeness"," != 0 | ",collapse ='')
str_input <- substr(str_input,0,nchar(str_input)-3)
str_input

Carbohydrate <- result_all[which(result_all$B_class == "Carbohydrate metabolism"),]
Lipid <- result_all[which(result_all$B_class == "Lipid metabolism"),]
# Carbohydrate <- result_all %>% filter_dt(B_class == "Carbohydrate metabolism")
# Lipid <- result_all %>% filter_dt(B_class == "Lipid metabolism")
Carbohydrate_Lipid <- rbind(Carbohydrate,Lipid)
#####################################################
#  复制str_input内容，到下方filter_dt中，不要带引号 #
#####################################################
Carbohydrate_Lipid <- Carbohydrate_Lipid %>% filter_dt(Carbohydrate_Lipid$Polynucleobacter_necessarius_Eae1_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_Eae2_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_Eae3_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_Eae5_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_Eco1_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_Eda1_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_Ewo1_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_Fsp1.4_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_STIR1_completeness != 0 | Carbohydrate_Lipid$Polynucleobacter_necessarius_amieti_completeness != 0 | Carbohydrate_Lipid$Fl_phosphoraccumulans_VTN8_completeness != 0 | Carbohydrate_Lipid$Ca_Devosia_euplotis_LIV5_completeness != 0 | Carbohydrate_Lipid$Ca_Devosia_symbiotica_Na2_completeness != 0 | Carbohydrate_Lipid$Ca_Protistobacter_heckmanni_Eae6_completeness != 0 | Carbohydrate_Lipid$Ca_Protistobacter_heckmanni_POH1_completeness != 0 | Carbohydrate_Lipid$Euplotes_octiocarinztus_completeness != 0 | Carbohydrate_Lipid$Euplotes_aediculatus_completeness != 0 | Carbohydrate_Lipid$Euplotes_amieti_completeness != 0 | Carbohydrate_Lipid$Euplotes_woodruffi_completeness != 0 | Carbohydrate_Lipid$Euplotes_vannus_completeness != 0 | Carbohydrate_Lipid$Tetrahymena_borealis_completeness != 0 | Carbohydrate_Lipid$Tetrahymena_elliotti_completeness != 0 | Carbohydrate_Lipid$Tetrahymena_malaccensis_completeness != 0 | Carbohydrate_Lipid$Tetrahymena_thermophila_completeness != 0 | Carbohydrate_Lipid$Paramecium_biaurelia_completeness != 0 | Carbohydrate_Lipid$Paramecium_bursaria_completeness != 0 | Carbohydrate_Lipid$Paramecium_caudatum_completeness != 0 | Carbohydrate_Lipid$Paramecium_tetraurelia_completeness != 0)

#  绘图数据
metabolize_Car_lipid <- Carbohydrate_Lipid %>% select_dt(-B_class,-m_name,-m_id)
metabolize_Car_lipid  <- as.data.frame(metabolize_Car_lipid)

#  自定义列名
row_labels_Car_lipid = structure(Carbohydrate_Lipid$m_name, names =rownames(Carbohydrate_Lipid))
rm(Carbohydrate,Lipid)

metabolize_Car_lipid_matrix <- metabolize_Car_lipid %>% as.matrix()

#  赋予列因子信息
rownames(metabolize_Car_lipid_matrix) = rownames(metabolize_Car_lipid)

p2 <- Heatmap(metabolize_Car_lipid_matrix,
        #viridis(alpha = 0.5,10, begin = 1, end = 0,option = "G"),
        mako(100, begin = 1, end = 0),
        #col = col_fun,
        name = "Completedness",
        column_title = "Carbohydrates and lipid metabolism",
        rect_gp = gpar(col = "white", lty = 1, lwd = 2),
        row_title_side = "right",
        row_dend_width = unit(1, "cm"),
        # 列不聚类
        cluster_columns = F,
        #column_dend_height = unit(0.5, "cm"),
        #row_names_centered = TRUE, 
        #column_names_centered = TRUE,
        column_names_rot = 90,
        row_labels = row_labels_Car_lipid[rownames(metabolize_Car_lipid)],
        row_names_max_width = max_text_width(row_labels_Car_lipid),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp= gpar(fontsize = 5),
        heatmap_legend_param = list(
          at = c(0,0.25,0.5,0.75,1),
          labels = c(0,0.25,0.5, 0.75,1),
          title = "Completedness",
          legend_height = unit(5, "cm"),
          title_position = "lefttop-rot"),
        #  标题字体调整
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        #  分割热图
        column_split = c(rep(("a_symbio_gene"), 11),rep(("b_symbio_gene_draft"), 4),rep(("c_eup_gene"), 5),rep(("d_cili_gene"), 8))
)

#  Cofactor and vitamin biosynthesis -------
#  输出多个筛选值，由于数量太多了，自动输出，复制粘贴
str_input <- paste0("cofa_vita$",names(result),"_completeness"," != 0 | ",collapse ='')
str_input <- substr(str_input,0,nchar(str_input)-3)
str_input

cofa_vita <- result_all %>% filter_dt(B_class == "Metabolism of cofactors and vitamins")
#####################################################
#  复制str_input内容，到下方filter_dt中，不要带引号 #
#####################################################
cofa_vita <- cofa_vita %>% filter_dt(cofa_vita$Polynucleobacter_necessarius_Eae1_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_Eae2_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_Eae3_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_Eae5_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_Eco1_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_Eda1_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_Ewo1_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_Fsp1.4_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_STIR1_completeness != 0 | cofa_vita$Polynucleobacter_necessarius_amieti_completeness != 0 | cofa_vita$Fl_phosphoraccumulans_VTN8_completeness != 0 | cofa_vita$Ca_Devosia_euplotis_LIV5_completeness != 0 | cofa_vita$Ca_Devosia_symbiotica_Na2_completeness != 0 | cofa_vita$Ca_Protistobacter_heckmanni_Eae6_completeness != 0 | cofa_vita$Ca_Protistobacter_heckmanni_POH1_completeness != 0 | cofa_vita$Euplotes_octiocarinztus_completeness != 0 | cofa_vita$Euplotes_aediculatus_completeness != 0 | cofa_vita$Euplotes_amieti_completeness != 0 | cofa_vita$Euplotes_woodruffi_completeness != 0 | cofa_vita$Euplotes_vannus_completeness != 0 | cofa_vita$Tetrahymena_borealis_completeness != 0 | cofa_vita$Tetrahymena_elliotti_completeness != 0 | cofa_vita$Tetrahymena_malaccensis_completeness != 0 | cofa_vita$Tetrahymena_thermophila_completeness != 0 | cofa_vita$Paramecium_biaurelia_completeness != 0 | cofa_vita$Paramecium_bursaria_completeness != 0 | cofa_vita$Paramecium_caudatum_completeness != 0 | cofa_vita$Paramecium_tetraurelia_completeness != 0)

#  自定义列名
row_labels_cofa_vita <- structure(cofa_vita$m_name,names = rownames(cofa_vita))
#  绘图数据
metabolize_cofa_vita <- cofa_vita %>% select_dt(-B_class,-m_name,-m_id) 
metabolize_cofa_vita <- as.matrix(metabolize_cofa_vita)
#  赋予列因子信息
rownames(metabolize_cofa_vita) <- rownames(cofa_vita)

p3 <- Heatmap(metabolize_cofa_vita,
        mako(100, begin = 1, end = 0),
        #viridis(alpha = 0.5,10, begin = 1, end = 0,option = "G"),
        #col = col_fun,
        name = "Completedness",
        column_title = "Cofactor and vitamin biosynthesis",
        rect_gp = gpar(col = "white", lty = 1, lwd = 2),
        row_title_side = "right",
        row_dend_width = unit(1, "cm"),
        # 列不聚类
        cluster_columns = F,
        #column_dend_height = unit(0.5, "cm"),
        #row_names_centered = TRUE, 
        #column_names_centered = TRUE,
        column_names_rot = 90,
        row_labels = row_labels_cofa_vita[rownames(metabolize_cofa_vita)],
        row_names_max_width = max_text_width(row_labels_cofa_vita),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp= gpar(fontsize = 5),
        heatmap_legend_param = list(
          at = c(0,0.25,0.5,0.75,1),
          labels = c(0,0.25,0.5, 0.75,1),
          title = "Completedness",
          legend_height = unit(5, "cm"),
          title_position = "lefttop-rot"),
        #  标题字体调整
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        #  分割热图
        #column_split = c(rep(("a_genomics"), 4),rep(("b_draft"), 5))
        column_split = c(rep(("a_symbio_gene"), 11),rep(("b_symbio_gene_draft"), 4),rep(("c_eup_gene"), 5),rep(("d_cili_gene"), 8))
)

#  Energy metabolism ------
#  输出多个筛选值，由于数量太多了，自动输出，复制粘贴
str_input <- paste0("Energy$",names(result),"_completeness"," != 0 | ",collapse ='')
str_input <- substr(str_input,0,nchar(str_input)-3)
str_input
#  数据
Energy <- result_all %>% filter_dt(B_class == "Energy metabolism") 
#####################################################
#  复制str_input内容，到下方filter_dt中，不要带引号 #
#####################################################
Energy <- Energy %>% filter_dt(Energy$Polynucleobacter_necessarius_Eae1_completeness != 0 | Energy$Polynucleobacter_necessarius_Eae2_completeness != 0 | Energy$Polynucleobacter_necessarius_Eae3_completeness != 0 | Energy$Polynucleobacter_necessarius_Eae5_completeness != 0 | Energy$Polynucleobacter_necessarius_Eco1_completeness != 0 | Energy$Polynucleobacter_necessarius_Eda1_completeness != 0 | Energy$Polynucleobacter_necessarius_Ewo1_completeness != 0 | Energy$Polynucleobacter_necessarius_Fsp1.4_completeness != 0 | Energy$Polynucleobacter_necessarius_STIR1_completeness != 0 | Energy$Polynucleobacter_necessarius_amieti_completeness != 0 | Energy$Fl_phosphoraccumulans_VTN8_completeness != 0 | Energy$Ca_Devosia_euplotis_LIV5_completeness != 0 | Energy$Ca_Devosia_symbiotica_Na2_completeness != 0 | Energy$Ca_Protistobacter_heckmanni_Eae6_completeness != 0 | Energy$Ca_Protistobacter_heckmanni_POH1_completeness != 0 | Energy$Euplotes_octiocarinztus_completeness != 0 | Energy$Euplotes_aediculatus_completeness != 0 | Energy$Euplotes_amieti_completeness != 0 | Energy$Euplotes_woodruffi_completeness != 0 | Energy$Euplotes_vannus_completeness != 0 | Energy$Tetrahymena_borealis_completeness != 0 | Energy$Tetrahymena_elliotti_completeness != 0 | Energy$Tetrahymena_malaccensis_completeness != 0 | Energy$Tetrahymena_thermophila_completeness != 0 | Energy$Paramecium_biaurelia_completeness != 0 | Energy$Paramecium_bursaria_completeness != 0 | Energy$Paramecium_caudatum_completeness != 0 | Energy$Paramecium_tetraurelia_completeness != 0)

#  自定义列名
row_labels_Energy = structure(Energy$m_name, names = rownames(Energy))
#  绘图数据
metabolize_Energy <- Energy %>% select_dt(-B_class,-m_name,-m_id) 
metabolize_Energy <- as.matrix(metabolize_Energy)
#  赋予列因子信息
rownames(metabolize_Energy) <- rownames(Energy)

p4 <- Heatmap(metabolize_Energy,
        mako(100, begin = 1, end = 0),
        #viridis(alpha = 0.5,10, begin = 1, end = 0,option = "G"),
        #col = col_fun,
        name = "Completedness",
        column_title = "Energy metabolism",
        rect_gp = gpar(col = "white", lty = 1, lwd = 2),
        row_title_side = "right",
        row_dend_width = unit(1, "cm"),
        # 列不聚类
        cluster_columns = F,
        #column_dend_height = unit(0.5, "cm"),
        #row_names_centered = TRUE, 
        #column_names_centered = TRUE,
        column_names_rot = 90,
        row_labels = row_labels_Energy[rownames(metabolize_Energy)],
        row_names_max_width = max_text_width(row_labels_Energy),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp= gpar(fontsize = 5),
        heatmap_legend_param = list(
          at = c(0,0.25,0.5,0.75,1),
          labels = c(0,0.25,0.5, 0.75,1),
          title = "Completedness",
          legend_height = unit(5, "cm"),
          title_position = "lefttop-rot"),
        #  标题字体调整
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        #  分割热图
        #column_split = c(rep(("a_genomics"), 4),rep(("b_draft"), 5))
        column_split = c(rep(("a_symbio_gene"), 11),rep(("b_symbio_gene_draft"), 4),rep(("c_eup_gene"), 5),rep(("d_cili_gene"), 8))
)

#  Resistance and other ------
#  输出多个筛选值，由于数量太多了，自动输出，复制粘贴
str_input <- paste0("result_1_7$",names(result),"_completeness"," != 0 | ",collapse ='')
str_input <- substr(str_input,0,nchar(str_input)-3)
str_input

#  合并剩余数据
result_1 <- result_all %>% filter_dt(B_class == "Biosynthesis of other secondary metabolites")
result_2 <- result_all %>% filter_dt(B_class == "Biosynthesis of terpenoids and polyketides")
result_3 <- result_all %>% filter_dt(B_class == "Gene set")
result_4 <- result_all %>% filter_dt(B_class == "Xenobiotics biodegradation")
result_5 <- result_all %>% filter_dt(B_class == "Nucleotide metabolism")
result_6 <- result_all %>% filter_dt(B_class == "Module set")
result_7 <- result_all %>% filter_dt(B_class == "Glycan metabolism")
result_1_7 <- rbind(result_1,result_2,result_3,result_4,result_5,result_6,result_7)
#####################################################
#  复制str_input内容，到下方filter_dt中，不要带引号 #
#####################################################
result_1_7 <- result_1_7 %>% filter_dt(result_1_7$Polynucleobacter_necessarius_Eae1_completeness != 0 | result_1_7$Polynucleobacter_necessarius_Eae2_completeness != 0 | result_1_7$Polynucleobacter_necessarius_Eae3_completeness != 0 | result_1_7$Polynucleobacter_necessarius_Eae5_completeness != 0 | result_1_7$Polynucleobacter_necessarius_Eco1_completeness != 0 | result_1_7$Polynucleobacter_necessarius_Eda1_completeness != 0 | result_1_7$Polynucleobacter_necessarius_Ewo1_completeness != 0 | result_1_7$Polynucleobacter_necessarius_Fsp1.4_completeness != 0 | result_1_7$Polynucleobacter_necessarius_STIR1_completeness != 0 | result_1_7$Polynucleobacter_necessarius_amieti_completeness != 0 | result_1_7$Fl_phosphoraccumulans_VTN8_completeness != 0 | result_1_7$Ca_Devosia_euplotis_LIV5_completeness != 0 | result_1_7$Ca_Devosia_symbiotica_Na2_completeness != 0 | result_1_7$Ca_Protistobacter_heckmanni_Eae6_completeness != 0 | result_1_7$Ca_Protistobacter_heckmanni_POH1_completeness != 0 | result_1_7$Euplotes_octiocarinztus_completeness != 0 | result_1_7$Euplotes_aediculatus_completeness != 0 | result_1_7$Euplotes_amieti_completeness != 0 | result_1_7$Euplotes_woodruffi_completeness != 0 | result_1_7$Euplotes_vannus_completeness != 0 | result_1_7$Tetrahymena_borealis_completeness != 0 | result_1_7$Tetrahymena_elliotti_completeness != 0 | result_1_7$Tetrahymena_malaccensis_completeness != 0 | result_1_7$Tetrahymena_thermophila_completeness != 0 | result_1_7$Paramecium_biaurelia_completeness != 0 | result_1_7$Paramecium_bursaria_completeness != 0 | result_1_7$Paramecium_caudatum_completeness != 0 | result_1_7$Paramecium_tetraurelia_completeness != 0)

#  自定义列名
row_labels_other = structure(result_1_7$m_name, names = rownames(result_1_7))
rm(result_1,result_2,result_3,result_4,result_5,result_6,result_7)
#  绘图数据
metabolize_other <- result_1_7 %>% select_dt(-B_class,-m_name,-m_id) 
metabolize_other <- as.matrix(metabolize_other)
#  行注释信息
rownames(metabolize_other) = rownames(result_1_7)

p5 <- Heatmap(metabolize_other,
        mako(100, begin = 1, end = 0),
        #viridis(alpha = 0.5,10, begin = 1, end = 0,option = "G"),
        #col = col_fun,
        name = "Completedness",
        column_title = "Resistance and other",
        rect_gp = gpar(col = "white", lty = 1, lwd = 2),
        row_title_side = "right",
        row_dend_width = unit(1, "cm"),
        # 列不聚类
        cluster_columns = F,
        #column_dend_height = unit(0.5, "cm"),
        #row_names_centered = TRUE, 
        #column_names_centered = TRUE,
        column_names_rot = 90,
        row_labels = row_labels_other[rownames(metabolize_other)],
        row_names_max_width = max_text_width(row_labels_other),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp= gpar(fontsize = 5),
        heatmap_legend_param = list(
          at = c(0,0.25,0.5,0.75,1),
          labels = c(0,0.25,0.5, 0.75,1),
          title = "Completedness",
          legend_height = unit(5, "cm"),
          title_position = "lefttop-rot"),
        #  标题字体调整
        column_title_gp = gpar(fontsize = 20, fontface = "bold"),
        #  分割热图
        column_split = c(rep(("a_symbio_gene"), 11),rep(("b_symbio_gene_draft"), 4),rep(("c_eup_gene"), 5),rep(("d_cili_gene"), 8))
)

#  保存pdf ------
##############################################
# 根据物种的多少适当调整height和width的大小。#
#        为了图片更美观调整图片的大小        #
#            可以保存多种图片格式            #
#        例如：png()  tiff()  svg()等        #
##############################################
pdf(file= "./headmap_result_red_0/Amina.pdf" , height=12,width=15)#新建一个PDF文件，设置名称、宽高及字体等
p1
dev.off()

pdf(file='./headmap_result_red_0/Car_lipid.pdf', height=15,width=15)#新建一个PDF文件，设置名称、宽高及字体等
p2
dev.off()

pdf(file='./headmap_result_red_0/Cofa_vite.pdf', height=11,width=15)#新建一个PDF文件，设置名称、宽高及字体等
p3
dev.off()

pdf(file='./headmap_result_red_0/Energy.pdf', height=12,width=15)#新建一个PDF文件，设置名称、宽高及字体等
p4
dev.off()

pdf(file='./headmap_result_red_0/Other.pdf', height=15,width=15)#新建一个PDF文件，设置名称、宽高及字体等
p5
dev.off()