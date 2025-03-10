#!/usr/bin/env Rscript

# TNM 是癌症分期系统中用于描述肿瘤的三个关键方面：
# T（Tumor）：原发肿瘤的大小和范围。
# 示例：T0（无原发肿瘤），T1-T4（根据肿瘤大小和范围划分）。
# TX 表示无法评估。
# N（Node）：区域淋巴结的受累情况。
# 示例：N0（无淋巴结受累），N1-N3（根据受累淋巴结的数量和范围划分）。
# NX 表示无法评估。
# M（Metastasis）：是否存在远处转移。
# 示例：M0（无远处转移），M1（有远处转移）。
# MX 表示无法评估。

library(argparse)
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)

# 创建参数解析器对象
parser <- ArgumentParser(description="TCGA数据下载和整理脚本",
                         epilog="\nRscript TCGA_downtools.r -d your_data_path -c TCGA_RNA_data -l TCGA_clinical_data")

# 添加参数
parser$add_argument("-d", "--data_path", help="保存数据的路径", type="character", default="your_data_path")
parser$add_argument("-c", "--count_tpm_dir", help="保存Counts和TPM数据的目录", type="character", default="TCGA_RNA_data")
parser$add_argument("-l", "--clinical_dir", help="保存临床数据的目录", type="character", default="TCGA_clinical_data")

# 解析命令行参数
args <- parser$parse_args()

# 如果需要安装BiocManager和TCGAbiolinks，可以取消注释以下两行
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("TCGAbiolinks")

# -----------------------------# TCGA 数据批量下载脚本# -----------------------------
rm(list = ls())
setwd(args$data_path)

# 创建目录
if (!dir.exists(args$count_tpm_dir)) { dir.create(args$count_tpm_dir) }
if (!dir.exists(args$clinical_dir)) { dir.create(args$clinical_dir) }

# 获取TCGA项目列表
projects <- getGDCprojects()
projects <- projects %>% as.data.frame() %>% select(project_id, tumor) %>% filter(grepl(pattern = "^TCGA", project_id))
projects <- projects[order(projects$tumor),]
projects

# 遍历每一个TCGA项目
for (i in 1:nrow(projects)) {
  project_id <- projects$project_id[i]
  tryCatch({
    # 查询基因表达数据（Counts和TPM）
    query_exp <- GDCquery(
      project = project_id,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    
    # 下载数据
    GDCdownload(query_exp)
    
    # 整理数据，生成SummarizedExperiment对象
    exprSet <- GDCprepare(query = query_exp)
    
    # 提取Counts和TPM数据
    exp_Count <- assay(exprSet, "unstranded")  # Counts数据
    exp_TPM <- assay(exprSet, "tpm_unstrand")  # TPM数据
    
    # 提取基因注释信息并合并基因ID和基因名
    ann <- rowRanges(exprSet)  # 提取基因注释信息
    ann <- as.data.frame(ann)
    
    # 保存RData
    save(exp_Count, exp_TPM, file = file.path(args$count_tpm_dir, paste0(project_id, ".Rdata")))
  }, error = function(e) {
    # 如果出错，跳过当前项目并打印错误信息
    message("Error in project ", project_id, ": ", e$message)
    return(NULL)  # 跳到下一个项目
  })
}

# 提取基因注释信息并合并基因ID和基因名
ann <- rowRanges(exprSet)  # 提取基因注释信息
ann <- as.data.frame(ann)
ann <- ann[rownames(exp_Count), ]  # 保留Counts数据中的基因
write.csv(ann, file = "tgca_gene_annotation.csv", row.names = FALSE)
message("所有项目已处理完成。")

# 基因注释 ----------------------------------------------------------------------------{
setwd(args$data_path)

# 基因注释文件
anno <- read.csv("tgca_gene_annotation.csv")
anno <- anno %>% dplyr::select(gene_id, gene_name)

# 获取RData
rdata_files <- list.files(path = args$count_tpm_dir, pattern = "*.Rdata", full.names = TRUE)

# 遍历每个 RData 文件
for (rdata_file in rdata_files) {
  load(rdata_file)
  
  # 提取项目名称作为文件名的一部分
  project_name <- tools::file_path_sans_ext(basename(rdata_file))
  
  # 注释，去重
  exp_TPM_annotated <- exp_TPM %>% 
    as.data.frame() %>%
    dplyr::mutate(gene_id = rownames(.)) %>%  # 添加 gene_id 列（原行名）
    dplyr::left_join(anno, by = "gene_id") %>%  # 合并基因注释数据框（按 gene_id 列）
    dplyr::select(gene_name, gene_id, everything()) %>%  # 将 gene_name 和 gene_id 列放在最前面
    dplyr::arrange(desc(rowMeans(.[, -c(1, 2)]))) %>%  # 按平均表达量降序排列，忽略 gene_name 和 gene_id 列
    dplyr::filter(!duplicated(gene_name))  # 去重，保留表达量最高的基因
  
  # 保存合并后的数据为新的 RData 文件
  save(exp_Count, exp_TPM_annotated, file = paste0("TCGA_RNA_data_anno/", project_name, ".Rdata"))
  message(paste("Processed and saved:", project_name))
}
# }

##### 2 下载indexed clinical 信息----------------------------------------------------------
clinical_dir <- args$clinical_dir

# 获取项目ID
projects <- getGDCprojects()
projects <- projects %>% as.data.frame() %>% dplyr::select(project_id, tumor) %>% dplyr::filter(grepl(pattern = "^TCGA", project_id))

# 创建目录
if (!dir.exists(clinical_dir)) { dir.create(clinical_dir) }

# 遍历每一个TCGA项目
for (i in 1:nrow(projects)) {
  project_id <- projects$project_id[i]
  
  # 查询临床数据
  clinical <- GDCquery_clinic(project = project_id, type = "clinical")  # fetch clinical indexed data (same as showed in the data portal).
  
  # 保存临床数据为 CSV 文件
  write.csv(clinical, file = file.path(clinical_dir, paste0(project_id, "_indexed_clinical.csv")), row.names = FALSE)
}

# 整理临床信息 ----------------------------------------------------------------------------
##### 步骤 1：读取已下载的临床数据 --------------------------------
clinical_dir <- args$clinical_dir  # 临床数据目录

# 获取所有临床数据文件
clinical_files <- list.files(clinical_dir, pattern = "*indexed_clinical.csv$", full.names = TRUE)
clinical_data0 <- lapply(clinical_files, read.csv)
lapply(clinical_data0, colnames)  # 检查列名

# 提取项目ID作为名字
names(clinical_data0) <- gsub("_indexed_clinical\\.csv$", "", basename(clinical_files))

# 定义所需的所有列
all_columns <- c('project','bcr_patient_barcode', 'vital_status', 'days_to_death', 'days_to_last_follow_up', 'days_to_birth', 'age_at_index', 'gender',
                 'ajcc_pathologic_stage','ajcc_pathologic_t','ajcc_pathologic_n','ajcc_pathologic_m',
                 'ajcc_clinical_t','ajcc_clinical_n','ajcc_clinical_m','updated_datetime')
table(duplicated(all_columns))

# 动态填充缺失列
clinical_data <- lapply(names(clinical_data0), function(project) {
  df <- clinical_data0[[project]]
  
  # 检查缺失列
  missing_cols <- setdiff(all_columns, colnames(df))
  
  # 如果有缺失列，则动态添加并赋值为 NA
  if (length(missing_cols) > 0) {
    df[missing_cols] <- NA  # 添加所有缺失的列并赋值为 NA
  }
  
  # 保留指定的列并确保列的顺序一致
  df <- df[, all_columns, drop = FALSE]
  return(df)
})
# 保留名字
names(clinical_data) <- names(clinical_data0)

##### 步骤 2：处理项目 --------------------------------
# 合并所有项目数据
meta <- do.call(rbind, clinical_data)
head(meta)

tmp <- meta[!is.na(meta$vital_status),]
table(tmp$vital_status)
rm(tmp)

##### 步骤 3：整合所有临床数据并处理列名---------------
# 提取癌症类型
table(meta$project)
meta$project <- stringr::str_split(meta$project, '-', simplify = TRUE)[, 2]

# 删除重复的患者记录
meta2 <- meta[!duplicated(meta$bcr_patient_barcode), ]

# 排序
table(is.na(meta2$bcr_patient_barcode))
meta2 <- meta2[order(meta2$project, meta2$bcr_patient_barcode), ]

# 设置患者ID为行名
rownames(meta2) <- meta2$bcr_patient_barcode

# 处理缺失值
meta2[meta2 == ""] <- NA

# 检查关键列数据
table(meta2$updated_datetime)
table(meta2$vital_status, useNA = "always")  # 检查生存状态

# Alive         Dead Not Reported         <NA>
#   7621         3643           16          148
Not_Reported <- meta2[meta2$vital_status == "Not Reported",]

table(is.na(meta2$bcr_patient_barcode), useNA = "always")
table(is.na(meta2$gender), useNA = "always")
table(is.na(meta2$age_at_index), useNA = "always")
table(is.na(meta2$days_to_birth), useNA = "always")

# 修改列名
colnames(meta2)
meta2 <- meta2 %>% dplyr::select(
  Cancer_Type = project,
  ID = bcr_patient_barcode,
  event = vital_status,
  death = days_to_death,
  last_followup = days_to_last_follow_up,
  age = age_at_index,
  gender = gender,
  Pathologic_Stage = ajcc_pathologic_stage,
  Pathologic_T = ajcc_pathologic_t,
  Pathologic_N = ajcc_pathologic_n,
  Pathologic_M = ajcc_pathologic_m,
  Clinical_T = ajcc_clinical_t,
  Clinical_N = ajcc_clinical_n,
  Clinical_M = ajcc_clinical_m,
  updated_datetime = updated_datetime
)
meta2 <- dplyr::select(meta2, -days_to_birth)
colnames(meta2)

#### 去除 "event", "death",  "last_followup", "age", .... "updated_datetime"  都为NA的样本
meta2 <- meta2 %>% filter(!if_all(c("event", "death", "last_followup", "age", "gender",
                                    "Pathologic_Stage", "Pathologic_T", "Pathologic_N", "Pathologic_M",
                                    "Clinical_T", "Clinical_N", "Clinical_M", "updated_datetime"), is.na))
table(meta2$Cancer_Type)

# ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP LAML  LGG LIHC LUAD LUSC MESO   OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC  UCS  UVM
#  92  412 1097  307   48  459   48  185  599  528  113  537  291  200  515  377  522  504   87  587  185  179  500  171  261  470  443  247  507  124  548   57   80

# 保存整理后的临床数据
save(meta2, file = file.path(args$clinical_dir, "TCGA_clinical_ALL_data.Rdata"))
write.csv(meta2, file.path(args$clinical_dir, "TCGA_clinical_ALL_data.csv"), row.names = FALSE)

##### 2.1 计算生存期，整理分期信息--------------------------------
load(file.path(args$clinical_dir, "TCGA_clinical_ALL_data.Rdata"))
table(meta2$Cancer_Type)

# calculate survival time
table(is.na(meta2$death), useNA = "always")
table(is.na(meta2$last_followup), useNA = "always")

# 定义time列
meta2$time <- ifelse(meta2$event == "Alive", meta2$last_followup, meta2$death)
table(is.na(meta2$time), useNA = "always")

# 将time列转换为月份数，保留两位小数
meta2$time <- round(as.numeric(meta2$time) / 30, 2)

# 定义event列，Alive=0，Dead=1
meta2$event <- case_when(
  meta2$event == 'Alive' ~ 0,
  meta2$event == 'Dead' ~ 1,
  meta2$event == 'Not Reported' ~ NA_real_,  # 或者使用 NA
  TRUE ~ NA_real_  # 默认情况下将其他所有值转换为 NA
)
table(meta2$event, useNA = "always")

# stage
# 手动标注分期
meta2 <- meta2 %>% mutate(
  stage = recode(
    Pathologic_Stage,
    "Stage 0" = "0",
    "Stage I" = "I",
    "Stage IA" = "I",
    "Stage IB" = "I",
    "Stage IS" = "I",
    "Stage II" = "II",
    "Stage IIA" = "II",
    "Stage IIB" = "II",
    "Stage IIC" = "II",
    "Stage III" = "III",
    "Stage IIIA" = "III",
    "Stage IIIB" = "III",
    "Stage IIIC" = "III",
    "Stage IV" = "IV",
    "Stage IVA" = "IV",
    "Stage IVB" = "IV",
    "Stage IVC" = "IV",
    "Stage X" = "X",  # 未知分期
    .default = NA_character_  # 如果没有匹配到，标记为 NA
  )
)

# 检查结果
table(meta2$stage, useNA = "always")

# 0    I   II  III   IV    X <NA>
# 7 2197 2239 1761  843   13 4220

# 重新标记 Pathologic_T
meta2 <- meta2 %>% mutate(
  Pathologic_T_recode = recode(
    Pathologic_T,
    "T0" = "T0",
    "T1" = "T1", "T1a" = "T1", "T1a1" = "T1", "T1b" = "T1", "T1b1" = "T1", "T1b2" = "T1", "T1c" = "T1",
    "T2" = "T2", "T2a" = "T2", "T2a1" = "T2", "T2a2" = "T2", "T2b" = "T2", "T2c" = "T2",
    "T3" = "T3", "T3a" = "T3", "T3b" = "T3",
    "T4" = "T4", "T4a" = "T4", "T4b" = "T4", "T4c" = "T4", "T4d" = "T4", "T4e" = "T4",
    "Tis" = "Tis", "TX" = "TX",
    .default = NA_character_
  )
)

# 重新标记 Pathologic_N
meta2 <- meta2 %>% mutate(
  Pathologic_N_recode = recode(
    Pathologic_N,
    "N0" = "N0", "N0 (i-)" = "N0", "N0 (i+)" = "N0", "N0 (mol+)" = "N0",
    "N1" = "N1", "N1a" = "N1", "N1b" = "N1", "N1c" = "N1", "N1mi" = "N1",
    "N2" = "N2", "N2a" = "N2", "N2b" = "N2", "N2c" = "N2",
    "N3" = "N3", "N3a" = "N3", "N3b" = "N3", "N3c" = "N3",
    "NX" = "NX",
    .default = NA_character_
  )
)

# 重新标记 Pathologic_M
meta2 <- meta2 %>% mutate(
  Pathologic_M_recode = recode(
    Pathologic_M,
    "M0" = "M0", "cM0 (i+)" = "M0",
    "M1" = "M1", "M1a" = "M1", "M1b" = "M1", "M1c" = "M1",
    "MX" = "MX",
    .default = NA_character_
  )
)

# 重新标记 Clinical_T
meta2 <- meta2 %>% mutate(
  Clinical_T_recode = recode(
    Clinical_T,
    "T1" = "T1", "T1a" = "T1", "T1b" = "T1", "T1c" = "T1",
    "T2" = "T2", "T2a" = "T2", "T2b" = "T2", "T2c" = "T2",
    "T3" = "T3", "T3a" = "T3", "T3b" = "T3",
    "T4" = "T4", "T4a" = "T4", "T4b" = "T4", "T4c" = "T4", "T4d" = "T4", "T4e" = "T4",
    "TX" = "TX",
    .default = NA_character_
  )
)

# 重新标记 Clinical_N
meta2 <- meta2 %>% mutate(
  Clinical_N_recode = recode(
    Clinical_N,
    "N0" = "N0",
    "N1" = "N1",
    "N2" = "N2", "N2a" = "N2", "N2b" = "N2", "N2c" = "N2",
    "N3" = "N3",
    "NX" = "NX",
    .default = NA_character_
  )
)

# 重新标记 Clinical_M
meta2 <- meta2 %>% mutate(
  Clinical_M_recode = recode(
    Clinical_M,
    "M0" = "M0",
    "M1" = "M1", "M1a" = "M1", "M1b" = "M1", "M1c" = "M1",
    "MX" = "MX",
    .default = NA_character_
  )
)

colnames(meta2)
table(meta2$Cancer_Type)

# ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP LAML  LGG LIHC LUAD LUSC MESO   OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC  UCS  UVM
#  92  412 1097  307   48  459   48  185  599  528  113  537  291  200  515  377  522  504   87  587  185  179  500  171  261  470  443  247  507  124  548   57   80

# 保存整理后的临床数据
save(meta2, file = file.path(args$clinical_dir, "TCGA_clinical_ALL_data.Rdata"))
write.csv(meta2, file.path(args$clinical_dir, "TCGA_clinical_ALL_data.csv"), row.names = FALSE)
