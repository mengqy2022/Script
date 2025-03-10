#!/usr/bin/env Rscript

library(argparse)
library(Seurat)            # 用于单细胞RNA-seq分析
library(dplyr)             # 数据操作
library(ggplot2)           # 可视化
library(clusterProfiler)   # 富集分析
library(org.Hs.eg.db)      # 基因注释数据库
library(Harmony)           # 用于批次效应去除
library(DoubletFinder)     # 双细胞预测

# 创建命令行解析器
parser <- ArgumentParser(description = 'Single-cell RNA-seq analysis pipeline',
                         epilog = '\nRscript single_cell_analysis.R --sample_paths /path/to/sample1 /path/to/sample2 --output_dir /path/to/output --resolution 0.6 --dims 1:20 --use_doublet_finder TRUE')


# 添加命令行参数
parser$add_argument('--sample_paths', type = 'character', nargs = '+', required = TRUE, help = 'Paths to the sample directories')
parser$add_argument('--output_dir', type = 'character', default = 'output', help = 'Directory to save the results')
parser$add_argument('--resolution', type = 'numeric', default = 0.5, help = 'Clustering resolution (default: 0.5)')
parser$add_argument('--dims', type = 'integer', default = 1:20, help = 'PCA dimensions for clustering (default: 1:20)')
parser$add_argument('--use_doublet_finder', type = 'logical', default = TRUE, help = 'Whether to use DoubletFinder for doublet detection (default: TRUE)')

args <- parser$parse_args()

# 1. 数据加载与质控函数
load_and_qc <- function(sample_paths, min_cells = 3, min_features = 200, max_mt_percent = 5, use_doublet_finder = TRUE) {
  samples <- lapply(sample_paths, function(path) {
    data <- Read10X(data.dir = path)
    seurat_obj <- CreateSeuratObject(counts = data, project = "SingleCell", min.cells = min_cells, min.features = min_features)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    
    # 去除低质量细胞和双细胞
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < max_mt_percent)
    
    # 如果使用DoubletFinder
    if (use_doublet_finder) {
      seurat_obj <- doublet_removal(seurat_obj)
    }
    
    return(seurat_obj)
  })
  
  return(samples)
}

# 2. 基于指标（nFeature_RNA 和 nCount_RNA）去除双细胞的函数
doublet_removal <- function(seurat_obj, min_features = 200, max_features = 2500, min_count = 500, max_count = 10000) {
  # 通过检查基因数和转录本数去除双细胞
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_features & nFeature_RNA < max_features)
  seurat_obj <- subset(seurat_obj, subset = nCount_RNA > min_count & nCount_RNA < max_count)
  
  # 可选：使用 DoubletFinder 去除预测的双细胞
  seurat_obj <- DoubletFinder::doubletFinder_v3(seurat_obj, PCs = 1:20, pN = 0.25, pK = 0.05, nExp = round(0.05 * ncol(seurat_obj)))
  
  # 根据DoubletFinder的预测，过滤掉双细胞
  seurat_obj <- subset(seurat_obj, subset = DF.classifications != "Doublet")
  
  return(seurat_obj)
}

# 3. 数据标准化和高变基因选择函数
normalize_and_find_variable_genes <- function(seurat_objs) {
  for (i in 1:length(seurat_objs)) {
    seurat_objs[[i]] <- NormalizeData(seurat_objs[[i]])
    seurat_objs[[i]] <- FindVariableFeatures(seurat_objs[[i]], selection.method = "vst", nfeatures = 2000)
  }
  return(seurat_objs)
}

# 4. 多样本整合函数
integrate_samples <- function(seurat_objs) {
  # 合并样本
  merged_samples <- merge(seurat_objs[[1]], y = seurat_objs[2:length(seurat_objs)], add.cell.ids = paste0("Sample", 1:length(seurat_objs)))
  merged_samples <- NormalizeData(merged_samples)
  merged_samples <- FindVariableFeatures(merged_samples, selection.method = "vst", nfeatures = 2000)
  
  return(merged_samples)
}

# 5. 批次效应去除函数
remove_batch_effects <- function(merged_samples) {
  merged_samples <- RunPCA(merged_samples, features = VariableFeatures(object = merged_samples))
  merged_samples <- RunHarmony(merged_samples, group.by.vars = "orig.ident")  # "orig.ident" 是批次信息
  
  return(merged_samples)
}

# 6. 降维和聚类函数
perform_clustering <- function(merged_samples, dims = 1:20, resolution = 0.5) {
  merged_samples <- FindNeighbors(merged_samples, dims = dims)
  merged_samples <- FindClusters(merged_samples, resolution = resolution)
  
  # 运行UMAP降维
  merged_samples <- RunUMAP(merged_samples, dims = dims)
  
  return(merged_samples)
}

# 7. 差异表达分析函数
find_differential_expression <- function(merged_samples) {
  # 对所有聚类进行差异表达分析
  cluster_ids <- unique(merged_samples$seurat_clusters)
  all_results <- list()
  
  for (i in 1:(length(cluster_ids) - 1)) {
    for (j in (i + 1):length(cluster_ids)) {
      cluster_1 <- cluster_ids[i]
      cluster_2 <- cluster_ids[j]
      
      cluster_markers <- FindMarkers(merged_samples, ident.1 = cluster_1, ident.2 = cluster_2)
      all_results[[paste0("Cluster_", cluster_1, "_vs_Cluster_", cluster_2)]] <- cluster_markers
    }
  }
  
  return(all_results)
}

# 8. 富集分析函数
perform_go_enrichment <- function(diff_genes) {
  # 选取显著差异基因（p值<0.05）
  gene_list <- rownames(diff_genes)[which(diff_genes$p_val_adj < 0.05)]
  
  # 转换基因ID为Entrez ID
  gene_symbols <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # 进行GO富集分析
  go_results <- enrichGO(gene = gene_symbols$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  return(go_results)
}

# 主分析函数：执行所有分析步骤
single_cell_analysis <- function(sample_paths, resolution = 0.5, dims = 1:20, output_dir = "output", use_doublet_finder = TRUE) {
  # 1. 数据加载与质控
  samples <- load_and_qc(sample_paths, use_doublet_finder = use_doublet_finder)
  
  # 2. 数据标准化与高变基因选择
  samples <- normalize_and_find_variable_genes(samples)
  
  # 3. 多样本整合
  merged_samples <- integrate_samples(samples)
  
  # 4. 批次效应去除
  merged_samples <- remove_batch_effects(merged_samples)
  
  # 5. 降维与聚类
  merged_samples <- perform_clustering(merged_samples, dims = dims, resolution = resolution)
  
  # 6. 差异表达分析
  cluster_markers_list <- find_differential_expression(merged_samples)
  
  # 7. 富集分析（对每个差异基因集进行富集分析）
  go_results_list <- list()
  for (key in names(cluster_markers_list)) {
    cluster_markers <- cluster_markers_list[[key]]
    go_results <- perform_go_enrichment(cluster_markers)
    go_results_list[[key]] <- go_results
  }
  
  # 保存结果
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # 保存 Seurat 对象、差异分析结果和富集分析结果
  saveRDS(merged_samples, file.path(output_dir, "merged_samples.rds"))
  saveRDS(cluster_markers_list, file.path(output_dir, "cluster_markers_list.rds"))
  saveRDS(go_results_list, file.path(output_dir, "go_results_list.rds"))
  
  # 可视化结果
  DimPlot(merged_samples, reduction = "umap", group.by = "seurat_clusters") + 
    ggsave(file.path(output_dir, "umap_plot.png"))
  
  # GO 富集分析可视化
  for (key in names(go_results_list)) {
    dotplot(go_results_list[[key]]) + 
      ggsave(file.path(output_dir, paste0(key, "_go_enrichment_plot.png")))
  }
}

# 执行主分析函数
single_cell_analysis(sample_paths = args$sample_paths, 
                     resolution = args$resolution, 
                     dims = args$dims, 
                     output_dir = args$output_dir,
                     use_doublet_finder = args$use_doublet_finder)