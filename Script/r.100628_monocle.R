rm(list = ls()); gc()
ORIGINAL_DIR <- "/data/nas1/huangyuqiong_OD/program156/"
output <- file.path(ORIGINAL_DIR, "10_monocle")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(output)

library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(clustree)
library(monocle)
library(ggsci)
library(readr)

future::plan("multisession", workers = 4) 
options(future.globals.maxSize = 10 * 1024^3) 
scobj <- readRDS("/data/nas1/huangyuqiong_OD/program156/07_annotation/subcluster_annotation.rds")
hubgene <- read_csv("/data/nas1/huangyuqiong_OD/program156/02_hub/hubgene.csv")
genes<-hubgene$x

scobj1 <-subset(scobj,idents =c("stroma"))
scobj1<- CreateSeuratObject(counts = GetAssayData(scobj1, assay = "RNA", layer = "counts"),
                            meta.data  = scobj1@meta.data)
scobj1 <- NormalizeData(scobj1)
scobj1 <- FindVariableFeatures(scobj1, selection.method = "vst")
scobj1 <- ScaleData(scobj1, features = rownames(scobj1))
scobj1 <- RunPCA(scobj1, features = VariableFeatures(object = scobj1),reduction.name = "pca")
ElbowPlot(scobj1)
scobj1<- RunUMAP(scobj1,dims = 1:30,reduction.name = "umap")
scobj1 <- FindNeighbors(scobj1, dims = 1:30, verbose = FALSE)
scobj1 <- FindClusters(scobj1, resolution = seq(0,0.5,0.05))
scobj1[["RNA"]]@layers$scale.data <- NULL
scobj1@meta.data$seurat_clusters <- scobj1@meta.data$RNA_snn_res.0.05
Idents(scobj1)<-scobj1@meta.data$seurat_clusters
scobj1@meta.data$celltype<-Idents(scobj1)
DimPlot(scobj1, reduction = "umap", label = T)

markers<-FindAllMarkers(scobj1)
top_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice(1:50) %>%
  ungroup() 
write.csv(top_markers,file=file.path(output,"top_marker.csv"))
saveRDS(scobj1 ,file = file.path(output,"mono_scobj.rds"))
VlnPlot(scobj1, features= c( "MIF","HSD17B2",
                            "ICAM4"),pt.size = 0,ncol=3)
ct <- scobj1@assays$RNA$counts
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  gene_id = row.names(ct),
  row.names = row.names(ct))
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
pd <- new("AnnotatedDataFrame",
          data=scobj1@meta.data)
cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)

cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds,cores=4)
cds <- detectGenes(cds, min_expr = 0.5)  
cds <- cds[fData(cds)$num_cells_expressed > 5, ]  
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(
  disp_table, 
  mean_expression >= 0.05 & dispersion_empirical > 1
)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.05)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
cds <- reduceDimension(
  cds, 
  reduction_method = "DDRTree", 
  max_components = 2,
  verbose = TRUE
)
cds <- orderCells(cds)
saveRDS(cds, file = file.path(output,"mono_monocle.rds") )

p1 = plot_cell_trajectory(cds, color_by = 'Pseudotime') 
p2 = plot_cell_trajectory(cds, color_by = 'celltype')  + scale_color_npg()
p3=plot_cell_trajectory(cds, color_by = 'State')
p<-p1+p2+p3
theme(text = element_text(family = "ArialMT"))
ggsave(filename = file.path(output,"30_cell_trajectory.pdf"),p,width=14,height=7)
ggsave(filename = file.path(output,"30_cell_trajectory.png"),p,width=14,height=7)


cds_HC<- cds[, pData(cds)$group == "HC"]
cds_HC@reducedDimS <- cds@reducedDimS[, colnames(cds_HC)]
p3<-plot_cell_trajectory(cds_HC, color_by = "group")+scale_color_manual(values = "#40E0D0")

cds_UC <- cds[, pData(cds)$group == "UC"]
cds_UC@reducedDimS <- cds@reducedDimS[, colnames(cds_UC)]
p4<-plot_cell_trajectory(cds_UC, color_by = "group")
p<-p3+p4
theme(text = element_text(family = "ArialMT"))
ggsave(filename = file.path(output,"31_group_trajectory.pdf"),p,width=14,height=7)
ggsave(filename = file.path(output,"31_group_trajectory.png"),p,width=14,height=7)

p3 <- plot_genes_in_pseudotime(cds[genes,], color_by = "Pseudotime")+
  theme(text = element_text(family = "ArialMT"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = file.path(output,"32_Genes_Jitterplot.pdf"), p3,width=7,height=7)
ggsave(filename = file.path(output,"32_Genes_Jitterplot.png"), p3,width=7,height = 7)



