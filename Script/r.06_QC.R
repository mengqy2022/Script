rm(list = ls()); gc()
ORIGINAL_DIR <- "/data/nas1/huangyuqiong_OD/program156/"
output <- file.path(ORIGINAL_DIR, "06_QC")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(output)
library(future)
plan("multisession", workers = 4) 
options(future.globals.maxSize = 50 * 1024^3)  
library(ggplot2)
library(ggrepel)
library(harmony)
library(Seurat)
library(patchwork)
scobj <- readRDS("/data/nas1/huangyuqiong_OD/program156/00_rawdata/seurat.rds")
scobj <- JoinLayers(scobj)
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-") 
pdf(file.path(output,"18a_scobj_feature.pdf"),w=18,h=7)
p <- VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3,
             pt.size = 0)
p[[1]] <- p[[1]] & scale_y_continuous(breaks = seq(0, max(scobj@meta.data$nFeature_RNA), by = 2500))
p[[2]] <- p[[2]] & scale_y_continuous(
  labels = function(x) paste0(round(x / 10000, 1), "W"),
  breaks = seq(0, max(scobj@meta.data$nCount_RNA), by = 10000)
)
p[[3]] <- p[[3]] & scale_y_continuous(breaks = seq(0, 100, by = 10))
print(p)
dev.off()
png(file.path(output,"18a_scobj_feature.png"),w=1500,h=500)
p <- VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
             ncol = 3,
             pt.size = 0)
p[[1]] <- p[[1]] & scale_y_continuous(breaks = seq(0, max(scobj@meta.data$nFeature_RNA), by = 2500))
p[[2]] <- p[[2]] & scale_y_continuous(
  labels = function(x) paste0(round(x / 10000, 1), "W"),
  breaks = seq(0, max(scobj@meta.data$nCount_RNA), by = 10000)
)
p[[3]] <- p[[3]] & scale_y_continuous(breaks = seq(0, 100, by = 10))
print(p)
dev.off()



scobj <- subset(scobj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30 & nCount_RNA<10000)   

pdf(file.path(output,"18b_scobj_feature_L.pdf"),w=18,h=7)
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        pt.size = 0)
dev.off()
png(file.path(output,"18b_scobj_feature_L.png"),w=1500,h=500)
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, 
        pt.size = 0)
dev.off()

saveRDS(scobj,"./seurat_filtered.rds")

scobj[["RNA"]]<-split(scobj[["RNA"]], f = scobj$orig.ident)
scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)

top10_hvg <- HVFInfo(scobj) %>%
  arrange(desc(variance.standardized)) %>%
  head(10) %>%
  rownames()
pdf(file.path(output,"19_VariableFeatures.pdf"))
plot <- VariableFeaturePlot(scobj, 
                            cols = c("gray", "red"), 
                            pt.size = 1.5,        
                            assay = "RNA")  
plot + 
  geom_text_repel(data = HVFInfo(scobj)[top10_hvg, ], 
                  aes(mean, variance.standardized, label = rownames(HVFInfo(scobj)[top10_hvg, ])),
                  size = 4, 
                  box.padding = 0.3, 
                  max.overlaps = Inf) + 
  ggtitle("Top 2000 Highly Variable Genes (HVGs)") + 
  theme_classic(base_size = 12)
dev.off()
png(file.path(output,"19_VariableFeatures.png"))
plot <- VariableFeaturePlot(scobj, 
                            cols = c("gray", "red"), 
                            pt.size = 1.5,        
                            assay = "RNA")  
plot + 
  geom_text_repel(data = HVFInfo(scobj)[top10_hvg, ], 
                  aes(mean, variance.standardized, label = rownames(HVFInfo(scobj)[top10_hvg, ])),
                  size = 4, 
                  box.padding = 0.3, 
                  max.overlaps = Inf) + 
  ggtitle("Top 2000 Highly Variable Genes (HVGs)") + 
  theme_classic(base_size = 12)
dev.off()
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, 
								features = VariableFeatures(object = scobj),
								reduction.name = "pca")
pdf(file.path(output,"20b_pca_ElbowPlot.pdf"))
ElbowPlot(scobj,ndims = 30)
dev.off()
png(file.path(output,"20b_pca_ElbowPlot.png"))
ElbowPlot(scobj,ndims = 30)
dev.off()

set.seed(123)
scobj <- JackStraw(scobj,dims = 30)
scobj <- ScoreJackStraw(scobj,dim=1:30)
pdf(file.path(output,"20a_pca_JackStrawPlot.pdf"))
JackStrawPlot(scobj,dims = 1:30)
dev.off()
png(file.path(output,"20a_pca_JackStrawPlot.png"))
JackStrawPlot(scobj,dims = 1:30)
dev.off()

scobj <- IntegrateLayers(
  object = scobj,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  group.by.vars = "orig.ident",
  verbose = FALSE)

scobj<- RunUMAP(scobj, reduction = "harmony", dims = 1:30,reduction.name = "umap")

scobj[["RNA"]]@layers$scale.data <- NULL

saveRDS(scobj,"./scobj_umap.rds")

scobj <- FindNeighbors(
  scobj, 
  reduction = "harmony",   
  dims = 1:30,            
  verbose = FALSE)

scobj <- FindClusters(scobj, resolution = seq(0,1,0.05))
saveRDS(scobj,file = "./seurat_unannotation.rds")

scobj@meta.data$seurat_clusters <- scobj@meta.data$RNA_snn_res.0.2
Idents(scobj)<-scobj@meta.data$seurat_clusters
p<-DimPlot(scobj, reduction = "umap", label = T)
ggsave(file.path(output,"21_dimplot_0.2.pdf"), p)
ggsave(file.path(output,"21_dimplot_0.2.png"), p)
saveRDS(scobj,file = "./suerat_res.0.2.rds")

DimPlot(scobj, group.by="group",reduction = "umap", label = T)

