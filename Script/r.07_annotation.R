rm(list = ls()); gc()
ORIGINAL_DIR <- "/data/nas1/huangyuqiong_OD/program156"
output <- file.path(ORIGINAL_DIR, "07_annotation")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(output)
library(stringi)
library(tidyverse)
library(ggplot2)
library(Seurat)
options(stringsAsFactors = F)
library(data.table)
library(tibble)
library(dplyr)
library(future)
plan("multisession", workers = 10) 
options(future.globals.maxSize = 20 * 1024^3) 

scobj <- readRDS("/data/nas1/huangyuqiong_OD/program156/06_QC/suerat_res.0.2.rds")
scobj <- JoinLayers(scobj)
allmarkers<-FindAllMarkers(scobj,min.pct=0.1,logfc.threshold=0.25,only.pos=TRUE)
saveRDS(allmarkers,file = "./all_markers.rds")
write.csv(allmarkers,"./all_marker.csv")
top_markers <- allmarkers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice(1:50) %>%
  ungroup() 
write.csv(top_markers,"./alltop50_markers.csv")

VlnPlot(scobj, features= c( "EPCAM","PHGR1","FABP1","KRT8",
                           "IGFBP7","COL3A1","PECAM1" ,"TPSAB1","TPSB2"
                           
),pt.size = 0,ncol=5)

VlnPlot(scobj, features= c( "CD3D","CD3E", 
                            "CD79A","MS4A1","MZB1"
                            
),pt.size = 0,ncol=5)
VlnPlot(scobj, features= c("HLA-DRA","HLA-DPB1",
                           "ITGAX","LYZ","CD14",
                            "CD68","CD163"
),pt.size = 0,ncol=5)
DimPlot(scobj,reduction = "umap", label = T)



scobj1<- RenameIdents(scobj,
                      "0"="T cells",
                      "2"="T cells",
                      "1"="epithelium",
                      "7"="epithelium",
                      "13"="epithelium",
                      "10"="epithelium",
                      "12"="B cells",
                      "5"="B cells",
                      "3"="plasma cells",
                      "6"="myeloid cells",
                      "9"="myeloid cells",
                      "4"="stroma",
                      "8"="mast cells",
                      "11"="endothelial cells"
)
DimPlot(scobj1, reduction = "umap", label = T)

p<-DimPlot(scobj1, reduction = "umap", label = T)
ggsave(file.path(output,"23_umap_annotation.pdf"), p, w=7,h=7)
ggsave(file.path(output,"23_umap_annotation.png"), p, w=7,h=7)


scobj1@meta.data$celltype = Idents(scobj1)
Idents(scobj1) <- factor(Idents(scobj1),levels = c("epithelium",
                                                   "stroma",
                                                   "B cells",
                                                   "plasma cells",
                                                   "T cells",
                                                   "myeloid cells",
                                                   "mast cells",
                                                   "endothelial cells"))
markers.to.plot <- c("EPCAM","PHGR1","FABP1","KRT8",
                     "IGFBP7","COL3A1",
                     "CD79A","MS4A1",
                     "MZB1",
                     "CD3D","CD3E", 
                     "HLA-DRA","HLA-DPB1",
                     "ITGAX","LYZ","CD14",
                     "CD68","CD163",
                     "TPSAB1","TPSB2",
                     "PECAM1")

p<-DotPlot(scobj1, features = markers.to.plot, dot.scale = 5) + 
  coord_flip()+
  RotatedAxis()+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))
ggsave(file.path(output,"22_dotplot_marker.pdf"), p,w=7,h=7)

ggsave(file.path(output,"22_dotplot_marker.png"), p,w=7,h=7)

saveRDS(scobj1,file = "./subcluster_annotation.rds")

