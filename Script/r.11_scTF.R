rm(list = ls()); gc()
ORIGINAL_DIR <- "11_scTF"
output <- file.path(ORIGINAL_DIR, "11_output")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(ORIGINAL_DIR)
options(future.globals.maxSize = 10 * 1024^3) 
future::plan("multisession", workers = 40) 

# 加载包 ---------------------------------------------------------------------
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)

scobj <- readRDS("10_scMetabolism/scobj_L_S_V_sub_macth.rds")
cell_types <- scobj$celltype
cell_types<- setdiff(cell_types, c("Myoid cells_1","Vascular smooth muscle cells_1","Sertoli cells_1"))
scobj1 <-subset(scobj,idents =cell_types)
scobj1<- CreateSeuratObject(counts = GetAssayData(scobj1, assay = "RNA", layer = "counts"),
                            meta.data  = scobj1@meta.data)
scobj1 <- NormalizeData(scobj1,normalization.method = "LogNormalize")
scobj1 <- FindVariableFeatures(scobj1, selection.method = "vst")
scobj1 <- ScaleData(scobj1, features = rownames(scobj))
scobj1 <- RunPCA(scobj1, features = VariableFeatures(object = scobj),reduction.name = "pca")
ElbowPlot(scobj1)
scobj1<- RunUMAP(scobj1,dims = 1:20,reduction.name = "umap")
scobj1 <- FindNeighbors(scobj1, dims = 1:20, verbose = FALSE)
scobj1 <- FindClusters(scobj1, resolution = seq(0.02, 0.4, 0.05))
scobj1[["RNA"]]@layers$scale.data <- NULL
saveRDS(scobj1,"./scobj_sub.rds")

scobj <- readRDS("11_scTF/scobj_sub.rds")
Idents(scobj) <- scobj$celltype
DimPlot(scobj,label = T)
count_matrix <- GetAssayData(scobj, assay = "RNA", layer = "counts")
write.csv(t(as.matrix(count_matrix)), file = "count__matrix.csv")


loom <- open_loom('out_SCENIC.loom')
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4]
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])
embeddings <- get_embeddings(loom)
close_loom(loom)
rownames(regulonAUC)

sub_regulonAUC <- regulonAUC[,match(colnames(scobj),colnames(regulonAUC))]
dim(sub_regulonAUC)
identical(colnames(sub_regulonAUC), colnames(scobj))
cellClusters <- data.frame(row.names = colnames(scobj),
                           seurat_clusters = as.character(scobj$seurat_clusters))

if ("celltype" %in% colnames(scobj@meta.data)) {
  if (is.factor(scobj$celltype)) {
    scobj$celltype <- droplevels(scobj$celltype)
    
    cat("更新后的分类分布（已删除空类型）:\n")
    print(table(scobj$celltype))
  } else {
    scobj$celltype <- droplevels(factor(scobj$celltype))
    warning("celltype 列已从字符型转为因子并清理空类型")
  }
} else {
  stop("scobj 中未找到 celltype 列")
}

scobj$celltype.group <- paste(scobj$celltype, scobj$group, sep = "_")
cellTypes <- data.frame(row.names = colnames(scobj),
                        celltype = scobj$celltype.group)
sub_regulonAUC[1:4,1:4]
save(sub_regulonAUC,cellTypes,cellClusters,scobj,
     file = 'for_rss_and_visual.Rdata')
selectedResolution <- "celltype" 
cellsPerGroup <- split(rownames(cellTypes),
                       cellTypes[,selectedResolution])
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]
dim(sub_regulonAUC)
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function
                                  (cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))


regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = 
                                            T, scale=T
))

dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled <- regulonActivity_byGroup_Scaled[complete.cases(regulonActivity_byGroup_Scaled), ]
regulonActivity_byGroup_Scaled <- regulonActivity_byGroup_Scaled[!apply(is.infinite(regulonActivity_byGroup_Scaled), 1, any), ]
regulon_avg_activity <- rowMeans(abs(regulonActivity_byGroup_Scaled))
top_regulons <- names(sort(regulon_avg_activity, decreasing = TRUE)[1:100])

pdf(file =file.path(output,'heatmap.pdf'),w=7,h=20)
pheatmap(regulonActivity_byGroup_Scaled[top_regulons,])
dev.off()
png(file =file.path(output,'heatmap.png'),w=7,h=20)
pheatmap(regulonActivity_byGroup_Scaled[top_regulons,])
dev.off()

regulon_avg_activity <- rowMeans(abs(regulonActivity_byGroup_Scaled))
top_regulons <- names(sort(regulon_avg_activity, decreasing = TRUE)[1:20])

if (exists("hubAUC") && nrow(hubAUC) > 0) {
  top_regulons <- rownames(hubAUC)
}

activity_matrix <- regulonActivity_byGroup_Scaled[top_regulons, ]

cor_matrix <- cor(t(activity_matrix), method = "pearson") 

adj_threshold <- 0.6 
adj_matrix <- ifelse(abs(cor_matrix) > adj_threshold & cor_matrix != 1, 1, 0)

network <- graph_from_adjacency_matrix(
  adj_matrix,
  mode = "undirected",
  weighted = NULL,
  diag = FALSE
)

V(network)$size <- scale(regulon_avg_activity[top_regulons]) * 5 + 8  
V(network)$label.cex <- 2 
V(network)$color <- "#4DB"  

set.seed(42) 
p <- ggraph(network, layout = "kk") + 
  geom_edge_link(alpha = 0.1, width = 0.3) +
  geom_node_point(aes(size = size), color = V(network)$color) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  labs(title = "Core TF Interaction Network")

ggsave(file.path(output, "TF_network.png"), plot = p, width =7, height = 7)
ggsave(file.path(output, "TF_network.pdf"), plot = p, width = 7, height = 7)

centrality <- data.frame(
  TF = V(network)$name,
  Degree = degree(network),
  Betweenness = betweenness(network),
  Eigenvector = eigen_centrality(network)$vector
)
write.csv(centrality, file.path(output, "network_centrality.csv"))
