rm(list = ls()); gc()
ORIGINAL_DIR <- "/data/nas1/huangyuqiong_OD/program156/"
output <- file.path(ORIGINAL_DIR, "09_cellchat")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(output)

library(Seurat)
library(ggsignif)
library(ggpubr)
library(cowplot)
library(dplyr)
library(NMF)
library(ggalluvial)
library(CellChat)
library(pheatmap)
library(qs)
library(patchwork)
library(future)
options(stringsAsFactors = FALSE)
library(gridExtra)
options(future.plan = NULL)
future::plan("multisession", workers = 10) 
options(future.globals.maxSize = 10 * 1024^3) 
scobj <- readRDS("/data/nas1/huangyuqiong_OD/program156/07_annotation/subcluster_annotation.rds")

HC.input <- subset(scobj,group == "HC")
HC.input <-HC.input[["RNA"]]$data
identityHC <- subset(scobj@meta.data,group == "HC")
identityHC <- subset(identityHC, select = "celltype")
cellchatHC <- createCellChat(object = HC.input, meta = identityHC,  group.by = "celltype")

CellChatDB <- CellChatDB.human
unique(CellChatDB$interaction$annotation)
cellchatHC@DB <- CellChatDB
cellchatHC <- subsetData(cellchatHC)

cellchatHC <- identifyOverExpressedGenes(cellchatHC)
cellchatHC<- identifyOverExpressedInteractions(cellchatHC)
cellchatHC<- projectData(cellchatHC, PPI.human)
cellchatHC@idents <- droplevels(cellchatHC@idents)
cellchatHC<- computeCommunProb(cellchatHC, raw.use = TRUE)

cellchatHC<- filterCommunication(cellchatHC, min.cells = 3)

cellchatHC <- computeCommunProbPathway(cellchatHC)

cellchatHC <- aggregateNet(cellchatHC)

saveRDS(cellchatHC,file.path(output,"cellchatHC.rds"))

UC.input <- subset(scobj,group == "UC")
UC.input <-UC.input[["RNA"]]$data
identityUC <- subset(scobj@meta.data,group == "UC")
identityUC <- subset(identityUC, select = "celltype")
cellchatUC <- createCellChat(object = UC.input, meta = identityUC,  group.by = "celltype")

cellchatUC@DB <- CellChatDB
cellchatUC <- subsetData(cellchatUC)

cellchatUC <- identifyOverExpressedGenes(cellchatUC)
cellchatUC<- identifyOverExpressedInteractions(cellchatUC)
cellchatUC<- projectData(cellchatUC, PPI.human)
cellchatUC@idents <- droplevels(cellchatUC@idents)
cellchatUC<- computeCommunProb(cellchatUC, raw.use = TRUE)
cellchatUC<- filterCommunication(cellchatUC, min.cells = 3)
cellchatUC<- computeCommunProbPathway(cellchatUC)
cellchatUC <- aggregateNet(cellchatUC)

saveRDS(cellchatUC,file.path(output,"cellchatUC.rds"))

object.list <- list(UC = cellchatUC, HC = cellchatHC)       
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
saveRDS(cellchat,file.path(output,"cellchat.rds"))




pdf(file.path(output,"26_count_celltype_group.pdf"),width=14,height = 14)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()
png(file.path(output,"26_count_celltype_group.png"),width = 2 * 480, height = 960)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

pdf(file.path(output,"27_weight_celltype_group.pdf"),width=14,height = 14)
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("weight of interactions - ", names(object.list)[i]))
}
dev.off()
png(file.path(output,"27_weight_celltype_group.png"),width = 2 * 480, height = 960)
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("weight of interactions - ", names(object.list)[i]))
}
dev.off()


pos.dataset = "UC"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, 
                                       thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "UC",
                              ligand.logFC = 0.2, receptor.logFC =0.2)

net.down <- subsetCommunication(cellchat, net = net, datasets = "UC",
                                ligand.logFC = -0.2, receptor.logFC = -0.2)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use ="stroma",
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", 
                                            names(object.list)[1]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        sources.use ="stroma", 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", 
                                            names(object.list)[1]))

png(file.path(output,"28_keytype_sources_DE.png"), width =1000, height = 500) 
gg1 + gg2
dev.off()
pdf(file.path(output,"28_keytype_sources_DE.pdf"), width =15, height = 7) 
gg1 + gg2
dev.off()

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        targets.use ="stroma",
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Up-regulated signaling in ", 
                                            names(object.list)[1]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, 
                        targets.use = "stroma", 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,
                        title.name = paste0("Down-regulated signaling in ", 
                                            names(object.list)[1]))

png(file.path(output,"29_keytype_targets_DE.png"), width =1000, height = 600) 
gg1 + gg2
dev.off()
pdf(file.path(output,"29_keytype_targets_DE.pdf"), width =15, height = 8) 
gg1 + gg2
dev.off()

