rm(list = ls()); gc()
ORIGINAL_DIR <- "/data/nas1/huangyuqiong_OD/program156/"
output <- file.path(ORIGINAL_DIR, "08_keytype")
if (!dir.exists(output)) {
  dir.create(output, recursive = TRUE)
}
setwd(output)

library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)
scobj <- readRDS("/data/nas1/huangyuqiong_OD/program156/07_annotation/subcluster_annotation.rds")

data <- as.data.frame(table(scobj@meta.data$sample.group,scobj@meta.data$celltype))
colnames(data) <- c("sample","CellType","Freq")
table(scobj@meta.data$celltype)
df <- data %>% 
  group_by(sample) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = Freq/Total) %>% 
  as.data.frame()

df$CellType  <- factor(df$CellType,levels =c("epithelium",
                                             "stroma",
                                             "B cells",
                                             "plasma cells",
                                             "T cells",
                                             "myeloid cells",
                                             "mast cells",
                                             "endothelial cells"))
write.csv(df,file = "./cell_percent.csv")

p<-ggplot(df, aes(x = sample, y = Percent, fill = CellType)) +
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 1, width = 0.95) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()+
  theme(axis.text.x  = element_text(angle = 45, hjust = 1) )
ggsave(file.path(output,"24_celltype_percent.pdf"), p,width=10,height=7)

ggsave(file.path(output,"24_celltype_percent.png"), p,width=10,height=7)



allmarkers2<-FindAllMarkers(scobj)
saveRDS(allmarkers2,file =file.path(output, "celltype_markers.rds"))
write.csv(allmarkers2,file.path(output, "celltype_marker.csv"))
markers <- allmarkers2 %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice(1:20) %>%
  ungroup() 
write.csv(markers,file.path(output, "top20_markers.csv"))

#hubgene####
library(ggplot2)
library(dplyr)
library(reshape2)
library(readr)
hubgene <- read_csv("/data/nas1/huangyuqiong_OD/program156/02_hub/hubgene.csv")
genes<-hubgene$x
DotPlot(scobj,features = genes)


p<-DotPlot(scobj, features = genes, dot.scale = 5) + 
  RotatedAxis()+
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8))+
  theme(text = element_text(family = "ArialMT"))
#ggsave(file.path(output,"dotplot_hub.pdf"), p,width=7,height=7)
#ggsave(file.path(output,"dotplot_hub.png"), p,width=7,height=7)

library(reshape2)
library(patchwork)
library(ggpubr)

expr_data <- FetchData(scobj, vars = c("group", "ident", genes)) 
colnames(expr_data)[1:2] <- c("group", "CellType") 

expr_data_long <- melt(expr_data, id.vars = c("group", "CellType"), variable.name = "Gene", value.name = "Expression")
group_colors <- c("UC" = "#A52A2A", "HC" = "#87CEEB") 
expr_data_long$group<-factor(expr_data_long$group,levels = c("HC","UC"))
plots_list <- list()
for(cell in unique(expr_data_long$CellType)){
  cellname <- cell
  p <- ggplot(subset(expr_data_long, CellType == cell), aes(x = group, y = Expression, fill = group)) +
    geom_violin(alpha = 0.8) + 
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) + 
    facet_wrap(~Gene, scales = "free_y", ncol = length(unique(expr_data_long$Gene))) + 
    scale_fill_manual(values = group_colors) +
    stat_compare_means(aes(group = group), method = "wilcox.test", label = "p.signif", size = 3,label.x = 1.3) + 
    theme_classic() +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 5),
      legend.position = "none"
    ) +
    labs(title = cell, x = "group", y = "Expression")+
    theme(text = element_text(family = "ArialMT"))
  
  plots_list[[cell]] <- p
}
final_plot <- wrap_plots(plots_list, ncol = 4)

ggsave(filename = file.path(output,"25_hubgene_exp_celltype.pdf"), plot =final_plot,width=14,height=7)
ggsave(filename = file.path(output,"25_hubgene_exp_celltype.png"), plot = final_plot,width=14,height=7)

