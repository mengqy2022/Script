# load("alldiff.RData")
# res <- res |> rownames_to_column("SYMBOL") |> arrange(desc(log2FoldChange))
# 
# ##给差异基因打标签，logFC > 1且 padj < 0.05认为是上调基因，logFC < -1且 padj < 0.05认为是下调基因
# df <- res |>  
#   mutate(significant = case_when(log2FoldChange > 1 & padj < 0.05 ~ "Up",
#                                  abs(log2FoldChange) < 1 | padj > 0.05 ~ "None",
#                                  log2FoldChange < -0.5 & padj < 0.05 ~ "Down"))
# df$significant |> as.factor()
# 
# head(df)

create_volcano_plot <- function(res_data_path, symbol_col_name = "SYMBOL", log2fc_col_name = "log2FoldChange") {
  # loading packages---------------------------
  library(ggrepel)
  library(ggfun)
  library(grid)
  library(tidyverse)
  
  # loading data-------------------------------
  res <- readRDS(res_data_path) %>%
    rownames_to_column(symbol_col_name) %>%
    arrange(desc(!! rlang::sym(log2fc_col_name)))
  
  ## 给差异基因打标签，logFC > 1且 padj < 0.05认为是上调基因，logFC < -1且 padj < 0.05认为是下调基因
  df <- res %>%
    mutate(significant = case_when(!! rlang::sym(log2fc_col_name) > 1 & padj < 0.05 ~ "Up",
                                   abs(!! rlang::sym(log2fc_col_name)) < 1 | padj > 0.05 ~ "None",
                                   !! rlang::sym(log2fc_col_name) < -1 & padj < 0.05 ~ "Down"))
  
  df$significant <- as.factor(df$significant)
  
  head(df)
  
  # plot--------------------------------
  p <- ggplot(data = df) + 
    geom_point(aes(x = !! rlang::sym(log2fc_col_name), y = -log10(padj), 
                   color = !! rlang::sym(log2fc_col_name),
                   size = -log10(padj))) + 
    geom_text_repel(data =  df %>%
                      tidyr::drop_na() %>%
                      dplyr::filter(significant != "None") %>%
                      dplyr::arrange(desc(-log10(padj))) %>%
                      dplyr::slice(1:10) %>%
                      dplyr::filter(significant == "Up"),
                    aes(x = !! rlang::sym(log2fc_col_name), y = -log10(padj), label = !! rlang::sym(symbol_col_name)),
                    nudge_x = 0.5,
                    nudge_y = 0.2,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    direction = "y",
                    hjust = "left",
                    max.overlaps = 200
    )+
    geom_text_repel(data =  df %>%
                      tidyr::drop_na() %>%
                      dplyr::filter(significant != "None") %>%
                      dplyr::filter(significant != "Up") %>%
                      dplyr::arrange(desc(-log10(padj))) %>%
                      dplyr::slice(1:10) %>%
                      dplyr::filter(significant == "Down"),
                    aes(x = !! rlang::sym(log2fc_col_name), y = -log10(padj), label = !! rlang::sym(symbol_col_name)),
                    box.padding = 0.5,
                    nudge_x = -0.2,
                    nudge_y = 0.2,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    segment.angle = 20,
                    direction = "y", 
                    hjust = "left",
                    max.overlaps = 200
    ) + 
    scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                          values = seq(0, 1, 0.2)) +
    scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                         values = seq(0, 1, 0.2)) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 4) + 
    scale_size(range = c(1,7)) + 
    ggtitle(label = "Volcano Plot") + 
    xlim(c(-8, 8)) +  # 根据你自己的数据选择合适范围!!!
    ylim(c(-1, 40)) + # 根据你自己的数据选择合适范围!!!
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.background = element_roundrect(color = "#808080", linetype = 1),
          axis.text = element_text(size = 13, color = "#000000"),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5)
    ) + 
    annotate(geom = "text", x = 4, y = 0, label = "p = 0.05", size = 5) + 
    coord_cartesian(clip = "off") + 
    annotation_custom(
      grob = grid::segmentsGrob(
        y0 = unit(-10, "pt"),
        y1 = unit(-10, "pt"),
        arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
        gp = grid::gpar(lwd = 3, col = "#74add1")
      ), 
      xmin = -8, 
      xmax = -1,
      ymin = 122,
      ymax = 122
    ) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Down",
        gp = grid::gpar(col = "#74add1")
      ),
      xmin = -8, 
      xmax = -1,
      ymin = 122,
      ymax = 122
    ) +
    annotation_custom(
      grob = grid::segmentsGrob(
        y0 = unit(-10, "pt"),
        y1 = unit(-10, "pt"),
        arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
        gp = grid::gpar(lwd = 3, col = "#d73027")
      ), 
      xmin = 8, 
      xmax = 1,
      ymin = 122,
      ymax = 122
    ) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Up",
        gp = grid::gpar(col = "#d73027")
      ),
      xmin = 8, 
      xmax = 1,
      ymin = 122,
      ymax = 122
    ) 
  
  print(p)
}

# 使用示例
# create_volcano_plot("alldiff.RDS", symbol_col_name = "SYMBOL", log2fc_col_name = "log2FoldChange")
