volcano_plot_custom_genes <- function(df, symbol_label, logFC_col="logFC", FDR_col="adj.P.Val", 
                                Symbol_col="Symbol", y_increased=20,
                                colours=c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142")) {
  library(ggplot2)
  library(ggrepel)
  library(ggthemes)

  ## Check that belongs to the column name is correct
  if (!logFC_col %in% colnames(df)) {
    cat(paste0("The '", logFC_col, "' column is not found in the dataset.\n"))
    logFC_col <- readline("Please enter the name of the 'logFC' column: ")
  }
  
  if (!FDR_col %in% colnames(df)) {
    cat(paste0("The '", FDR_col, "' column is not found in the dataset.\n"))
    FDR_col <- readline("Please enter the name of the 'FDR' column: ")
  }
  
  if (!Symbol_col %in% colnames(df)) {
    cat(paste0("The '", Symbol_col, "' column is not found in the dataset.\n"))
    Symbol_col <- readline("Please enter the name of the 'Symbol' column: ")
  }
  
  print(paste0("Column name of logFC_col is: ", logFC_col))
  print(paste0("Column name of FDR_col is: ", FDR_col))
  print(paste0("Column name of Symbol_col is: ", Symbol_col))
  
  y_upper <- ceiling(max(-log10(df[[sym(FDR_col)]])))
  #logFC_min <- floor(min(df[[sym(logFC_col)]]))
  #logFC_max <- ceiling(max(df[[sym(logFC_col)]]))
  x_abs <- max(abs(c(floor(min(df[[sym(logFC_col)]])), ceiling(max(df[[sym(logFC_col)]])))))
  
  ## When not entered "symbol_label"
  if (missing(symbol_label) || is.null(symbol_label)) {
    symbol_labels <- readline("Please enter the symbol labels: ")
    ## Convert the user input to a character vector
    symbol_label <- eval(parse(text=symbol_labels))
    ## If user enters a single symbol label as string, convert it to a vector
    if (length(symbol_label) == 1) {
      symbol_label <- c(symbol_label)
    }
  }
  
  plot <- ggplot(data = df) + 
    geom_point(aes(x = !!sym(logFC_col), y = -log10(!!sym(FDR_col)), color = !!sym(logFC_col), size = -log10(!!sym(FDR_col)))) +  
    geom_text_repel(data = df %>% filter(Symbol %in% symbol_label),
                    aes(x = !!sym(logFC_col), y = -log10(!!sym(FDR_col)), label = !!sym(Symbol_col)),  
                    nudge_x = 0.5, nudge_y = 0.2, segment.curvature = -0.1, segment.ncp = 3,
                    direction = "y", hjust = "left",
                    max.overlaps = 200)+
    scale_color_gradientn(colours = colours, values = seq(0, 1, 0.2)) +
    scale_fill_gradientn(colours = colours, values = seq(0, 1, 0.2)) +
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 4) + 
    scale_size(range = c(1,6)) + 
    ggtitle(label = "Volcano Plot") + 
    xlim(-x_abs, x_abs) + 
    ylim(c(-1, y_upper+y_increased)) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.background = element_roundrect(color = "#808080", linetype = 1),
          axis.text = element_text(size = 13, color = "#000000"),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5)) + 
    coord_cartesian(clip = "off") + 
    annotation_custom(grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#74add1")), 
      xmin = -x_abs, xmax = -1, ymin = y_upper+(y_increased+1), ymax = y_upper+(y_increased+1)) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Down",
        gp = grid::gpar(col = "#74add1")),
      xmin = -x_abs, xmax = -1, ymin = y_upper+(y_increased+1), ymax = y_upper+(y_increased+1)) +
    annotation_custom(
      grob = grid::segmentsGrob(
        y0 = unit(-10, "pt"),
        y1 = unit(-10, "pt"),
        arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
        gp = grid::gpar(lwd = 3, col = "#d73027")), 
      xmin = x_abs, xmax = 1, ymin = y_upper+(y_increased+1), ymax = y_upper+(y_increased+1)) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Up",
        gp = grid::gpar(col = "#d73027")),
      xmin = x_abs, xmax = 1, ymin = y_upper+(y_increased+1), ymax = y_upper+(y_increased+1))
  
  return(plot)
}
