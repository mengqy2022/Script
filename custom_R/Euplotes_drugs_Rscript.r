#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Email     : 15877464851@163.com
# @Time      : 2024/10/16

options("repos" = "https://mirrors.ustc.edu.cn/CRAN/")
options(BioC_mirror = "https://mirrors.ustc.edu.cn/bioc/")

# 检查并加载包，若未安装则安装
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
      install.packages(pkg, ask = FALSE)
      library(pkg, character.only = TRUE)
    }
  }
}

base_packages <- c("argparse", "BiocManager", "R6", "readxl", "ggrepel", "ggpubr", "ggprism", "ggsci", "rstatix", "ggthemr")
load_packages(base_packages)

# 安装和加载Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("tidyverse", ask = FALSE, force = TRUE)

library(R6)
library(tidyverse)

# 定义 PlottingClass
PlottingClass <- R6Class("PlottingClass",
  public = list(
    data = NULL,
    color_palette = NULL,

    initialize = function(data, var, color_palette) {
      self$data <- self$factor_reorder(data, var)
      self$color_palette <- color_palette
    },

    factor_reorder = function(data, var) {
      if (!"Concentration" %in% colnames(data)) {
        stop("错误：数据中缺少 'Concentration' 列。")
      }
      
      if (!is.character(var) || length(var) < 1) {
        stop("错误：var 参数必须是长度大于0的字符向量。")
      }
      
      if (!all(var %in% unique(data$Concentration))) {
        stop("错误：var 参数中的某些水平在 'Concentration' 列中不存在。")
      }
      
      data$Concentration <- factor(data$Concentration, levels = var)
      return(data)
    },

    set_theme = function() {
      theme_prism() +
        theme(strip.text = element_text(size = 18),
              axis.line = element_line(color = "black", linewidth = 0.4),
              axis.text = element_text(color = "black", size = 15),
              axis.title = element_text(color = "black", size = 20),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_line(linewidth = 0.2, color = "#e5e5e5"))
    },

    bubble = function(title = "Penicillin treatment") {
      plot <- ggplot(self$data, aes(x = Days, y = Concentration)) +
        geom_point(aes(colour = Concentration, size = Totals), alpha = .7) +
        geom_text(aes(label = Totals), vjust = -1.1, size = 4, fontface = 'bold') +
        scale_size_continuous(range = c(3, 10)) +
        labs(title = title) +
        guides(colour = 'none', size = 'none') +
        self$set_theme() +
        scale_color_manual(values = self$color_palette) +
        scale_x_continuous(breaks = seq(0, 6, by = 1))
      
      if ("Times" %in% colnames(self$data)) {
        plot <- plot + facet_grid(Times ~ Repeat, scales = 'free') +
          theme(strip.text.x = element_text(size = 15),
                strip.text.y = element_text(size = 15))
      } else {
        plot <- plot + facet_grid(Repeat, scales = 'free')
      }
      
      return(plot)
    },

    plot_smooth = function(days_col = "Days", concentration_col = "Concentration", totals_col = "Totals", times_col = "Times") {
      Record_1 <- self$data %>%
        group_by(!!sym(concentration_col), !!sym(days_col)) %>%
        summarise(mean_totals = as.integer(mean(!!sym(totals_col), na.rm = TRUE)), .groups = 'drop')

      best_in_class <- Record_1 %>%
        group_by(!!sym(concentration_col)) %>%
        slice_max(mean_totals, n = 1, with_ties = FALSE) %>%
        distinct(mean_totals, .keep_all = TRUE)

      start <- Record_1 %>%
        filter(!!sym(days_col) == 0)
  
      n <- length(unique(self$data$Repeat)) * length(unique(self$data$Times))

      plot <- ggplot(Record_1, aes(x = !!sym(days_col), y = mean_totals, color = !!sym(concentration_col))) +
        geom_smooth(aes(linetype = !!sym(concentration_col)), se = TRUE, alpha = .05) +
        geom_point(size = 3, alpha = .8) +
        self$set_theme() +
        scale_x_continuous(breaks = seq(0, 6, by = 1)) +
        geom_label_repel(data = best_in_class, aes(label = mean_totals), show.legend = FALSE) +
        geom_label_repel(data = start, aes(label = mean_totals), show.legend = FALSE) +
        xlab("Days") +
        ylab("Cell numbers") +
        annotate("text", x = 6, y = -15, label = paste("n =", n), colour = "black", fontface = 'bold') +
        scale_color_manual(values = self$color_palette) +
        guides(color = guide_legend(override.aes = list(size = 6))) +
        theme(title = element_text(size = 15), 
              plot.title = element_text(hjust = .5),
              legend.position = "top") +
        labs(title = "Growth curve")

      return(plot)
    },

    boxplot_pvalue = function(step_increase, day) {
      df_p_val <- self$data %>% 
        wilcox_test(Totals ~ Concentration, ref.group = "Control") %>% 
        adjust_pvalue(p.col = "p", method = "bonferroni") %>%
        add_significance(p.col = "p.adj") %>%
        add_xy_position(step.increase = step_increase)

      Record_1 <- self$data %>% filter(Days == day)

      ggplot(Record_1, aes(x = Concentration, y = Totals)) +
        stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.4), width = 0.4) +
        geom_boxplot(aes(fill = Concentration), position = position_dodge(width = 0.1), outlier.shape = NA) +
        geom_jitter(width = 0.2, shape = 21, color = "grey20", size = 2, fill = "white", stroke = 1, show.legend = FALSE) +
        self$set_theme() +
        ylab("Cell numbers") +
        scale_fill_manual(values = self$color_palette) +
        scale_x_discrete(guide = "prism_bracket") +
        scale_y_continuous(guide = "prism_offset_minor") +
        stat_pvalue_manual(df_p_val, label = "p.adj.signif", 
                           label.size = 4, 
                           hide.ns = FALSE, 
                           bracket.size = 1) +
        theme(title = element_text(size = 15), plot.title = element_text(hjust = .5), legend.position = "none") +
        labs(title = "Growth variation")
    }
  )
)

# 创建绘图函数
create_plots <- function(data, concentrations, file_name_prefix, titles, color_palette, step_increase = 0.12, day = 6, plot_class, output_dir) {
  plot <- PlottingClass$new(data, var = concentrations, color_palette = color_palette)

  bubble_plot <- plot$bubble(title = titles)
  ggthemr::ggthemr("flat dark")
  
  smooth_plot <- plot$plot_smooth()
  box_plot <- plot$boxplot_pvalue(step_increase, day)

  if ("Times" %in% colnames(data)) {
    n_times <- length(unique(data$Times))
  } else {
    n_times <- 1
  }
  
  if (plot_class == "smooth_box") {
    combined_plot <- smooth_plot + box_plot + plot_layout(widths = c(5, 3))
    ggsave(paste0(file_name_prefix, "_smooth_box.pdf"), combined_plot, width = 35, height = 20, units = "cm")

  } else if (plot_class == "bubble_smooth_box") {
    box_plot <- box_plot + coord_flip()
    combined_plot <- bubble_plot / box_plot + plot_layout(heights = c(5, 3))
    ggsave(file.path(output_dir, paste0(file_name_prefix, "_bubble.pdf")), combined_plot, width = 35, height = 20, units = "cm")
    ggsave(file.path(output_dir, paste0(file_name_prefix, "_plot.pdf")), smooth_plot, width = 25, height = 20, units = "cm")
  }

  return(list(bubble_plot = bubble_plot, smooth_plot = smooth_plot, box_plot = box_plot))
}

# 创建命令行参数解析器
parser <- ArgumentParser(description = "绘制生长曲线图 (版本 1.0.0 - 最初开发版本)")
parser$add_argument("-d", "--data", required = TRUE, help = "输入数据文件路径 (CSV或Excel格式)")
parser$add_argument("--concentrations", required = TRUE, help = "浓度水平，逗号分隔")
parser$add_argument("--file_name_prefix", required = TRUE, help = "输出文件前缀")
parser$add_argument("--titles", required = TRUE, help = "图表标题")
parser$add_argument("--color_palette", required = TRUE, help = "颜色配置，逗号分隔")
parser$add_argument("--step_increase", type = "numeric", default = 0.12, help = "步长增加值")
parser$add_argument("--day", type = "numeric", default = 6, help = "指定的天数")
parser$add_argument("--plot_class", required = TRUE, help = "图形类别 ('smooth_box' 或 'bubble_smooth_box')")
parser$add_argument("--output_dir", default = ".", help = "输出文件夹路径 (默认为当前目录)")
parser$add_argument("--sheet", help = "Excel文件中要读取的工作表名称（可选）")  

args <- parser$parse_args()

# 读取数据
if (grepl("\\.csv$", args$data)) {
  data <- read.csv(args$data)
} else if (grepl("\\.xlsx$", args$data) || grepl("\\.xls$", args$data)) {
  if (is.null(args$sheet)) {
    data <- read_excel(args$data)  # 默认读取第一个工作表
  } else {
    data <- read_excel(args$data, sheet = args$sheet)  # 指定工作表读取
  }
} else {
  stop("错误：不支持的文件格式，仅支持 CSV 和 Excel 文件。")
}

# 解析颜色配置并转换为向量
color_palette <- strsplit(args$color_palette, ",")[[1]]

create_plots(data, strsplit(args$concentrations, ",")[[1]], args$file_name_prefix, args$titles, color_palette, args$step_increase, args$day, args$plot_class, args$output_dir)

# 结束信息
cat("绘图完成，输出文件保存在：", args$output_dir, "\n")
