# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Email     : 15877464851@163.com
# @Time      : 2024/09/23 10:00

library(R6)
library(tidyverse)
library(ggrepel)
library(ggprism)
library(rstatix)
library(patchwork)
library(ggpubr)

PlottingClass <- R6Class("PlottingClass",
  public = list(
    data = NULL,

    # 初始化函数
    initialize = function(data, var) {
      self$data <- self$factor_reorder(data, var)
    },
    
    
    
    # 因子重排
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
      

      # 按数据输入的变量顺序排序
      data$Concentration <- factor(data$Concentration, levels = unique(var))
      return(data)
    },

    # 设置主题
    set_theme = function() {
      theme_prism() +
        theme(strip.text = element_text(size = 18),
              axis.line = element_line(color = "black", linewidth = 0.4),
              axis.text = element_text(color = "black", size = 15),
              axis.title = element_text(color = "black", size = 20),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_line(linewidth = 0.2, color = "#e5e5e5"))
    },

    # 绘制泡泡图
    bubble = function(title = "Penicillin treatment",
                      fill_color = c("red", "green", "#9ED8DB", "#3498DB")) {
      max_days <- max(self$data$Days, na.rm = TRUE)
      plot <- ggplot(self$data, aes(x = Days, y = Concentration)) +
        geom_point(aes(colour = Concentration, size = Totals), alpha = .7) +
        geom_text(aes(label = Totals), vjust = -1.1, size = 4, fontface = 'bold') +
        scale_size_continuous(range = c(3, 10)) +
        labs(title = title) +
        guides(colour = 'none', size = 'none') +
        theme(title = element_text(size = 15), plot.title = element_text(hjust = .5),
              axis.text = element_text(color = "black", size = 15),
              axis.title = element_text(color = "black", size = 20)) +
        scale_color_manual(values = fill_color) +  # 使用自定义颜色
        scale_x_continuous(breaks = seq(0, max_days + 4, by = 1))

      if ("Times" %in% colnames(self$data)) {
        plot <- plot + facet_grid(Times ~ Repeat, scales = 'free') +
          theme(strip.text.x = element_text(size = 15),
                strip.text.y = element_text(size = 15))
      } else {
        plot <- plot + facet_grid(Repeat, scales = 'free')
      }

      return(plot)
    },

    # 绘制平滑图
    plot_smooth = function(days_col = "Days", concentration_col = "Concentration", 
                           totals_col = "Totals", times_col = "Times",
                           fill_color = c("red", "green", "#9ED8DB", "#3498DB")) {
      # 聚合计算均值
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

      max_days <- max(Record_1[[days_col]], na.rm = TRUE)

      plot <- ggplot(Record_1, aes(x = !!sym(days_col), y = mean_totals, color = !!sym(concentration_col))) +
        geom_smooth(aes(linetype = !!sym(concentration_col)), se = TRUE, alpha = .05) +
        geom_point(size = 3, alpha = .8) +
        self$set_theme() +
        scale_x_continuous(breaks = seq(0, max_days, by = 1)) +
        geom_label_repel(data = best_in_class, aes(label = mean_totals), show.legend = FALSE) +
        geom_label_repel(data = start, aes(label = mean_totals), show.legend = FALSE) +
        xlab("Days") +
        ylab("Cell numbers") +
        annotate("text", x = max_days, y = -15, label = paste("n =", n), colour = "black", fontface = 'bold') +
        scale_color_manual(values = fill_color) + 
        guides(color = guide_legend(override.aes = list(size = 6), ncol = 6)) +
        theme(title = element_text(size = 15), 
              plot.title = element_text(hjust = .5),
              legend.position = "top") +
        labs(title = "Growth curve")

      return(plot)
    },

    # 绘制箱线图
    boxplot_pvalue = function(step_increase, days_col = "Days", ref_group = "Control", fill_color = c("red", "green", "#9ED8DB", "#3498DB")) {
      
      max_days <- max(self$data[[days_col]], na.rm = TRUE)
      
      Record_1 <- self$data %>% filter(Days == max_days)

      df_p_val <- Record_1 %>%
        t_test(Totals ~ Concentration, ref.group = ref_group, paired = TRUE) %>%
        adjust_pvalue(p.col = "p", method = "bonferroni") %>%
        add_significance(p.col = "p.adj") %>%
        add_xy_position(step.increase = step_increase)
      
      # 检查 Totals 是否都大于 0
      all_positive <- all(Record_1$Totals > 0)
      
      # 创建 ggplot 对象
      p <- ggplot(Record_1, aes(x = Concentration, y = Totals)) +
        stat_boxplot(geom = "errorbar", position = position_dodge(width = 0.4), width = 0.4) +
        geom_boxplot(aes(fill = Concentration), position = position_dodge(width = 0.1), outlier.shape = NA) +
        geom_jitter(width = 0.2, shape = 21, color = "grey20", size = 2, fill = "white", stroke = 1, show.legend = FALSE) +
        self$set_theme() +
        ylab("Cell numbers") +
        scale_fill_manual(values = fill_color) +
        scale_x_discrete(guide = "prism_bracket") +
        stat_pvalue_manual(df_p_val, label = "p.adj.signif", 
                          label.size = 4, 
                          hide.ns = FALSE, 
                          bracket.size = 1) +
        guides(color = guide_legend(override.aes = list(size = 6), ncol = 6)) +
        theme(title = element_text(size = 15), 
              plot.title = element_text(hjust = .5), 
              legend.position = "top") +
        labs(title = "Growth variation")
      
      # 根据检查结果调整 y 轴
      if (all_positive) {
        p <- p + scale_y_continuous(limits = c(0, NA), guide = "prism_offset_minor")
      } else {
        p <- p + scale_y_continuous(guide = "prism_offset_minor")
      }
      
      # 打印或绘制图形
      print(p)
    }
  )
)

create_plots <- function(data, concentrations, file_name_prefix, titles, 
                         step_increase = 0.12, 
                         plot_class, 
                         color = c("red", "green", "#9ED8DB", "#3498DB"), 
                         ref_group = "Control") { 
  # 创建绘图对象
  plot <- PlottingClass$new(data, var = concentrations)

  # 生成泡泡图
  bubble_plot <- plot$bubble(title = titles, fill_color = color)
  ggthemr::ggthemr("flat dark")
  
  # 生成平滑图
  smooth_plot <- plot$plot_smooth(fill_color = color)

  # 运行统计测试并生成箱线图，传递 ref_group 参数
  box_plot <- plot$boxplot_pvalue(step_increase, fill_color = color, ref_group = ref_group)

  # 保存图片
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
    ggsave(paste0(file_name_prefix, "_bubble.pdf"), combined_plot, width = 35, height = 25, units = "cm")
    ggsave(paste0(file_name_prefix, "_plot.pdf"), smooth_plot, width = 20, height = 20, units = "cm")
  }

  return(list(bubble_plot = bubble_plot, smooth_plot = smooth_plot, box_plot = box_plot))
}


sort_size <- function(data) {
  factor_sort <- reorder(data$Concentration, data$Totals, FUN = mean)
  factor_sort <- unique(factor(factor_sort, levels = rev(levels(factor_sort))))
  factor_sort <- trimws(as.character(unique(levels(factor_sort))))
  return(factor_sort)
}


