library(shiny)
library(readxl)
library(tidyverse)
library(ggthemr)
library(zip)
library(rstatix)
library(DT)  # 用于展示数据表格

# 定义 PlottingClass
PlottingClass <- R6::R6Class("PlottingClass",
  public = list(
    data = NULL,

    initialize = function(data, var) {
      self$data <- self$factor_reorder(data, var)
    },

    factor_reorder = function(data, var) {
      if (!"Concentration" %in% colnames(data)) {
        stop("错误：数据中缺少 'Concentration' 列。")
      }
      if (!is.character(var) || length(var) == 0) {
        stop("错误：var 参数必须是长度大于0的字符向量。")
      }
      if (!all(var %in% unique(data$Concentration))) {
        stop("错误：var 参数中的某些水平在 'Concentration' 列中不存在。")
      }
      data$Concentration <- factor(data$Concentration, levels = var)
      return(data)
    },

    set_theme = function() {
      ggplot2::theme_minimal() +
        theme(
          strip.text = element_text(size = 18),
          axis.line = element_line(color = "black", linewidth = 0.4),
          axis.text = element_text(color = "black", size = 15),
          axis.title = element_text(color = "black", size = 20),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.2, color = "#e5e5e5")
        )
    },

    bubble = function(title = "Penicillin treatment",
                      fill_color = c("red", "green", "#9ED8DB", "#3498DB")) {
      max_days <- max(self$data$Days, na.rm = TRUE)
      plot <- ggplot(self$data, aes(x = Days, y = Concentration)) +
        geom_point(aes(colour = Concentration, size = Totals), alpha = .7) +
        geom_text(aes(label = Totals), vjust = -1.1, size = 4, fontface = 'bold') +
        scale_size_continuous(range = c(3, 10)) +
        labs(title = title) +
        guides(colour = 'none', size = 'none') +
        self$set_theme() +
        scale_color_manual(values = fill_color) + 
        scale_x_continuous(breaks = seq(0, max_days + 4, by = 1)) +
        facet_grid(reformulate(ifelse("Times" %in% colnames(self$data), "Times", "Repeat")), scales = 'free')

      return(plot)
    },

    plot_smooth = function(days_col = "Days", concentration_col = "Concentration", 
                           totals_col = "Totals", fill_color = c("#3498DB", "black")) {
      Record_1 <- self$data %>%
        group_by(!!sym(concentration_col), !!sym(days_col)) %>%
        summarise(mean_totals = mean(!!sym(totals_col), na.rm = TRUE), .groups = 'drop')

      ggplot(Record_1, aes(x = !!sym(days_col), y = mean_totals, color = !!sym(concentration_col))) +
        geom_smooth(se = TRUE, alpha = .05) +
        geom_point(size = 3, alpha = .8) +
        self$set_theme() +
        scale_color_manual(values = fill_color)
    },

    boxplot_pvalue = function(ref_group = "Control", fill_color = c("red", "green", "#9ED8DB", "#3498DB")) {
      df_p_val <- self$data %>%
        wilcox_test(Totals ~ Concentration, ref.group = ref_group)

      ggplot(self$data, aes(x = Concentration, y = Totals, fill = Concentration)) +
        geom_boxplot() +
        self$set_theme() +
        scale_fill_manual(values = fill_color)
    }
  )
)

# 创建示例数据
create_demo_data <- function() {
  data.frame(
    Concentration = c("Control", "100 U", "200 U", "500 U"),
    Days = rep(1:6, each = 4),
    Totals = sample(1:100, 24, replace = TRUE),
    Repeat = rep(1:3, each = 8),
    Times = rep(1:2, each = 12)
  )
}

# 定义UI
ui <- fluidPage(
  titlePanel("绘图应用"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "上传数据文件 (CSV 或 Excel)", 
                accept = c(".csv", ".xlsx")),
      textInput("concentrations", "输入浓度 (用逗号分隔)", value = "浓度1,浓度2"),
      textInput("file_name_prefix", "文件名前缀", value = "plot"),
      textInput("titles", "图表标题", value = "Penicillin treatment"),
      numericInput("step_increase", "步长增加", value = 0.12, step = 0.01),
      numericInput("day", "选择的天数", value = 6),
      selectInput("plot_class", "选择图表类型", 
                  choices = c("smooth_box", "bubble_smooth_box")),
      uiOutput("sheet_ui"),  # 动态生成的Sheet选择UI
      numericInput("color_count", "输入颜色数量", value = 4, min = 1, max = 10),
      uiOutput("color_inputs"),  # 动态生成的颜色输入框
      numericInput("image_width", "图像宽度 (英寸)", value = 10, min = 1),  # 图像宽度输入框
      numericInput("image_height", "图像高度 (英寸)", value = 6, min = 1),  # 图像高度输入框
      actionButton("generate", "生成图表"),
      actionButton("demo", "使用演示数据")  # 演示数据按钮
    ),
    
    mainPanel(
      plotOutput("bubblePlot"),
      plotOutput("smoothPlot"),
      plotOutput("boxPlot"),
      downloadButton("download", "下载图表"),
      textOutput("error_message"),  # 显示错误消息
      DTOutput("dataTable")  # 数据表格输出
    )
  )
)

# 定义服务器逻辑
server <- function(input, output, session) {
  
  sheets <- reactive({
    req(input$file)
    if (grepl("\\.xlsx$", input$file$name)) {
      excel_sheets(input$file$datapath)
    } else {
      NULL
    }
  })

  # 动态生成Sheet选择下拉菜单
  output$sheet_ui <- renderUI({
    req(sheets())
    selectInput("sheet", "选择Sheet", choices = sheets())
  })

  # 动态生成颜色输入框
  output$color_inputs <- renderUI({
    lapply(1:input$color_count, function(i) {
      textInput(paste0("color", i), paste("颜色", i), value = ifelse(i <= 4, c("red", "green", "#9ED8DB", "#3498DB")[i], "#FFFFFF"))
    })
  })

  # 准备数据
  prepare_data <- function(data) {
    output$dataTable <- renderDT({ datatable(data) })
    return(data)
  }

  plots <- eventReactive(input$generate, {
    req(input$file)
    data <- if (grepl("\\.xlsx$", input$file$name)) {
      req(input$sheet)
      read_excel(input$file$datapath, sheet = input$sheet)
    } else {
      read.csv(input$file$datapath)
    }
    
    data <- prepare_data(data)

    concentrations <- gsub("，", ",", input$concentrations)  # 替换中文逗号为英文逗号
    if (trimws(concentrations) == "") {
      output$error_message <- renderText("错误：浓度输入不能为空！")
      return(NULL)
    } else {
      output$error_message <- renderText("")  # 清空错误消息
    }

    concentrations <- unlist(strsplit(concentrations, ","))

    fill_colors <- sapply(1:input$color_count, function(i) input[[paste0("color", i)]])
    plot_obj <- PlottingClass$new(data, var = concentrations)

    # 生成图表
    return(list(
      bubble_plot = plot_obj$bubble(title = input$titles, fill_color = fill_colors),
      smooth_plot = plot_obj$plot_smooth(fill_color = fill_colors),
      box_plot = plot_obj$boxplot_pvalue(ref_group = "Control", fill_color = fill_colors)
    ))
  })

  # 使用演示数据的逻辑
  observeEvent(input$demo, {
    demo_data <- create_demo_data()
    prepare_data(demo_data)
    concentrations <- unique(demo_data$Concentration)
    plot_obj <- PlottingClass$new(demo_data, var = concentrations)
    fill_colors <- c("red", "green", "#9ED8DB", "#3498DB")

    output$bubblePlot <- renderPlot({ plot_obj$bubble(fill_color = fill_colors) })
    output$smoothPlot <- renderPlot({ plot_obj$plot_smooth(fill_color = fill_colors) })
    output$boxPlot <- renderPlot({ plot_obj$boxplot_pvalue(fill_color = fill_colors) })
  })

  output$bubblePlot <- renderPlot({
    req(plots())
    plots()$bubble_plot
  })

  output$smoothPlot <- renderPlot({
    req(plots())
    plots()$smooth_plot
  })

  output$boxPlot <- renderPlot({
    req(plots())
    plots()$box_plot
  })

  output$download <- downloadHandler(
    filename = function() {
      paste(input$file_name_prefix, "plots.zip", sep = "_")
    },
    content = function(file) {
      zip_file <- tempfile(fileext = ".zip")
      files_to_save <- c(
        paste0(input$file_name_prefix, "_smooth_box.pdf"),
        paste0(input$file_name_prefix, "_bubble.pdf"),
        paste0(input$file_name_prefix, "_plot.pdf")
      )

      # 使用用户输入的宽度和高度保存图像
      ggsave(paste0(input$file_name_prefix, "_smooth_box.pdf"), plots()$smooth_plot, 
             width = input$image_width, height = input$image_height, units = "in", dpi = 300)
      ggsave(paste0(input$file_name_prefix, "_bubble.pdf"), plots()$bubble_plot, 
             width = input$image_width, height = input$image_height, units = "in", dpi = 300)
      ggsave(paste0(input$file_name_prefix, "_plot.pdf"), plots()$box_plot, 
             width = input$image_width, height = input$image_height, units = "in", dpi = 300)

      zip::zip(zip_file, files = files_to_save)
      file.rename(zip_file, file)
    }
  )
}

# 运行Shiny应用
shinyApp(ui = ui, server = server)
