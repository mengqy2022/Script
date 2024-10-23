library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)

# 定义一个类
PCRAnalysis <- R6::R6Class("PCRAnalysis",
  public = list(
    data = NULL,
    sample = NULL,
    internal_reference = NULL,
    threshold = NULL,
    baseline_value = NULL,
    groups_mean = NULL,
    
    # 构造函数
    initialize = function(data, sample) {
      self$data <- data
      self$sample <- sample
    },
    
    # 文件处理函数
    file_handle = function() {
      data <- self$data[, 1:(1 + self$sample * 3)]
      data <- data %>% pivot_longer(!Cycle, names_to = "sample_id", values_to = "RFU")
      lab_data <- data %>% group_by(sample_id) %>% summarise(sample_id = unique(sample_id), max = max(RFU))

      # 基线计算
      baseline <- list()
      
      for (i in 3:15) {
        baseline[i] <- list(data[which(data$Cycle == i),])
      }
      
      baseline <- do.call(rbind, baseline)
      self$baseline_value <- baseline %>% group_by(sample_id) %>% summarise(RFU = sum(RFU) / length(3:15))
      self$baseline_value <- sum(self$baseline_value$RFU) / self$sample * 3
      
      # 阈值
      self$threshold <- sd(baseline$RFU) * 10
      
      # 绘图
      ggplot(data = data, aes(x = Cycle, y = RFU, colour = sample_id)) +
        geom_line(linewidth = 1, show.legend = F) +
        geom_hline(aes(yintercept = self$threshold), size = 1.5, color = "firebrick") +
        geom_hline(aes(yintercept = self$baseline_value), size = 1, color = "black", linetype = "dashed") +
        scale_x_continuous(expand = c(0, 1)) +
        theme_bw()
    },
    
    # 计算过程
    calculation_process = function(int_number, control_number) {
      # 设置内参
      self$internal_reference <- self$data[(seq(1, dim(self$data)[1], by = 3)[int_number]:(int_number * 3)),] %>% arrange(Cq)

      # 设置对照
      Control <- self$data[(seq(1, dim(self$data)[1], by = 3)[control_number]):(control_number * 3),] %>% arrange(Cq)
      Control$Cq <- Control$Cq - self$internal_reference$Cq
      Control_mean <- mean(Control$Cq)

      # 实验组
      Experimental_1 <- self$data[-((seq(1, dim(self$data)[1], by = 3)[control_number]):(control_number * 3)),]
      Experimental_2 <- self$data[-((seq(1, dim(self$data)[1], by = 3)[int_number]:(int_number * 3))),]
      Experimental <- Experimental_1 %>% inner_join(Experimental_2)

      # 计算样本基因量
      data_result <- data.frame()

      for (i in 1:(dim(Experimental)[1] / 3)) {
        start <- (seq(1, dim(Experimental)[1], by = 3))[i]
        a <- Experimental[c(start:(i * 3)),]
        a <- a %>% arrange(Cq)

        # 每个样本CQ减去内参CQ
        a[2] <- a[2] - self$internal_reference[2]
        # 再减去对照组的平均CQ
        a[2] <- a[2] - Control_mean
        # 归一化
        a[2] <- 2 ^ -a[2]

        data_result <- data_result %>% rbind(a) %>% arrange(Sample)
      }

      # 计算每组均值
      sample_mean <- data.frame()
      self$groups_mean <- data.frame()

      for (i in 1:(dim(data_result)[1] / 3)) {
        start <- (seq(1, dim(data_result)[1], by = 3))[i]
        a <- data_result[c(start:(i * 3)),]
        a[2] <- mean(a$Cq)

        sample_mean <- sample_mean %>% rbind(a)
        self$groups_mean <- self$groups_mean %>% rbind(a)
      }
    },

    # 其他绘图方法及计算方法可依次添加...
  )
)
