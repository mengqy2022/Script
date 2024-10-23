file_handle <- function(data, sample) {
  
  data <- data[,1:(1+sample*3)]
  data <- data %>% pivot_longer(!Cycle,names_to = "sample_id", values_to = "RFU")
  lab_data <- data %>% group_by(sample_id) %>% summarise(sample_id = unique(sample_id), max = max(RFU))

  #  基线计算
  baseline <<- list()
  
  for (i in 3:15) {
    baseline[i] <- list(data[which(data$Cycle == i),])
  }
  
  baseline <- do.call(rbind, baseline)
  baseline_value <<- baseline %>% group_by(sample_id) %>% summarise(RFU = sum(RFU)/length(3:15))
  baseline_value <<- sum(baseline_value$RFU)/sample*3
  
  #  阈值
  threshold <<- sd(baseline$RFU)*10
  
  #  绘图
  ggplot(data = data, aes(x = Cycle, y = RFU,colour = sample_id)) +
    geom_line(linewidth = 1,show.legend = F)+
    geom_hline(aes(yintercept = threshold),size = 1.5,color = "firebrick") +
    geom_hline(aes(yintercept = baseline_value),size = 1,color = "black", linetype = "dashed")+
    #annotate("text", x = 40, y = 559.6111, label = "A1") +
    #geom_label_repel(data = lab_data, aes(label = max),show.legend = F)+
    scale_x_continuous(expand = c(0,1))+
    #scale_y_continuous(expand = c(0,0))+
    theme_bw()
  
}

Calculation_process <- function(data, int_number, control_number) {
  
  #  设置内参
  internal_reference <<- data[(seq(1, dim(data)[1], by =3)[int_number]:(int_number*3)),] %>% arrange(Cq)
  
  #  设置对照
  # 对照CQ均值：平均值（对照-内置）
  Control <<- data[(seq(1, dim(data)[1], by =3)[control_number]):(control_number*3),] %>% arrange(Cq)
  Control$Cq <<- Control$Cq - internal_reference$Cq
  Control_mean <<- mean(Control$Cq)
  
  #  实验组
  Experimental_1 <- CQ[-((seq(1, dim(CQ)[1], by =3)[control_number]):(control_number*3)),]
  Experimental_2 <- CQ[-((seq(1, dim(CQ)[1], by =3)[int_number]:(int_number*3))),]
  Experimental <<- Experimental_1 %>% inner_join(Experimental_2)
  
  # 计算样本基因量
  
  data_result <<- data.frame()
  
  for (i in 1:(dim(Experimental)[1]/3)) {
    
    start <- (seq(1, dim(Experimental)[1], by =3))[i]
    
    a <- Experimental[c(start:(i*3)),]
    
    a <- a %>% arrange(Cq)
    
    # 每个样本CQ减去内参CQ
    a[2] <- a[2] - internal_reference[2]
    # 再减去对照组的平均CQ
    a[2] <- a[2] - Control_mean
    # 归一化
    a[2] <- 2 ^ -a[2]

    data_result <<- data_result %>% rbind(a) %>% arrange(Sample)
    
  }
  
  # 计算每组均值
  sample_mean <- data.frame()
  groups_mean <<- data.frame()
  
  for (i in 1:(dim(data_result)[1]/3)) {
    
    b <- 3*i
    
    start <- (seq(1, dim(data_result)[1], by =3))[i]
    
    a <- data_result[c(start:(i*3)),]
    
    a[2] <- mean(a$Cq)
    
    sample_mean <- sample_mean %>% rbind(a)
    
    print(sample_mean[(seq(1, dim(data_result)[1], by =3))[i],2])
    
    b <- sample_mean[(seq(1, dim(data_result)[1], by =3))[i],]
    
    groups_mean <<- groups_mean %>% rbind(b)
    
  }
}

sample_six_plot_bar <- function(data) {
  
  ggplot(data, aes(Sample, Cq, fill = Sample))+
    geom_bar(stat="summary", fun=mean, position="dodge", color = "black") +
    theme_prism(palette = "candy_bright",
                base_fontface = "plain", # 字体样式，可选 bold, plain, italic
                base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
                base_size = 16,  # 图形的字体大小
                base_line_size = 0.8, # 坐标轴的粗细
                axis_text_angle = 0) + # 可选值有 0，45，90，270
    geom_jitter(width = 0.2,shape = 21, color="grey20",size=2, fill="white", size=2, stroke=1,show.legend = F) +
    stat_summary(geom = "errorbar", fun.data = 'mean_sd', width = 0.2, show.legend = F)+
    scale_fill_brewer() + 
    ylab("Quantification of symbiotic bacteria") +
    geom_signif(comparisons = list(c("DFMO_1", "size_1"),
                                   c("size_1", "size_2"),
                                   c("size_2", "size_3"),
                                   c("size_3", "size_4"),
                                   
                                   c("DFMO_1", "size_2"),
                                   c("size_2", "size_4"),
                                   
                                   c("DFMO_1", "size_3"),
                                   
                                   c("size_1", "size_3"),
                                   
                                   c("size_1", "size_4"),
                                   
                                   c("DFMO_1", "size_4")),
                map_signif_level = T,
                y_position = c(1.5,1.5,1.5,1.5, 1.6,1.6, 1.7, 1.8, 1.9, 2.0),
                size=0.8,
                color="black")
  
}

Calculation_process_unicellular <- function(data, int_number) {
  
  #  设置内参
  internal_reference <<- data[(seq(1, dim(data)[1], by =3)[int_number]:(int_number*3)),] %>% arrange(Cq)
  internal_reference_mean <<- mean(internal_reference$Cq)

  Experimental <- data

  # 计算样本基因量

  data_result <<- data.frame()

  for (i in 1:(dim(Experimental)[1]/3)) {

    start <- (seq(1, dim(Experimental)[1], by =3))[i]

    a <- Experimental[c(start:(i*3)),]

    a <- a %>% arrange(Cq)

    # 再减去内参的平均CQ
    a[2] <- a[2] - internal_reference_mean
    # 归一化
    a[2] <- 2 ^ -a[2]

    data_result <<- data_result %>% rbind(a) %>% arrange(Sample)

  }
  
  # 计算每组均值
  sample_mean <- data.frame()
  groups_mean <<- data.frame()
  
  for (i in 1:(dim(data_result)[1]/3)) {
    
    b <- 3*i
    
    start <- (seq(1, dim(data_result)[1], by =3))[i]
    
    a <- data_result[c(start:(i*3)),]
    
    a[2] <- mean(a$Cq)
    
    sample_mean <- sample_mean %>% rbind(a)
    
    print(sample_mean[(seq(1, dim(data_result)[1], by =3))[i],2])
    
    b <- sample_mean[(seq(1, dim(data_result)[1], by =3))[i],]
    
    groups_mean <<- groups_mean %>% rbind(b)
    
  }
}

sample_five_plot_bar <- function(data) {
  
  ggplot(data, aes(Sample, Cq, fill = Sample))+
    geom_bar(stat="summary", fun=mean, position="dodge", color = "black") +
    theme_prism(palette = "candy_bright",
                base_fontface = "plain", # 字体样式，可选 bold, plain, italic
                base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
                base_size = 13,  # 图形的字体大小
                base_line_size = 0.8, # 坐标轴的粗细
                axis_text_angle = 0) + # 可选值有 0，45，90，270
    geom_jitter(width = 0.2,shape = 21, color="grey20",size=2, fill="white", size=2, stroke=1,show.legend = F) +
    stat_summary(geom = "errorbar", fun.data = 'mean_sd', width = 0.2, show.legend = F)+
    scale_fill_brewer() + 
    ylab("Quantification of symbiotic bacteria") +
    xlab("") +
    geom_signif(comparisons = list(#c("DFMO_1", "size_1"),
                                   c("size_1", "size_2"),
                                   c("size_2", "size_3"),
                                   c("size_3", "size_4"),
                                   c("size_4", "size_5"),
                                   
                                   #c("DFMO_1", "size_2"),
                                   c("size_2", "size_4"),
                                  
                                   c("size_2", "size_5"),
                                   
                                   c("size_1", "size_4"),
                                   
                                   c("size_1", "size_3"),
                                   c("size_3", "size_5"),
                                   
                                   #c("DFMO_1", "size_3"),
                                   
                                   #c("DFMO_1", "size_4"),
                                   
                                   c("size_1", "size_5")), 
                map_signif_level = T, 
                #y_position = c(4.4,4.4,4.4,4.4,4.4, 4.6,4.6, 4.8, 5.0, 5.2, 5.4,5.4, 5.6, 5.8, 6.0),
                y_position = c(4.4,4.4,4.4,4.4, 4.6, 4.8, 5.0, 5.2,5.2, 5.4),
                #y_position = c(2.4,2.4,2.4,2.4, 2.6, 2.8, 3.0, 3.2, 3.4,3.4, 3.6),
                size=0.8,
                textsize = 3,
                #tip_length = 2,
                color="black")
  
}

sample_three_plot_bar <- function(data) {

  ggplot(data, aes(Sample, Cq, fill = Sample))+
    geom_bar(stat="summary", fun=mean, position="dodge", color = "black") +
    theme_prism(palette = "candy_bright",
                base_fontface = "plain", # 字体样式，可选 bold, plain, italic
                base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
                base_size = 13,  # 图形的字体大小
                base_line_size = 0.8, # 坐标轴的粗细
                axis_text_angle = 0) + # 可选值有 0，45，90，270
    geom_jitter(width = 0.2,shape = 21, color="grey20",size=2, fill="white", size=2, stroke=1,show.legend = F) +
    stat_summary(geom = "errorbar", fun.data = 'mean_sd', width = 0.2, show.legend = F)+
    scale_fill_brewer() + 
    ylab("Quantification of symbiotic bacteria") +
    xlab("") +
    geom_signif(comparisons = list(
      c("size_1", "size_2"),
      c("size_2", "size_3"),
      c("size_1", "size_3")), 
      map_signif_level = T, 
      y_position = c(8.4,8.4, 9),
      size=0.8,
      textsize = 3,
      color="black")

}

