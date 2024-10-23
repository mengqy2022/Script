theme_custom_1 <- theme_prism( border = TRUE,base_size = 5)+
  theme(strip.text.x= element_text(size = 8),
        ##'@X轴和Y轴字体标题大小
        title = element_text(size = 8),
        legend.box.spacing = unit(1,"cm"),
        ##'@右上标签字体大小
        legend.text = element_text(size = 6),
        ##'@右上标签标题大小
        legend.title = element_text(size = 8),
        ##'@X轴和Y轴字体大小
        axis.text.y = element_text(size = 6, angle = 0, vjust = 0.2),
        axis.text.x = element_text(size = 6,angle = 45),
        panel.grid = element_line(color = "gray",
                                  linewidth = 0.15,
                                  linetype = 2),
        panel.spacing = unit(1, "lines"),
        plot.caption = element_text(size = 8)
  )