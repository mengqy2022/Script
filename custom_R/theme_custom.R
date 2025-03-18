#  ggplot2自定义主题
#  2024/01/08
#  mqy


theme_custom_1 <- function (base_size = 12, base_family = "Arial") {
  
  half_line <- base_size/2
  
  theme(
    line = element_line(color = "black", linewidth = .5,
                        linetype = 1, lineend = "butt"),
    rect = element_rect(fill = "white", color = "black",
                        linewidth = .5, linetype = 1),
    text = element_text(family = base_family, face = "plain",
                        color = "black", size = base_size,
                        lineheight = .9, hjust = .5, vjust = .5,
                        angle = 0, margin = margin(), debug = FALSE),
    axis.line = element_blank(),
    axis.line.x = NULL,
    axis.line.y = NULL,
    axis.text = element_text(size = base_size * 1.1, color = "gray30"),
    axis.text.x = element_text(margin = margin(t = .8 * half_line/2),
                               vjust = 1),
    axis.text.x.top = element_text(margin = margin(b = .8 * half_line/2),
                                   vjust = 0),
    axis.text.y = element_text(margin = margin(r = .8 * half_line/2),
                               hjust = 1),
    axis.text.y.right = element_text(margin = margin(l = .8 * half_line/2),
                                     hjust = 0),
    axis.ticks = element_line(color = "gray30", linewidth = .7),
    axis.ticks.length = unit(half_line / 1.5, "pt"),
    axis.ticks.length.x = NULL,
    axis.ticks.length.x.top = NULL,
    axis.ticks.length.x.bottom = NULL,
    axis.ticks.length.y = NULL,
    axis.ticks.length.y.left = NULL,
    axis.ticks.length.y.right = NULL,
    axis.title.x = element_text(margin = margin(t = half_line),
                                vjust = 1, size = base_size * 1.3,
                                face = "bold"),
    axis.title.x.top = element_text(margin = margin(b = half_line),
                                    vjust = 0),
    axis.title.y = element_text(angle = 90, vjust = 1,
                                margin = margin(r = half_line),
                                size = base_size * 1.3, face = "bold"),
    axis.title.y.right = element_text(angle = -90, vjust = 0,
                                      margin = margin(l = half_line)),
    legend.background = element_rect(color = NA),
    legend.spacing = unit(.4, "cm"),
    legend.spacing.x = NULL,
    legend.spacing.y = NULL,
    legend.margin = margin(.2, .2, .2, .2, "cm"),
    legend.key = element_rect(fill = "gray95", color = "white"),
    legend.key.size = unit(1.2, "lines"),
    legend.key.height = NULL,
    legend.key.width = NULL,
    legend.text = element_text(size = rel(.8)),
    legend.text.align = NULL,
    legend.title = element_text(hjust = 0),
    legend.title.align = NULL,
    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    legend.box = NULL,
    legend.box.margin = margin(0, 0, 0, 0, "cm"),
    legend.box.background = element_blank(),
    legend.box.spacing = unit(.4, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "gray30",
                                fill = NA, linewidth = .7),
    panel.grid.major = element_line(color = "gray90", linewidth = 1),
    panel.grid.minor = element_line(color = "gray90", linewidth = .5,
                                    linetype = "dashed"),
    panel.spacing = unit(base_size, "pt"),
    panel.spacing.x = NULL,
    panel.spacing.y = NULL,
    panel.ontop = FALSE,
    strip.background = element_rect(fill = "white", color = "gray30"),
    strip.text = element_text(color = "black", size = base_size),
    strip.text.x = element_text(margin = margin(t = half_line,
                                                b = half_line)),
    strip.text.y = element_text(angle = -90,
                                margin = margin(l = half_line,
                                                r = half_line)),
    strip.text.y.left = element_text(angle = 90),
    strip.placement = "inside",
    strip.placement.x = NULL,
    strip.placement.y = NULL,
    strip.switch.pad.grid = unit(0.1, "cm"),
    strip.switch.pad.wrap = unit(0.1, "cm"),
    plot.background = element_rect(color = NA),
    plot.title = element_text(size = base_size * 1.8, hjust = .5,
                              vjust = 1, face = "bold",
                              margin = margin(b = half_line * 1.2)),
    plot.title.position = "panel",
    plot.subtitle = element_text(size = base_size * 1.3,
                                 hjust = .5, vjust = 1,
                                 margin = margin(b = half_line * .9)),
    plot.caption = element_text(size = rel(0.9), hjust = 1, vjust = 1,
                                margin = margin(t = half_line * .9)),
    plot.caption.position = "panel",
    plot.tag = element_text(size = rel(1.2), hjust = .5, vjust = .5),
    plot.tag.position = "topleft",
    plot.margin = margin(base_size, base_size, base_size, base_size),
    complete = TRUE
  )
}

theme_custom_2 <- theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "none", # 删除图例
        axis.text.x = element_text(size=10,angle=70,hjust=1),
        axis.text.y = element_text(size=10),
        panel.grid.major.y = element_line(color=1,size=0.4,linetype=2),
        panel.grid.minor.y = element_line(color="black",size=0.25,linetype=3),
        strip.text = element_text(size = 12),
        strip.placement = 'outside'
  )

theme_custom_3 <- theme_prism( border = TRUE,base_size = 5)+
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

theme_custom_4 <- function(){
  theme_test() +
    theme(axis.title.x = element_blank(),
          axis.line = element_line(color = "#3D4852"),
          axis.ticks = element_line(color = "#3D4852"),
          panel.grid.major.y = element_line(color = "#DAE1E7"),
          panel.grid.major.x = element_blank(),
          plot.margin = unit(rep(0.2,4),"cm"),
          axis.text = element_text(size = 10, color = "#22292F"),
          axis.title = element_text(size = 10, hjust = 1),
          axis.title.y = element_blank(),
          axis.text.y = element_text(margin = margin(r = 5)),
          axis.text.x = element_text(margin = margin(t = 5)),
          legend.position = "non")
}