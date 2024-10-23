theme_custom_2 <- function(){
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