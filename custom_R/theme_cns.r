theme_cns <- function() {
  
  # 设置整个图片的背景为空，这样即可获得透明背景
  theme(plot.background=element_blank(),
        #  这里设置线的宽度为0.5mm，按理来说0.5mm对应1.415pt，
        # 但不知为何，这里的0.5mm粗细的线条在AI中打开后为1.07pt。
        panel.grid=element_blank(),
        panel.background=element_blank(),
        panel.border=element_rect(color="black",linewidth=0.5,fill=NA),
        #  设置横纵坐标轴为空，
        axis.line=element_blank(),
        #  指定坐标轴刻度线为黑色，宽度0.5mm
        axis.ticks=element_line(color="black",linewidth=0.5),
        #  指定坐标轴的文字为黑色，字号7pt
        axis.text=element_text(color="black",size=7),
        #  指定坐标轴标题文字为黑色，字号7pt
        axis.text=element_text(color="black",size=7),
        #  指定图片标题为黑色，字号7pt
        plot.title=element_text(color="black",size=7),
        #  指定legend的背景设为空
        legend.background=element_blank(),
        legend.key=element_blank(),
        #  指定legend的字号和颜色
        legend.text=element_text(color="black",size=7),
        legend.title=element_text(color="black",size=7),
        )
}
