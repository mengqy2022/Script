#############################
#目的函数，绘制火山图
#直接输入差异分析的结果中的logFC以及P\adj.P的列名
#cut_P、cut_FC为阈值
#top为绘制基因的数目
#x.num用于调整xlim,y.num
#highlight,高亮你关注的基因，输入为基因名或者list
#orderByFC等于TRUE时绘制,按logFC排序，FALSE时按Pvalue排序
#background默认颜色为c("#74add1","#d73027"),去除设定为null即可
#############################

# deg <- read.table("DEGs_mRNA_all_TCGA.csv",sep = ",",check.names = F,
#                   stringsAsFactors = F,header = T,row.names = 1)

suppressMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggfun)
  library(grid)
})
Volcanol_1<-function(deg,
                     logFC="log2FoldChange",
                     p="padj",
                     cut_P=0.05,
                     cut_FC=1,
                     x.num=0,
                     y.num=0,
                     top=10,
                     highlight=NULL,
                     orderByFC=TRUE,
                     background=c("#74add1","#d73027")){
  
  deg$Symbol<-rownames(deg)
  deg$P<-deg[,p]
  deg$FC<-deg[,logFC]
  deg$change = ifelse(deg$P>cut_P,'Stable', 
                      ifelse( deg$FC >cut_FC,'Up', 
                              ifelse( deg$FC< -cut_FC ,'Down','Stable') ))
  max_lfc<-max(deg[deg$change=="Up",logFC])
  min_lfc<-min(deg[deg$change=="Down",logFC])
  if(max_lfc>3){max_lfc<- -max_lfc}
  if(min_lfc< -3){min_lfc <- -min_lfc}
  if(y.num==0){
    y.num=max(-log10(deg$P))
  }
  if(x.num==0){
    x.num=max(abs(deg$FC))
  }
  if(is.null(background)){
    if(orderByFC==TRUE){
      Up1<-deg %>% 
        tidyr::drop_na() %>%
        dplyr::filter(change != "Stable") %>%
        dplyr::arrange(desc(FC)) %>%
        dplyr::slice(1:(top/2)) %>%
        dplyr::filter(change == "Up")    
      Down1<-deg %>% 
        tidyr::drop_na() %>%
        dplyr::filter(change != "Stable") %>%
        dplyr::arrange(desc(FC)) %>%
        dplyr::slice((dim(.)[1]-(top/2-1)):dim(.)[1]) %>%
        dplyr::filter(change == "Down")
      
      p1<-ggplot(data = deg) + 
        geom_point(aes(x = FC, y = -log10(P), 
                       color = FC,
                       size = -log10(P))) +
        labs(y= paste0("-log10(",p,")"),x=logFC,color = logFC, size = paste0("-log10(",p,")")) + 
        geom_point(data =  deg %>%
                     tidyr::drop_na() %>%
                     dplyr::filter(change != "Stable") %>%
                     dplyr::arrange(desc(FC)) %>%
                     dplyr::slice(c(1:(top/2),(dim(.)[1]-(top/2-1)):dim(.)[1])),
                   aes(x = FC, y = -log10(P),
                       size = -log10(P)),
                   shape = 21, show.legend = F, color = "#000000")+
        geom_text_repel(data =  Up1,
                        aes(x = FC, y = -log10(P), label = Symbol),
                        box.padding = 0.5,
                        nudge_y = 0.2,
                        #                     hjust = 0,
                        #                    segment.size = 0.5,
                        #                    segment.linetype = 1,
                        #                    max.overlaps = 100,
                        #                    max.iter = 1000000,
                        #                    max.time = 10,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        nudge_x =  Up1[,logFC]+(max_lfc/2),
                        min.segment.length = 0,
                        fontface = "bold",
                        family = "Times",
                        direction = "y", 
                        hjust = "left"
        ) + 
        geom_text_repel(data =  Down1,
                        aes(x = FC, y = -log10(P), label =Symbol),
                        box.padding = 0.5,
                        nudge_x =  Down1[,logFC]+(min_lfc/2),
                        min.segment.length = 0,
                        fontface = "bold",
                        family = "Times",
                        nudge_y = 0.2,
                        #                  hjust = 0,
                        #                 segment.size = 0.5,
                        #                 segment.linetype = 1,
                        #                 max.overlaps = 100,
                        #                 max.iter = 1000000,
                        #                 max.time = 10,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        direction = "y", 
                        hjust = "right"           
        )+ 
        scale_color_gradientn(colours = c("#00ffef", "#00cc99","#ffffbf", "#f8d568", "#fb9902"),
                              values = seq(0, 1, 0.2)) +
        scale_fill_gradientn(colours = c("#00ffef", "#00cc99","#ffffbf", "#f8d568", "#fb9902"),
                             values = seq(0, 1, 0.2)) +
        geom_vline(xintercept = c(-cut_FC, cut_FC), linetype = 2) +
        geom_hline(yintercept = -log10(cut_P), linetype = 4) + 
        scale_size(range = c(1,7)) + 
        xlim(c(-x.num, x.num)) + 
        ylim(c(0, y.num+1)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(),
              legend.background = element_roundrect(color = "#808080", linetype = 1),
              axis.text = element_text(size = 13, color = "#000000"),
              axis.title = element_text(size = 15),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)
        ) + 
        #   annotate(geom = "text", x = -x.num+log2(x.num), y =-log10(0.05), label = paste0("p = ",cut_P), size =6) + 
        coord_cartesian(clip = "off") + 
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
            gp = grid::gpar(lwd = 3, col = "#74add1")
          ), 
          xmin = -x.num, 
          xmax = -1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label = paste0("Down ( N=",table(deg$change)[1],")"),
            gp = grid::gpar(col = "#74add1")
          ),
          xmin = -x.num, 
          xmax = -1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
            gp = grid::gpar(lwd = 3, col = "#d73027")
          ), 
          xmin = x.num, 
          xmax = 1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label =  paste0("Up ( N=",table(deg$change)[3],")"),
            gp = grid::gpar(col = "#d73027")
          ),
          xmin = x.num, 
          xmax = 1,
          ymin = y.num+1,
          ymax = y.num+1
        )
    }else{
      Up2<- deg %>% 
        tidyr::drop_na() %>%
        dplyr::filter(change != "Stable") %>%
        dplyr::arrange(desc(-log10(P))) %>%
        dplyr::slice(1:top) %>%
        dplyr::filter(change == "Up")
      Down2<-deg %>% 
        tidyr::drop_na() %>%
        dplyr::filter(change != "Stable") %>%
        dplyr::arrange(desc(-log10(P))) %>%
        dplyr::slice(1:top) %>%
        dplyr::filter(change == "Down")
      p1<-ggplot(data = deg) + 
        geom_point(aes(x = FC, y = -log10(P), 
                       color = FC,
                       size = -log10(P))) +
        labs(y= paste0("-log10(",p,")"),x=logFC,color = logFC, size = paste0("-log10(",p,")")) + 
        geom_point(data =  deg %>%
                     tidyr::drop_na() %>%
                     dplyr::filter(change != "Stable") %>%
                     dplyr::arrange(desc(-log10(P))) %>%
                     dplyr::slice(1:top),
                   aes(x = FC, y = -log10(P),
                       size = -log10(P)),
                   shape = 21, show.legend = F, color = "#000000")+
        geom_text_repel(data =Up2,
                        aes(x = FC, y = -log10(P), label = Symbol),
                        box.padding = 0.5,
                        nudge_x =  Up2[,logFC]+(max_lfc/10),
                        #                   nudge_y = 0.5,
                        nudge_y = 0.5,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        direction = "y",
                        hjust = "left"
        ) + 
        geom_text_repel(data =Down2  ,
                        aes(x = FC, y = -log10(P), label =Symbol),
                        nudge_x =  Down2[,logFC]-(min_lfc/10),
                        box.padding = 0.5,
                        nudge_y = 0.5,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        direction = "y", 
                        hjust = "right"
        )+scale_color_gradientn(colours = c("#00ffef", "#00cc99","#ffffbf", "#f8d568", "#fb9902"),
                                values = seq(0, 1, 0.2)) +
        scale_fill_gradientn(colours = c("#00ffef", "#00cc99","#ffffbf", "#f8d568", "#fb9902"),
                             values = seq(0, 1, 0.2)) +
        
        geom_vline(xintercept = c(-cut_FC, cut_FC), linetype = 2) +
        geom_hline(yintercept = -log10(cut_P), linetype = 4) + 
        scale_size(range = c(1,7)) + 
        xlim(c(-x.num, x.num)) + 
        ylim(c(0, y.num+1)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(),
              legend.background = element_roundrect(color = "#808080", linetype = 1),
              axis.text = element_text(size = 13, color = "#000000"),
              axis.title = element_text(size = 15),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)
        ) + 
        annotate(geom = "text", x = -x.num+log2(x.num), y =-log10(0.05), label = paste0("p = ",cut_P), size =6) + 
        coord_cartesian(clip = "off") + 
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
            gp = grid::gpar(lwd = 3, col = "#74add1")
          ), 
          xmin = -x.num, 
          xmax = -1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label = paste0("Down ( N=",table(deg$change)[1],")"),
            gp = grid::gpar(col = "#74add1")
          ),
          xmin = -x.num, 
          xmax = -1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
            gp = grid::gpar(lwd = 3, col = "#d73027")
          ), 
          xmin = x.num, 
          xmax = 1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label =  paste0("Up ( N=",table(deg$change)[3],")"),
            gp = grid::gpar(col = "#d73027")
          ),
          xmin = x.num, 
          xmax = 1,
          ymin = y.num+1,
          ymax = y.num+1
        )
    }
    if(is.null(highlight)){return(p1)
    }else{
      highlight1 <- data.frame(Symbol = highlight)
      hlight_data <- merge(highlight1,deg)
      p2 <- p1 +
        geom_point(size = 4,data = hlight_data,aes(x=FC,y=-log10(P)),shape = 21,stroke = 0.53,
                   color = "green") +
        geom_text_repel(
          data = hlight_data[which(hlight_data[,logFC] > 0),],
          aes(x=FC,y=-log10(P),label = Symbol), color="#E43535",
          nudge_y      =0.5,
          size=5,
          direction    = "y",
          hjust        = 0,
          segment.size = 0.5,
          segment.linetype = 1,
          max.overlaps = 100,
          max.iter = 1000000,
          max.time = 10,
          nudge_x =  hlight_data[which(hlight_data[,logFC] > 0),][,logFC]+max_lfc,
          min.segment.length = 0,
          fontface = "bold",family = "Times")+
        geom_text_repel(
          data = hlight_data[which(hlight_data[,logFC] < 0),],
          aes(x=FC,y=-log10(P),label = Symbol), color="#E43535",
          nudge_y      =- 0.5,
          size=5,
          direction    = "y",
          hjust        = 0,
          segment.size = 0.5,
          segment.linetype = 1,
          max.overlaps = 100,
          max.iter = 1000000,
          max.time = 10,
          nudge_x =  hlight_data[which(hlight_data[,logFC] < 0),][,logFC]+min_lfc,
          min.segment.length = 0,
          fontface = "bold",family = "Times"
        )
      return(p2)
    }
  }
  else{
    if(orderByFC==TRUE){
      Up1<-deg %>% 
        tidyr::drop_na() %>%
        dplyr::filter(change != "Stable") %>%
        dplyr::arrange(desc(FC)) %>%
        dplyr::slice(1:(top/2)) %>%
        dplyr::filter(change == "Up")    
      Down1<-deg %>% 
        tidyr::drop_na() %>%
        dplyr::filter(change != "Stable") %>%
        dplyr::arrange(desc(FC)) %>%
        dplyr::slice((dim(.)[1]-(top/2-1)):dim(.)[1]) %>%
        dplyr::filter(change == "Down")
      
      p1<-ggplot(data = deg) + 
        geom_point(aes(x = FC, y = -log10(P), 
                       color = FC,
                       size = -log10(P))) +
        labs(y= paste0("-log10(",p,")"),x=logFC,color = logFC, size = paste0("-log10(",p,")")) + 
        geom_point(data =  deg %>%
                     tidyr::drop_na() %>%
                     dplyr::filter(change != "Stable") %>%
                     dplyr::arrange(desc(FC)) %>%
                     dplyr::slice(c(1:(top/2),(dim(.)[1]-(top/2-1)):dim(.)[1])),
                   aes(x = FC, y = -log10(P),
                       size = -log10(P)),
                   shape = 21, show.legend = F, color = "#000000")+
        geom_text_repel(data =  Up1,
                        aes(x = FC, y = -log10(P), label = Symbol),
                        box.padding = 0.5,
                        nudge_y = 0.2,
                        #                     hjust = 0,
                        #                    segment.size = 0.5,
                        #                    segment.linetype = 1,
                        #                    max.overlaps = 100,
                        #                    max.iter = 1000000,
                        #                    max.time = 10,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        nudge_x =  Up1[,logFC]+(max_lfc/2),
                        min.segment.length = 0,
                        fontface = "bold",
                        family = "Times",
                        direction = "y", 
                        hjust = "left"
        ) + 
        geom_text_repel(data =  Down1,
                        aes(x = FC, y = -log10(P), label =Symbol),
                        box.padding = 0.5,
                        nudge_x =  Down1[,logFC]+(min_lfc/2),
                        min.segment.length = 0,
                        fontface = "bold",
                        family = "Times",
                        nudge_y = 0.2,
                        #                  hjust = 0,
                        #                 segment.size = 0.5,
                        #                 segment.linetype = 1,
                        #                 max.overlaps = 100,
                        #                 max.iter = 1000000,
                        #                 max.time = 10,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        direction = "y", 
                        hjust = "right"           
        )+
        scale_color_gradientn(colours = c("#00ffef", "#00cc99","#ffffbf", "#f8d568", "#fb9902"),
                              values = seq(0, 1, 0.2)) +
        scale_fill_gradientn(colours = c("#00ffef", "#00cc99","#ffffbf", "#f8d568", "#fb9902"),
                             values = seq(0, 1, 0.2)) +
        
        geom_vline(xintercept = c(-cut_FC, cut_FC), linetype = 2) +
        geom_hline(yintercept = -log10(cut_P), linetype = 4) + 
        scale_size(range = c(1,7)) + 
        xlim(c(-x.num, x.num)) + 
        ylim(c(0, y.num+1)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(),
              legend.background = element_roundrect(color = "#808080", linetype = 1),
              axis.text = element_text(size = 13, color = "#000000"),
              axis.title = element_text(size = 15),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)
        ) + 
        #   annotate(geom = "text", x = -x.num+log2(x.num), y =-log10(0.05), label = paste0("p = ",cut_P), size =6) + 
        coord_cartesian(clip = "off") + 
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
            gp = grid::gpar(lwd = 3, col = "#74add1")
          ), 
          xmin = -x.num, 
          xmax = -1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label = paste0("Down ( N=",table(deg$change)[1],")"),
            gp = grid::gpar(col = "#74add1")
          ),
          xmin = -x.num, 
          xmax = -1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
            gp = grid::gpar(lwd = 3, col = "#d73027")
          ), 
          xmin = x.num, 
          xmax = 1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label =  paste0("Up ( N=",table(deg$change)[3],")"),
            gp = grid::gpar(col = "#d73027")
          ),
          xmin = x.num, 
          xmax = 1,
          ymin = y.num+1,
          ymax = y.num+1
        )+
        annotate("rect", xmin = cut_FC, xmax = x.num, ymin =-log10(cut_P), ymax = y.num, alpha = 0.2,fill=background[2]) +
        annotate("rect", xmin = -cut_FC, xmax = -x.num, ymin = -log10(cut_P), ymax = y.num, alpha = 0.2,fill=background[1]) 
    }else{
      Up2<- deg %>% 
        tidyr::drop_na() %>%
        dplyr::filter(change != "Stable") %>%
        dplyr::arrange(desc(-log10(P))) %>%
        dplyr::slice(1:top) %>%
        dplyr::filter(change == "Up")
      Down2<-deg %>% 
        tidyr::drop_na() %>%
        dplyr::filter(change != "Stable") %>%
        dplyr::arrange(desc(-log10(P))) %>%
        dplyr::slice(1:top) %>%
        dplyr::filter(change == "Down")
      p1<-ggplot(data = deg) + 
        geom_point(aes(x = FC, y = -log10(P), 
                       color = FC,
                       size = -log10(P))) +
        labs(y= paste0("-log10(",p,")"),x=logFC,color = logFC, size = paste0("-log10(",p,")")) + 
        geom_point(data =  deg %>%
                     tidyr::drop_na() %>%
                     dplyr::filter(change != "Stable") %>%
                     dplyr::arrange(desc(-log10(P))) %>%
                     dplyr::slice(1:top),
                   aes(x = FC, y = -log10(P),
                       size = -log10(P)),
                   shape = 21, show.legend = F, color = "#000000")+
        geom_text_repel(data =Up2,
                        aes(x = FC, y = -log10(P), label = Symbol),
                        box.padding = 0.5,
                        nudge_x =  Up2[,logFC]+(max_lfc/10),
                        #                   nudge_y = 0.5,
                        nudge_y = 0.5,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        direction = "y",
                        hjust = "left"
        ) + 
        geom_text_repel(data =Down2  ,
                        aes(x = FC, y = -log10(P), label =Symbol),
                        nudge_x =  Down2[,logFC]-(min_lfc/10),
                        box.padding = 0.5,
                        nudge_y = 0.5,
                        segment.curvature = -0.1,
                        segment.ncp = 3,
                        segment.angle = 20,
                        direction = "y", 
                        hjust = "right"
        )+
        scale_color_gradientn(colours = c("#00ffef", "#00cc99","#ffffbf", "#f8d568", "#fb9902"),
                              values = seq(0, 1, 0.2)) +
        scale_fill_gradientn(colours = c("#00ffef", "#00cc99","#ffffbf", "#f8d568", "#fb9902"),
                             values = seq(0, 1, 0.2)) +
        
        geom_vline(xintercept = c(-cut_FC, cut_FC), linetype = 2) +
        geom_hline(yintercept = -log10(cut_P), linetype = 4) + 
        scale_size(range = c(1,7)) + 
        xlim(c(-x.num, x.num)) + 
        ylim(c(0, y.num+1)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(),
              legend.background = element_roundrect(color = "#808080", linetype = 1),
              axis.text = element_text(size = 13, color = "#000000"),
              axis.title = element_text(size = 15),
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5)
        ) + 
        annotate(geom = "text", x = -x.num+log2(x.num), y =-log10(0.05), label = paste0("p = ",cut_P), size =6) + 
        coord_cartesian(clip = "off") + 
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
            gp = grid::gpar(lwd = 3, col = "#74add1")
          ), 
          xmin = -x.num, 
          xmax = -1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label = paste0("Down ( N=",table(deg$change)[1],")"),
            gp = grid::gpar(col = "#74add1")
          ),
          xmin = -x.num, 
          xmax = -1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::segmentsGrob(
            y0 = unit(-10, "pt"),
            y1 = unit(-10, "pt"),
            arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
            gp = grid::gpar(lwd = 3, col = "#d73027")
          ), 
          xmin = x.num, 
          xmax = 1,
          ymin = y.num+1,
          ymax = y.num+1
        ) +
        annotation_custom(
          grob = grid::textGrob(
            label =  paste0("Up ( N=",table(deg$change)[3],")"),
            gp = grid::gpar(col = "#d73027")
          ),
          xmin = x.num, 
          xmax = 1,
          ymin = y.num+1,
          ymax = y.num+1
        )+
        annotate("rect", xmin = cut_FC, xmax = x.num, ymin =-log10(cut_P), ymax = y.num, alpha = 0.2,fill=background[2]) +
        annotate("rect", xmin = -cut_FC, xmax = -x.num, ymin = -log10(cut_P), ymax = y.num, alpha = 0.2,fill=background[1]) 
    }
    if(is.null(highlight)){return(p1)
    }else{
      highlight1 <- data.frame(Symbol = highlight)
      hlight_data <- merge(highlight1,deg)
      p2 <- p1 +
        geom_point(size = 4,data = hlight_data,aes(x=FC,y=-log10(P)),shape = 21,stroke = 0.53,
                   color = "green") +
        geom_text_repel(
          data = hlight_data[which(hlight_data[,logFC] > 0),],
          aes(x=FC,y=-log10(P),label = Symbol), color="#E43535",
          nudge_y      =0.5,
          size=5,
          direction    = "y",
          hjust        = 0,
          segment.size = 0.5,
          segment.linetype = 1,
          max.overlaps = 100,
          max.iter = 1000000,
          max.time = 10,
          nudge_x =  hlight_data[which(hlight_data[,logFC] > 0),][,logFC]+max_lfc,
          min.segment.length = 0,
          fontface = "bold",family = "Times")+
        geom_text_repel(
          data = hlight_data[which(hlight_data[,logFC] < 0),],
          aes(x=FC,y=-log10(P),label = Symbol), color="#E43535",
          nudge_y      =- 0.5,
          size=5,
          direction    = "y",
          hjust        = 0,
          segment.size = 0.5,
          segment.linetype = 1,
          max.overlaps = 100,
          max.iter = 1000000,
          max.time = 10,
          nudge_x =  hlight_data[which(hlight_data[,logFC] < 0),][,logFC]+min_lfc,
          min.segment.length = 0,
          fontface = "bold",family = "Times"
        )+
        annotate("rect", xmin = cut_FC, xmax = x.num, ymin =-log10(cut_P), ymax = y.num, alpha = 0.2,fill=background[2]) +
        annotate("rect", xmin = -cut_FC, xmax = -x.num, ymin = -log10(cut_P), ymax = y.num, alpha = 0.2,fill=background[1]) 
      return(p2)
    }
  }
}

# p <- Volcanol_1(na.omit(deg),
#                 logFC="log2FoldChange",
#                 p="padj",
#                 cut_P=0.05,
#                 cut_FC=1,
#                 x.num=0,
#                 y.num=0,
#                 top=10,
#                 highlight=NULL,
#                 orderByFC=TRUE,
#                 background=c("#74add1","#d73027"))
# p