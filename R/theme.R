library(ggplot2)

mainTheme <- theme(panel.grid.major = element_blank() ,
                   text = element_text(size=13) ,
                   axis.title.y = element_text(margin = margin(t = 0, r = 18, b = 0, l = 0)) ,
                   axis.title.x = element_text(margin = margin(t = 18, r = 0, b = 0, l = 0)) ,
                   legend.title = element_blank() ,
                   legend.margin=margin(c(0.3,1,0.3,1), unit='lines') ,
                   legend.background = element_rect(fill="white", linetype="solid",  colour ="#979797"))

theme_map <- theme(plot.background = element_rect(fill = "#FFFFFF", color = NA),
                   panel.background = element_rect(fill = "#F8F8F8", color = NA),
                   panel.border = element_rect(colour = "black", fill  =NA, linewidth = 0.1),
                   text = element_text(size=8),
                   axis.text.x =  element_text(size=6),
                   axis.text.y =  element_text(size=6))


theme_plot <- theme(plot.background = element_rect(fill = "#FFFFFF", color = NA),
                    panel.background = element_rect(fill = "#F8F8F8", color = NA),
                    panel.border = element_rect(colour = "black", fill  =NA, linewidth = 0.1),
                    panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    text = element_text(size=8),
                    axis.text.x =  element_text(size=6),
                    axis.text.y =  element_text(size=6))


myColors <- c("#6FBBE8","#A1ECD8","#F6F9AB","#FCB46D","#B21414")