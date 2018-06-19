library("ggrepel")
library(png)
library(grid)
library(magick)
library(magrittr) # For piping the logo

theme_inlay <- function () { 
  theme(plot.margin =margin(0.0,0.0,0.0,0.0,"cm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        plot.title = element_text(size=9),
        axis.title = element_text(size=9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.4,color='black'),
        aspect.ratio = 1.0,
        panel.background = element_blank(), # bg of the panel
        plot.background = element_blank() # bg of the plot
  )
}

# SMALL: PCA of different sets into SU with facet = set
g <- ggplot(data=projHomoIntoSu2[seq.int(1L,length(projHomoIntoSu2$PC1),1L),],
            aes(x=PC1, y=PC2)) +
  geom_point(data=scores2[seq.int(1L,length(scores2$PC1),10L),], size = 0.1, color = 'black',shape=1) +
  geom_point(size = 0.1, color='red',shape=1) +
  # ggtitle("PCA of low (black) and high (red) symmetry")
  theme_light()+
  theme_inlay()+
  scale_x_discrete(position = "top")
  
  
  

#--------------------
# big
img <- readPNG('hom.png')
hom <- rasterGrob(img, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank") +
  annotation_custom(hom, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_point()

dat <- NULL
dat$x <- seq(1,4)*1.7+0.1
# dat$y <- c(3,5,7,9)*0.8+1
dat$y <- log(c(40,360,1140,4000000))
dat$label <- c('Homoleptic', 'Strong Symmetry', 'Weak Symmetry', 'Heteroleptics' )
dat$label2 <- c('40', '360', '1140', '4 %.% 10^6') 
dat <- data.frame(dat)

h <- ggplot(dat, aes(x, y)) + geom_point(shape = 1, size = 1.4, stroke = 0.4) +
  theme_light() +
  #ylim(0,11) +
  xlim(1, 8) +
  xlab('symmetry class') +
  ylab('log(space size)') +
  
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=0.5,color='black'),
        axis.title = element_text(size=10),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0))
        ) +
  
  geom_text(aes(label = label),
            color = "gray20",
            data = dat,
            position = position_nudge(y = 0.5),
            size = 3
  ) +
  geom_text(aes(label = label2),
            parse = TRUE,
            color = "gray20",
            data = dat,
            position = position_nudge(y = -0.5),
            size = 3
  )
  

  # geom_segment(aes(x = 2.93, y = 5.65, xend = 4.2, yend = 6.82), color="gray", linetype="dashed", size=0.5) +
  # geom_segment(aes(x = 2.93, y = 5.25, xend = 4.2, yend = 0.25  ), color="gray", linetype="dashed", size=0.5)
  # ggtitle("Size of TMC space as a function of resultion")

cairo_pdf(file="tmc_viz.pdf",width = 3.3, height = 2.95)
vp <- viewport(width = 0.5, height = 0.5, x = 0.75, y = 0.36)
print(h)
print(g,vp=vp)
dev.off()

vp <- viewport(width=0.5, height =  0.5, x = 0.75, y = 0.36)
vp2 <- viewport(width=0.5, height =  0.5, x = 0.1, y = 0.1)

print(h)
print(g,vp=vp)
print(hom, vp=vp2)
