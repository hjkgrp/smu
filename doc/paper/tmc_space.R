library("ggrepel")

# SMALL: PCA of different sets into SU with facet = set
g <- ggplot(data=projHomoIntoSu2[seq.int(1L,length(projHomoIntoSu2$PC1),1L),],
            aes(x=PC1, y=PC2)) +
  geom_point(data=scores2[seq.int(1L,length(scores2$PC1),10L),], size = 0.1, color = 'black',shape=1) +
  geom_point(size = 0.15, color='red',shape=0) +
  ggtitle("PCA of low (black) and high (red) symmetry") + theme(plot.margin =margin(0.05,0.05,0.05,0.05,"cm")) +
  theme_light()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        plot.title = element_text(size=12),
        axis.title = element_text(size=12)
        )+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.25,color='black'),
        panel.background = element_rect(fill='gray98'))
  
  

#--------------------

dat <- NULL
dat$x <- seq(1,4)
dat$y <- c(2,4,6,8)
dat$label <- c('Homoleptic', 'Strong Symmetry', 'Weak Symmetry', 'Complete Heteroleptics' )
dat <- data.frame(dat)

# big
h <- ggplot(dat, aes(x, y)) + geom_point(shape = 1, size = 2, stroke = 1) + theme_light() + ylim(0,10.5) + xlim(1, 8)  +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1.25,color='black')) +
  xlab('Resultion') +
  ylab('log(space size)') +
  geom_text(aes(label = label),
            color = "gray20",
            data = dat,
            position = position_nudge(y = 0.4)
  ) +
  geom_segment(aes(x = 2.63, y = 4.55, xend = 4.2, yend = 6.9), color="gray", linetype="dashed", size=0.5) +
  geom_segment(aes(x = 2.63, y = 4.15, xend = 4.2, yend = 0.1  ), color="gray", linetype="dashed", size=0.5)+
  ggtitle("Size of TMC space as a function of resultion")


  

vp <- viewport(width=0.47, height =  0.6, x = 0.72, y = 0.36)
print(h)
print(g,vp=vp)

