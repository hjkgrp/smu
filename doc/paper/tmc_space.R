t1 <- textGrob(expression("Concentration of " * phantom(bold("affluence")) * "and" * phantom(bold("poverty")) * " nationwide"),
               x = 0.5, y = 1.1, gp = gpar(col = "black"))

t2 <- textGrob(expression(phantom("Concentration of ") * bold("affluence") * phantom(" and poverty nationwide")),
               x = 0.5, y = 1.1, gp = gpar(col = "#EEB422"))

t3 <- textGrob(expression(phantom("Concentration of affluence and ") * bold("poverty") * phantom(" nationwide")),
               x = 0.5, y = 1.1, gp = gpar(col = "#238E68"))


# PCA of different sets into SU with facet = set
g <- ggplot(data=projHomoIntoSu2[seq.int(1L,length(projHomoIntoSu2$PC1),1L),],
            aes(x=PC1, y=PC2)) +
  geom_point(size = 0.1, color='red')
g <- g + geom_point(data=scores2[seq.int(1L,length(scores2$PC1),10L),], size = 0.1, color = 'black')
# g <- g + guides(color = guide_legend(override.aes = list(size = 5)))
g <- g + ggtitle("PCA of" ,bold(low), "(black) and high (red) symmetry") +
labs(color='Convergence') + theme_light()+theme(plot.margin =margin(0.05,0.05,0.05,0.05,"cm")) +
  theme(axis.text.x=element_blank(),
axis.text.y=element_blank(),axis.ticks=element_blank())

dat <- NULL
dat$x <- seq(1,5)
dat$y <- c(2,4,6,8,10)
dat <- data.frame(dat)

h <- ggplot(dat, aes(x, y)) + geom_point() + theme_light() + ylim(0,10) + xlim(1, 8)  +
  theme(axis.title.x='ebi',
        axis.title.y='er',
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank())
  

vp <- viewport(width=0.47, height =  0.6, x = 0.72, y = 0.36)
print(h)
print(g,vp=vp)
