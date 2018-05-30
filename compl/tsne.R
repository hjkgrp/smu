# make T-SNE
tsne_data <- tsne(split_[,1:39], k=2,perplexity = 40, max_iter=2500, epoch=150)
all_tsne <- tsne_data

# make DF from T-SNE data
plot_df <- data.frame(all_tsne)
colnames(plot_df) <- c('x','y')
# add in things to color by
plot_df$split = c(4,4,4,4,3,4,3,2,2,4,4,4,4,4,4,4,3,4,4,4,4,2,4,4,3,4,2,2,3,2,4,3,4,3,3,3,2,3,3,3,4,3,4,4,2,4,3,2,2,4,3,4,4,3,4,4,4,4,4,4,4,4,3,4,2,3,4,3,4,2,3,3,3,4,4,4,3,4,3,3,4,3,2,3,4,4,2,3)
# labels so we can check on it
plot_df$ID <- split_[,40]
# here comes the akima
require(akima)
pts.grid <- interp(as.data.frame(plot_df)$x, 
                   as.data.frame(plot_df)$y, 
                   as.data.frame(plot_df)$split,
                   duplicate = 'mean')
pts.grid2 <- expand.grid(x=pts.grid$x, y=pts.grid$y)
pts.grid2$z <- as.vector(pts.grid$z)

# plot 
graphics.off()
# x11()
p <- ggplot(aes(x=x,y= y),size=8, data = plot_df)
p <- p + geom_raster(data=na.omit(pts.grid2), aes(x=x, y=y,fill=z),size=0,alpha=1,interpolate=TRUE) 
p <- p + geom_point(aes(color=split),size=3) + scale_fill_gradient2(low="darkblue",high= "firebrick",mid="white",midpoint=0)
p <- p + theme_normal() + theme(legend.position = c(0.1,0.15)) +xlab('') + ylab('')
p <- p  + theme(axis.text.x = element_blank(),axis.text.y = element_blank())
print(p)
