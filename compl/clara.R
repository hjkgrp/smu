require(cluster)
library(cluster)
new_descriptors <- validChildren
clara_predictors <- new_descriptors
rownames(clara_predictors) <- clara_predictors$runs
# load(file="selected_all.R")
load(file="selected_rf_41.Rda")
rf39vars <- selected_rf_41[!(selected_rf_41 %in% c('alpha','ox'))] # bc we dont have alpha and ox in the racs here
clara_predictors$alpha<-0.2
clara_predictors <- clara_predictors[,rf39vars]
clara_predictors$alpha<-NULL

clara_predictors <- scale(clara_predictors,center=TRUE,scale=TRUE)
clara_predictors_center <- attr(clara_predictors, 'scaled:center')
clara_predictors_scale <- attr(clara_predictors, 'scaled:scale')
clara_predictors <- as.data.frame(clara_predictors)


nc = 200
start_time <- Sys.time()
pa <-clara(clara_predictors,nc,rngR = TRUE) # rngR uses the R random number generator
end_time <- Sys.time()
# pa200 <- pa
# save(pa200,file="pa200.Rdata")

cluster_info <- data.frame(pa$clusinfo)
meds  <- as.data.frame(t(apply(pa$medoids, 1, function(r)r*clara_predictors_scale+ clara_predictors_center)))

library(plyr)
library(dplyr)

meds$metal=factor(sqrt(meds$mc.Z.0.all))
meds$ox=factor(meds$ox)
meds$metal=revalue(meds$metal,c("24"='Cr',"25"="Mn","26"="Fe","27"="Co"))



# Isolation Plot
graphics.off()
x11()
prop_plot <- ggplot(aes(x = size, color = max_diss ,y = isolation), size = 8, data = cluster_info)
prop_plot <- prop_plot + geom_point(size = 4) + theme_light() + guides(color = guide_colorbar(title = "max diss."))
prop_plot <- prop_plot + theme(legend.title = element_text(), legend.position = c(0.8,0.25))
prop_plot <- prop_plot + scale_color_continuous(low="gray70", high="blue")

print(prop_plot)

cairo_pdf('clustal1.pdf',width=2*3.33,height = 4.00)
print(prop_plot)
dev.off()

# Metal distribution
x11()
prop_plot <- ggplot(aes(x = metal),size=8, data = meds)
# prop_plot <- ggplot(aes(x = metal, fill=ox, group=ox),size=8, data = meds)
prop_plot <- prop_plot +  geom_bar(position = position_dodge()) + theme_light()
print(prop_plot)

graphics.off()
cairo_pdf('clustal2.pdf',width=2*3.33,height = 4.00)
print(prop_plot)
dev.off()


## make PCA for comparison
PCA = prcomp(clara_predictors, scale. = F, center = F,tol = sqrt(.Machine$double.eps))

dfPCA <- as.data.frame(PCA$x)
dfPCA$cluster <- (pa$clustering)

graphics.off()
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
na_clusc <- scale_color_gradientn(colours = myPalette(100), limits=c(0,nc),guide = guide_colorbar(title='cluster'))
na_clusd <- scale_color_discrete(guide = FALSE)
naf_clusc <- scale_fill_gradientn(colours = myPalette(100), limits=c(0,nc),guide = guide_colorbar(title='cluster'))

# PCA of first two components
find_hull <- function(df) df[chull(df$PC1, df$PC2), ]
hull<-find_hull(dfPCA[rownames(meds),])


x11()
prop_plot <- ggplot(aes(x = PC1,y = PC2, color = cluster),size = 2, data = dfPCA)
prop_plot <- prop_plot +  geom_point(shape = 16, size = 3, alpha = 0.25) + na_clusc
prop_plot <- prop_plot +  geom_point(aes(fill = cluster), color = 'black', shape = 24, size = 4, stroke = 2,
                                     alpha = 1, data = dfPCA[rownames(dfPCA) %in% rownames(pa$medoids),])+
                                     na_clusc + naf_clusc
prop_plot <- prop_plot + geom_polygon(data = hull, alpha = 1,fill=NA,size=1.5)+geom_point(data=hull,size=8,color='black') #+
  #geom_point(data=dfPCA[rownames(meds[2,]),],size=8,color='white')
prop_plot <- prop_plot + xlab('PC 1') + ylab('PC 2') + 
  theme(text = element_text(size=20,family = 'sans'),
    panel.background = element_rect(fill='white'), axis.text.x = element_text(face = "bold", color = "black"),
    axis.title.x=element_text(face = "bold", color = "black"), axis.line = element_line(colour = "black"),
    axis.text.y=element_text(face = "bold", color = "black"), 
    axis.title.y=element_text(face = "bold", color = "black"),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(prop_plot)

cairo_pdf('clustv1_medoid_hull.pdf',width=2*3.5,height = 4.00)
print(prop_plot)
dev.off()


# PCA of second two components
x11()
prop_plot <- ggplot(aes(x=PC3,y= PC4,color=cluster),size=8, data = dfPCA)
prop_plot <- prop_plot +  geom_point(shape=16,size=5,alpha=0.25)  + na_clusc
prop_plot <- prop_plot +  geom_point(aes(fill=cluster),color='black',shape=24,size=4,stroke=2,
                                     alpha=1,data=dfPCA[rownames(dfPCA) %in% rownames(pa$medoids),])+
  na_clusc +naf_clusc

# add the convex hull from pc1 and pc2 and project it here
#prop_plot <- prop_plot + geom_polygon(data = hull, alpha = 1,fill=NA,size=1.5)+geom_point(data=hull,size=8,color='black')+
  geom_point(data=dfPCA[rownames(meds[2,]),],size=8,color='white')
prop_plot <- prop_plot +     xlab('PC 3') + ylab('PC 4') +    theme(text = element_text(size=20,family = 'sans'),panel.background = element_rect(fill='white'),axis.text.x=element_text(face = "bold", color = "black"),
                                                                    axis.title.x=element_text(face = "bold", color = "black"), axis.line = element_line(colour = "black"),
                                                                    axis.text.y=element_text(face = "bold", color = "black"), 
                                                                    axis.title.y=element_text(face = "bold", color = "black"),
                                                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank())

print(prop_plot)

cairo_pdf('clustv2.pdf',width=2*3.5,height = 4.00)
print(prop_plot)
dev.off()



#grid.draw(ggarrange(plots=list(fgp, fgp1)))



#write.csv(rownames(meds),file="distance_clara_100.csv")