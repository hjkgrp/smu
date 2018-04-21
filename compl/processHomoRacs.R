#rm(list=ls())

# read in the homoRACs
racs <- read.csv(file = "homolepticRACs.csv", header = TRUE)
load(file="selected_all.R")
load(file="selected_rf_41.Rda")
# remove alpha and ox from rf41 because we only want topology related racs for the fingerprint.
rf39vars <- selected_rf_41[!(selected_rf_41 %in% c('alpha','ox'))]
rf153vars <- selected_all[!(selected_all %in% c('alpha','ox'))]

rf39racs <- racs[,rf39vars]
rf153racs <- racs[,rf153vars]

# apply the var on cols (meaning of '2') of rf39racs. which ones have it zero?
# which(apply(rf153racs, 2, var)==0)
rf133racs <- rf153racs[ , apply(rf153racs, 2, var) != 0]


plot(rf133racs[,c(1:6)], pch='.')

############################
# PCA with function prcomp #
############################
pca1 = prcomp(rf133racs,scale. = TRUE) #scale. = TRUE

# sqrt of eigenvalues
eigendecay = pca1$sdev
plot(eigendecay)

# loadings
head(pca1$rotation)

# PCs (aka scores) [PC1, PC2, ... data points]
scores = as.data.frame(pca1$x)

require(plyr)
require(dplyr)

# mc.Z.0.all is the metal charge and since autocorrelation goes Z1*Z2, we sqrt it
scores$metal <- revalue(as.factor(sqrt(rf133racs$mc.Z.0.all)),c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
# the sqrt of Z1Z2 of a ligand centric first atom is the CA
scores$cai <- revalue(as.factor(sqrt(rf133racs$lc.Z.0.eq)),c('6'='C','7'='N','8'='O','15'='P','16'='S'))
# electronegativity 
scores$enm <- rf133racs$D_mc.chi.2.all
# proxy for number of atoms
scores$nats <- rf133racs$f.chi.0.all




# graphics.off()
# x11()
# PC1 vs PC5 separates the metals
# PC4 vs PC2 is interesting: cai and esp f.Z.0.ax
g<-ggplot(data=scores,aes(x=PC4,y=PC2,color=rf133racs$f.Z.0.ax, shape=cai)) +geom_point(size=2) +theme_light() 
  #scale_color_gradient2(high='firebrick',low='dodgerblue',mid='gray',midpoint=-15)
print(g)

ggplot(data = scores, aes(x = PC2, y = PC4, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "darkcyan", alpha = 0.8, size = 2)


plot(scores[1:6], pch='.')

