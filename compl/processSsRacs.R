rm(list=ls())

# read in the RACs
ss_racs <- read.csv(file = "ssRACs.csv", header = TRUE)
fo_racs <- read.csv(file = "foRACs.csv", header = TRUE)
homo_racs <- read.csv(file = "homolepticRACs.csv", header = TRUE)

load(file="selected_all.R")
load(file="selected_rf_41.Rda")

# remove alpha and ox from rf41 because we only want topology related racs for the fingerprint.
rf39vars <- selected_rf_41[!(selected_rf_41 %in% c('alpha','ox'))]
rf153vars <- selected_all[!(selected_all %in% c('alpha','ox'))]

ss_rf39racs <- ss_racs[,rf39vars]
ss_rf153racs <- ss_racs[,rf153vars]

fo_rf39racs <- fo_racs[,rf39vars]
fo_rf153racs <- fo_racs[,rf153vars]

homo_rf39racs <- homo_racs[,rf39vars]
homo_rf153racs <- homo_racs[,rf153vars]

split_results_39racs <- split_results[,rf39vars]
split_results_153racs <- split_results[,rf153vars]
props <- (split_results[,!(colnames(split_results) %in% (rf153vars))])

# Apply the variance function 'var' on all cols ('2' means cols) of rf39racs. Throw out the one's with zero variance for PCA later.
# which(apply(rf153racs, 2, var)==0)
# apply the same deletion in the ss and fo even though they have actually more, to get same dimensionality
homo_rf133racs <- homo_rf153racs[ , apply(homo_rf153racs, 2, var) != 0]
ss_rf133racs   <- ss_rf153racs[ , apply(homo_rf153racs, 2, var) != 0]
fo_rf133racs   <- fo_rf153racs[ , apply(homo_rf153racs, 2, var) != 0]
split_results_133racs <- split_results_153racs[ , apply(homo_rf153racs, 2, var) != 0]

props$goodConvergence <- as.factor(props$ox_2_HS_flag_oct + props$ox_2_LS_flag_oct + props$ox_3_LS_flag_oct + props$ox_3_HS_flag_oct)

###############
# Standardize #
###############

# find moments of standardization
homo_rf133racs_data <- scale(homo_rf133racs, center=TRUE, scale=TRUE)
homo_rf133racs_center <-attr(homo_rf133racs_data, 'scaled:center')
homo_rf133racs_scale <- attr(homo_rf133racs_data, 'scaled:scale')
homo_rf133racs_data <- as.data.frame(homo_rf133racs_data)

# standardize ss and fo data with homo moments
ss_rf133racs_data <- scale(ss_rf133racs, center=homo_rf133racs_center, scale=homo_rf133racs_scale)
ss_rf133racs_data <- as.data.frame(ss_rf133racs_data)
fo_rf133racs_data <- scale(fo_rf133racs, center=homo_rf133racs_center, scale=homo_rf133racs_scale)
fo_rf133racs_data <- as.data.frame(fo_rf133racs_data)

# standardize results with homo moments
split_results_133racs_data <- scale(split_results_133racs, center=homo_rf133racs_center, scale=homo_rf133racs_scale)
split_results_133racs_data <- as.data.frame(split_results_133racs_data)

############################
# PCA with function prcomp #
############################
pca1 = prcomp(homo_rf133racs_data,scale. = FALSE) #scale. = TRUE

# sqrt of eigenvalues
eigendecay = pca1$sdev
plot(eigendecay)

# PCs (aka scores) [PC1, PC2, ... data points]
scores = as.data.frame(pca1$x)

scores$set <- 'homo'
# mc.Z.0.all is the metal charge and since autocorrelation goes Z1*Z2, we sqrt it
scores$metal <- revalue(as.factor(sqrt(homo_rf133racs$mc.Z.0.all)),c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
# the sqrt of Z1Z2 of a ligand centric first atom is the CA. below we have a cont spectr of atoms bc avg in ss is an atom in R
#scores$cai <- revalue(as.factor(sqrt(homo_rf133racs$lc.Z.0.eq)),c('6'='C','7'='N','8'='O','15'='P','16'='S'))
scores$cai <- (4*sqrt(homo_rf133racs$lc.Z.0.eq)+2*sqrt(homo_rf133racs$lc.Z.0.ax))/6

# projection of ss and fo into homo coordinates and reassignments
projSsIntoHomo <- as.data.frame(predict(pca1, ss_rf133racs_data))
projSsIntoHomo$set <- 'ss'
projSsIntoHomo$metal <- revalue(as.factor(sqrt(ss_rf133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projSsIntoHomo$cai <- (4*sqrt(ss_rf133racs$lc.Z.0.eq) + 2*sqrt(ss_rf133racs$lc.Z.0.ax))/6

projFoIntoHomo <- as.data.frame(predict(pca1, fo_rf133racs_data))
projFoIntoHomo$set <- 'fo'
projFoIntoHomo$metal <- revalue(as.factor(sqrt(fo_rf133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projFoIntoHomo$cai <- (4*sqrt(fo_rf133racs$lc.Z.0.eq) + 2*sqrt(fo_rf133racs$lc.Z.0.ax))/6

projSplitIntoHomo <- as.data.frame(predict(pca1, split_results_133racs_data))
projSplitIntoHomo$set <- 'calcHomo'
projSplitIntoHomo$metal <- revalue(as.factor(sqrt(split_results_133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projSplitIntoHomo$cai <- (4*sqrt(split_results_133racs$lc.Z.0.eq) + 2*sqrt(split_results_133racs$lc.Z.0.ax))/6

totalScores <- rbind(projSsIntoHomo, scores)

require(plyr)
require(dplyr)

# # electronegativity 
# totalScores$enm <- homo_rf133racs$D_mc.chi.2.all
# # proxy for number of atoms
# totalScores$nats <- homo_rf133racs$f.chi.0.all

# SS:
# PC1 vs PC2: metals not separated but cai sep in a patch
# f.Z.0.ax
# 4v6, 3v5, 1v2 show patches
# 7v8 cluster nicely
# HOMO: 
# PC1 vs PC5 separates the metals
# PC4 vs PC2 is interesting: cai and esp f.Z.0.ax
g <- ggplot(data=projSsIntoHomo[seq.int(1L,length(projSsIntoHomo$PC1),100L),], aes(x=PC1, y=PC2, color='black')) +
  geom_point(size = 0.1)
#g <- g + geom_point(data=projFoIntoHomo[seq.int(1L,length(projFoIntoHomo$PC1),100L),], aes(x=PC1, y=PC2, color='black'), size=0.1)
g <- g + geom_point(data=scores, aes(x=PC1, y=PC2), size=0.1)
g <- g + geom_point(data=projSplitIntoHomo, aes(x=PC1, y=PC2), size=1, shape=props$goodConvergence )
# http://www.perbang.dk/rgbgradient/
# g <- g + scale_colour_gradientn(colours = c('black', 'blue', 'red', '#D306DB', '#D721BC','#DC3C9E', '#E15880','#E57361', '#EA8E43' ,'#EFAA25','yellow' ))
g <- g + scale_colour_gradientn()

# g <- g + geom_point(data=projSplitIntoHomo, aes(x=PC1, y=PC2))
g <- g + theme_light() #+ facet_wrap('set')
# g <- g + scale_fill_manual(values=c("Cr", "Mn", "Fe", "Co"))
# g <- g + scale_size_manual(values = c("ss"=0.1, "homo"=3))
# g <- g + scale_shape_manual(values = c("ss"=1, "homo"=1))

# cairo_pdf(file="pca.pdf",width = 6, height = 5)
print(g)
# dev.off()

# g <- ggplot(data=totalScores, aes(x=PC1, y=PC2)) + stat_bin2d(binwidth = c(4,1)) + theme_light() +geom_point(data=scores,aes(color=metal))
# print(g)





