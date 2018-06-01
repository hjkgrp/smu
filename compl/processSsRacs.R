# rm(list=ls())
require(plyr)
require(dplyr)
require(ggplot2)
HF_to_kcalmol = 627.5095
loadSet <-function(name){
  nameStr = deparse(substitute(name))
  if (nameStr != 'ho'){
    pipeStr <- paste("awk 'BEGIN{i=0}{i++;if (i%100==1) print $1}' < RAC_", nameStr, '.csv',sep='')
    name <- read.csv(pipe(pipeStr), header = TRUE)
  } else {
    name <- read.csv(file = paste("RAC_",nameStr,".csv", sep=''), header = TRUE)
  }
  name$set <- nameStr
  return(name)
}

standardizeSet <-function(rf133racs, center, scale){
  # standardize data with homo moments
  rf133racs_data <- scale(rf133racs, center=center, scale=scale)
  rf133racs_data <- as.data.frame(rf133racs_data)
  return(rf133racs_data)
}

homo_racs <- loadSet(ho)
homo_racs$set[1621:2212] <- 'ho_b'

ssadc_racs <- loadSet(ssadc)
ssadc_b_racs <- loadSet(ssadc_b)
ssac_b_racs <- loadSet(ssac_b)
ssracs <- rbind(ssadc_racs, ssadc_b_racs)
ss_racs <- rbind(rbind(ssadc_racs, ssadc_b_racs), ssac_b_racs)

fo_racs <- loadSet(fo)

ft_m_racs <-loadSet(ft) # contains only mono
ft_b_racs <- loadSet(ft_b)
ft_racs <- rbind(ft_m_racs, ft_b_racs)

ws_m_racs <- loadSet(ws) # contains random sample of ws monodentates (150k)
ws_b_racs <- loadSet(ws_b) # contains random sample of ws bidentates
ws_racs <- rbind(ws_m_racs, ws_b_racs)


superset <- do.call("rbind", list(homo_racs, ss_racs, fo_racs, ft_racs, ws_racs))

# read in descriptors
load(file="selected_all.R")
load(file="selected_rf_41.Rda")

# remove alpha and ox from "rf41" and "all" because we only want topology related racs for the fingerprint.
rf39vars <- selected_rf_41[!(selected_rf_41 %in% c('alpha','ox'))]
rf153vars <- selected_all[!(selected_all %in% c('alpha','ox'))]

homo_rf39racs <- homo_racs[,rf39vars]
homo_rf153racs <- homo_racs[,rf153vars]

ssadc_rf39racs <- ssadc_racs[,rf39vars]
ssadc_rf153racs <- ssadc_racs[,rf153vars]
ssadc_b_rf39racs <- ssadc_b_racs[,rf39vars]
ssadc_b_rf153racs <- ssadc_b_racs[,rf153vars]
ssac_b_rf39racs <- ssac_b_racs[,rf39vars]
ssac_b_rf153racs <- ssac_b_racs[,rf153vars]
ss_rf39racs <- ss_racs[,rf39vars]
ss_rf153racs <- ss_racs[,rf153vars]

fo_rf39racs <- fo_racs[,rf39vars]
fo_rf153racs <- fo_racs[,rf153vars]

ft_m_rf39racs <- ft_m_racs[,rf39vars]
ft_m_rf153racs <- ft_m_racs[,rf153vars]
ft_b_rf39racs <- ft_b_racs[,rf39vars]
ft_b_rf153racs <- ft_b_racs[,rf153vars]
ft_rf39racs <- ft_racs[,rf39vars]
ft_rf153racs <- ft_racs[,rf153vars]

ws_m_rf39racs <- ws_m_racs[,rf39vars]
ws_b_rf39racs <- ws_b_racs[,rf39vars]
ws_m_rf153racs <- ws_m_racs[,rf153vars]
ws_b_rf153racs <- ws_b_racs[,rf153vars]
ws_rf39racs <- ws_racs[,rf39vars]
ws_rf153racs <- ws_racs[,rf153vars]

superset_rf39racs <- superset[,rf39vars]
superset_rf153racs <- superset[,rf153vars]

split_results_39racs <- split_results[,rf39vars]
split_results_153racs <- split_results[,rf153vars]
props_conv <- results[,!(colnames(results) %in% (rf153vars))]
props_conv$goodConvergence <- results$betterConvergence
props <- (split_results[,!(colnames(split_results) %in% (rf153vars))])
props$goodConvergence <- as.factor(props$ox_2_HS_flag_oct + props$ox_2_LS_flag_oct + props$ox_3_LS_flag_oct + props$ox_3_HS_flag_oct)


# Apply the variance function 'var' on all cols ('2' means cols) of rf39racs. Throw out the one's with zero variance for PCA later.
# which(apply(rf153racs, 2, var)==0)
# apply the same deletion in the ss and fo even though they have actually more, to get same dimensionality
# 153
superset_rf133racs <- superset_rf153racs[ , apply(superset_rf153racs, 2, var) != 0]
superset_setnames <-  superset_rf153racs[ , apply(superset_rf153racs, 2, var) != 0]
homo_rf133racs <- homo_rf153racs[ , apply(superset_rf153racs, 2, var) != 0]
ss_rf133racs   <- ss_rf153racs[ , apply(superset_rf153racs, 2, var) != 0]
ft_rf133racs <- ft_rf153racs[ , apply(superset_rf153racs, 2, var) != 0]
fo_rf133racs   <- fo_rf153racs[ , apply(superset_rf153racs, 2, var) != 0]
ws_rf133racs   <- ws_rf153racs[ , apply(superset_rf153racs, 2, var) != 0]
split_results_rf133racs <- split_results_153racs[ , apply(superset_rf153racs, 2, var) != 0]

# 39
superset_rf39racs <- superset_rf39racs[ , apply(superset_rf39racs, 2, var) != 0]
homo_rf39racs <- homo_rf39racs[ , apply(superset_rf39racs, 2, var) != 0]
ss_rf39racs   <- ss_rf39racs[ , apply(superset_rf39racs, 2, var) != 0]
ft_rf39racs <- ft_rf39racs[ , apply(superset_rf39racs, 2, var) != 0]
fo_rf39racs   <- fo_rf39racs[ , apply(homo_rf39racs, 2, var) != 0]
ws_rf39racs   <- ws_rf39racs[ , apply(homo_rf39racs, 2, var) != 0]
split_results_rf39racs <- split_results_39racs[ , apply(superset_rf39racs, 2, var) != 0]


###############
# Standardize #
###############

# find moments of standardization for the superset
# 155
superset_rf133racs_data <- scale(superset_rf133racs, center=TRUE, scale=TRUE)
superset_rf133racs_center <-attr(superset_rf133racs_data, 'scaled:center')
superset_rf133racs_scale <- attr(superset_rf133racs_data, 'scaled:scale')
superset_rf133racs_data <- as.data.frame(superset_rf133racs_data)
# 39
superset_rf39racs_data <- scale(superset_rf39racs, center=TRUE, scale=TRUE)
superset_rf39racs_center <-attr(superset_rf39racs_data, 'scaled:center')
superset_rf39racs_scale <- attr(superset_rf39racs_data, 'scaled:scale')
superset_rf39racs_data <- as.data.frame(superset_rf39racs_data)

# stdize sets
# 155
homo_rf133racs_data <- standardizeSet(homo_rf133racs, superset_rf133racs_center, superset_rf133racs_scale)
ss_rf133racs_data <- standardizeSet(ss_rf133racs, superset_rf133racs_center, superset_rf133racs_scale)
fo_rf133racs_data <- standardizeSet(fo_rf133racs, superset_rf133racs_center, superset_rf133racs_scale)
ft_rf133racs_data <- standardizeSet(ft_rf133racs, superset_rf133racs_center, superset_rf133racs_scale)
ws_rf133racs_data <- standardizeSet(ws_rf133racs, superset_rf133racs_center, superset_rf133racs_scale)
split_results_rf133racs_data <- standardizeSet(split_results_rf133racs, superset_rf133racs_center, superset_rf133racs_scale)
# 39
homo_rf39racs_data <- standardizeSet(homo_rf39racs, superset_rf39racs_center, superset_rf39racs_scale)
ss_rf39racs_data <- standardizeSet(ss_rf39racs, superset_rf39racs_center, superset_rf39racs_scale)
fo_rf39racs_data <- standardizeSet(fo_rf39racs, superset_rf39racs_center, superset_rf39racs_scale)
ft_rf39racs_data <- standardizeSet(ft_rf39racs, superset_rf39racs_center, superset_rf39racs_scale)
ws_rf39racs_data <- standardizeSet(ws_rf39racs, superset_rf39racs_center, superset_rf39racs_scale)
split_results_rf39racs_data <- standardizeSet(split_results_39racs, superset_rf39racs_center, superset_rf39racs_scale)

############################
# PCA with function prcomp #
############################
pca1 = prcomp(superset_rf133racs_data,scale. = FALSE) #scale is FALSE because we standardize above
pca2 = prcomp(superset_rf39racs_data,scale. = FALSE) #scale is FALSE because we standardize above

# use only descriptors with < 3 length:
pca3 = prcomp(superset_rf39racs_data[,],scale. = FALSE)

# sqrt of eigenvalues
eigendecay1 = pca1$sdev
eigendecay2 = pca2$sdev
plot(eigendecay2)

# PCs (aka scores) [PC1, PC2, ... data points]
scores = as.data.frame(pca1$x)
scores2 = as.data.frame(pca2$x)

# >>>>>>>>>>155
# mc.Z.0.all is the metal charge and since autocorrelation goes Z1*Z2, we sqrt it
scores$metal <- revalue(as.factor(sqrt(superset_rf133racs$mc.Z.0.all)),c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
# the sqrt of Z1Z2 of a ligand centric first atom is the CA. below we have a cont spectr of atoms bc avg in ss is an atom in R
#scores$cai <- revalue(as.factor(sqrt(homo_rf133racs$lc.Z.0.eq)),c('6'='C','7'='N','8'='O','15'='P','16'='S'))
scores$cai <- (4*sqrt(superset_rf133racs$lc.Z.0.eq)+2*sqrt(superset_rf133racs$lc.Z.0.ax))/6

# projection of sets into SU(perset) coordinates and reassignments
projHomoIntoSu <- as.data.frame(predict(pca1, homo_rf133racs_data))
projHomoIntoSu$set <- 'homo'
projHomoIntoSu$metal <- revalue(as.factor(sqrt(homo_rf133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projHomoIntoSu$cai <- (4*sqrt(homo_rf133racs$lc.Z.0.eq) + 2*sqrt(homo_rf133racs$lc.Z.0.ax))/6
# is it right that we dont use the normalized data for the last 2 lines???

projSsIntoSu <- as.data.frame(predict(pca1, ss_rf133racs_data))
projSsIntoSu$set <- 'ss'
projSsIntoSu$metal <- revalue(as.factor(sqrt(ss_rf133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projSsIntoSu$cai <- (4*sqrt(ss_rf133racs$lc.Z.0.eq) + 2*sqrt(ss_rf133racs$lc.Z.0.ax))/6
# is it right that we dont use the normalized data for the last 2 lines???

projFoIntoSu <- as.data.frame(predict(pca1, fo_rf133racs_data))
projFoIntoSu$set <- 'fo'
projFoIntoSu$metal <- revalue(as.factor(sqrt(fo_rf133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projFoIntoSu$cai <- (4*sqrt(fo_rf133racs$lc.Z.0.eq) + 2*sqrt(fo_rf133racs$lc.Z.0.ax))/6

projFtIntoSu <- as.data.frame(predict(pca1, ft_rf133racs_data))
projFtIntoSu$set <- 'ft'
projFtIntoSu$metal <- revalue(as.factor(sqrt(ft_rf133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projFtIntoSu$cai <- (4*sqrt(ft_rf133racs$lc.Z.0.eq) + 2*sqrt(ft_rf133racs$lc.Z.0.ax))/6

projWsIntoSu <- as.data.frame(predict(pca1, ws_rf133racs_data))
projWsIntoSu$set <- 'ws'
projWsIntoSu$metal <- revalue(as.factor(sqrt(ws_rf133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projWsIntoSu$cai <- (4*sqrt(ws_rf133racs$lc.Z.0.eq) + 2*sqrt(ws_rf133racs$lc.Z.0.ax))/6

projSplitIntoSu <- as.data.frame(predict(pca1, split_results_rf133racs_data))
projSplitIntoSu$set <- 'calcHomo'
projSplitIntoSu$metal <- revalue(as.factor(sqrt(split_results_rf133racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projSplitIntoSu$cai <- (4*sqrt(split_results_rf133racs$lc.Z.0.eq) + 2*sqrt(split_results_rf133racs$lc.Z.0.ax))/6
# <<<<<<<<<<<<<<<<155

#>>>>>>>>>>>>>>>>>39
# mc.Z.0.all is the metal charge and since autocorrelation goes Z1*Z2, we sqrt it
scores2$metal <- revalue(as.factor(sqrt(superset_rf39racs$mc.Z.0.all)),c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
# the sqrt of Z1Z2 of a ligand centric first atom is the CA. below we have a cont spectr of atoms bc avg in ss is an atom in R
#scores$cai <- revalue(as.factor(sqrt(homo_rf39racs$lc.Z.0.eq)),c('6'='C','7'='N','8'='O','15'='P','16'='S'))
scores2$cai <- (4*sqrt(superset_rf39racs$lc.Z.0.eq)+2*sqrt(superset_rf39racs$lc.Z.0.ax))/6
#scores2$set = 'SU'
scores2$set = superset_setnames

# projection of sets into SU(perset) coordinates and reassignments
projHomoIntoSu2 <- as.data.frame(predict(pca2, homo_rf39racs_data))
projHomoIntoSu2$set <- 'homo'
projHomoIntoSu2$metal <- revalue(as.factor(sqrt(homo_rf39racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projHomoIntoSu2$cai <- (4*sqrt(homo_rf39racs$lc.Z.0.eq) + 2*sqrt(homo_rf39racs$lc.Z.0.ax))/6
projHomoIntoSu2$set <- homo_racs$set
# is it right that we dont use the normalized data for the last 2 lines???

projSsIntoSu2 <- as.data.frame(predict(pca2, ss_rf39racs_data))
projSsIntoSu2$set <- 'ss'
projSsIntoSu2$metal <- revalue(as.factor(sqrt(ss_rf39racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projSsIntoSu2$cai <- (4*sqrt(ss_rf39racs$lc.Z.0.eq) + 2*sqrt(ss_rf39racs$lc.Z.0.ax))/6
# is it right that we dont use the normalized data for the last 2 lines???
projSsIntoSu2$set <- ss_racs$set

projFoIntoSu2 <- as.data.frame(predict(pca2, fo_rf39racs_data))
projFoIntoSu2$set <- 'fo'
projFoIntoSu2$metal <- revalue(as.factor(sqrt(fo_rf39racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projFoIntoSu2$cai <- (4*sqrt(fo_rf39racs$lc.Z.0.eq) + 2*sqrt(fo_rf39racs$lc.Z.0.ax))/6
projFoIntoSu2$set <- fo_racs$set

projFtIntoSu2 <- as.data.frame(predict(pca2, ft_rf39racs_data))
projFtIntoSu2$set <- 'ft'
projFtIntoSu2$metal <- revalue(as.factor(sqrt(ft_rf39racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projFtIntoSu2$cai <- (4*sqrt(ft_rf39racs$lc.Z.0.eq) + 2*sqrt(ft_rf39racs$lc.Z.0.ax))/6
projFtIntoSu2$set <- ft_racs$set

projWsIntoSu2 <- as.data.frame(predict(pca2, ws_rf39racs_data))
projWsIntoSu2$set <- 'ws'
projWsIntoSu2$metal <- revalue(as.factor(sqrt(ws_rf39racs$mc.Z.0.all)), c('24'='Cr','25'='Mn','26'='Fe','27'='Co'))
projWsIntoSu2$cai <- (4*sqrt(ws_rf39racs$lc.Z.0.eq) + 2*sqrt(ws_rf39racs$lc.Z.0.ax))/6
projWsIntoSu2$set <- ws_racs$set

projSplitIntoSu2 <- as.data.frame(predict(pca2, split_results_rf39racs_data))
projSplitIntoSu2$set <- 'splitHomo'              
projSplitIntoSu2$metal <- split_results$metal
projSplitIntoSu2$cai <- (4*sqrt(split_results$lc.Z.0.eq) + 2*sqrt(split_results$lc.Z.0.ax))/6
#<<<<<<<<<<<<<< 39

# totalScores2 <- rbind(projSsIntoSu2, scores2)

# # electronegativity 
# # proxy for number of atoms
# totalScores$nats <- homo_rf133racs$f.chi.0.all
cai_colors <- projSsIntoHomo$cai[seq.int(1L,length(projSsIntoHomo$cai),100L)]
g <- ggplot(data=projSsIntoHomo[seq.int(1L,length(projSsIntoHomo$PC1),100L),],
            aes(x=PC1, y=PC2, color=cai_colors)) +
     geom_point(size = 0.1)
# g <- g + geom_point(data=scores, aes(x=PC1, y=PC2, color='red'), size=0.1)
g <- g + geom_point(data=projSplitIntoHomo, aes(x=PC1, y=PC2), size=1, color=props$goodConvergence )
# http://www.perbang.dk/rgbgradient/
# g <- g + scale_colour_gradientn(colours = c('black', 'blue', 'red', '#D306DB', '#D721BC','#DC3C9E', '#E15880','#E57361', '#EA8E43' ,'#EFAA25','yellow' ))
g <- g + scale_colour_gradientn()
g <- g + theme_light() #+ facet_wrap('set')
# g <- g + scale_fill_manual(values=c("Cr", "Mn", "Fe", "Co"))
# g <- g + scale_size_manual(values = c("ss"=0.1, "homo"=3))
# g <- g + scale_shape_manual(values = c("ss"=1, "homo"=1))
# cairo_pdf(file="pca_monoANDbi_SsIntoHomo.pdf",width = 6, height = 5)
print(g)
# dev.off()
g <- ggplot(data=totalScores, aes(x=PC1, y=PC2)) + stat_bin2d(binwidth = c(4,1)) + theme_light() + geom_point(data=scores,aes(color=metal))
print(g)

# PCA of calc'd splitting energies into SU colored by cai/metal/ene
plot_df <- projSplitIntoSu2
plot_df$split<-split_results$split
plot_df$split<-plot_df$split*HF_to_kcalmol 
plot_df$ox<-split_results$ox
#plot_df_2 <- scores2[seq.int(1L,length(scores2$PC1),10L),]
#plot_df_2$split <- df$krrsplit
g <- ggplot(data=plot_df[plot_df$ox==2 & plot_df$split < 150, ],
            aes(x=PC1, y=PC2, fill = split,shape=metal)) +scale_fill_gradient2(high='firebrick',low='dodgerblue',mid='white',midpoint =0)+
        geom_point(data=scores2[seq.int(1L,length(scores2$PC1),10L),], shape=15,size = 1.15,alpha=0.5,fill='purple',aes(color=cai)) +
        geom_point(size = 5,color='black') + theme_light()
g <- g + ggtitle("PCA of all sets with homoleptic splitting energy results superimposed") + 
  theme(title = element_text(hjust = 0.5)) + scale_shape_manual(values= c("Cr"=21,"Co"=22,"Fe"=23,"Mn"=24))+
  labs(fill='spliting energy')  + labs(color='connecting atom') + theme(panel.background = element_rect(fill='white',color='black'))
g <- g + scale_colour_gradientn(colours = c('black', 'blue', 'red', '#D306DB', 
                                            '#D721BC','#DC3C9E', '#E15880','#E57361', '#EA8E43' ,'#EFAA25','yellow' ))
  
# cairo_pdf(file="pca_rf39_SplitIntoSU_convergence.pdf",width = 6, height = 5)
print(g)
# dev.off()

# PCA of SS into SU colored by metal/cai
g <- ggplot(data=projSsIntoSu2[seq.int(1L,length(projSsIntoSu2$PC1),10L),],
            aes(x=PC1, y=PC2, color = cai)) +
  geom_point(size = 2)
g <- g + geom_point(data=scores2[seq.int(1L,length(scores2$PC1),10L),], size = 0.2)
g <- g + guides(color = guide_legend(override.aes = list(size = 5)))
g <- g + ggtitle("PCA of all sets with SS projected into it") + 
  labs(color='Convergence')  
# cairo_pdf(file="pca_rf39_SsIntoSU_cai.pdf",width = 6, height = 5)
print(g)
# dev.off()

# PCA of different sets into SU with facet = set
g <- ggplot(data=projHomoIntoSu2[seq.int(1L,length(projHomoIntoSu2$PC1),1L),],
            aes(x=PC1, y=PC2, color=cai)) +
  geom_point(size = 2, color='red')
g <- g + geom_point(data=scores2[seq.int(1L,length(scores2$PC1),10L),], size = 0.2)
g <- g + guides(color = guide_legend(override.aes = list(size = 5)))
# g <- g + ggtitle("PCA of all sets with SS projected into it") + 
  labs(color='Convergence')  #+ facet_wrap('set') 
cairo_pdf(file="pca_rf39_HomoIntoSU_cai_span.pdf",width = 6, height = 5)
print(g)
dev.off()

# -- t-SNE

require(tsne)
colors = rainbow(length(unique(split_[,40])))
names(colors) = unique(split_[,40])
ecb = function(x,y){ plot(x,t='n'); text(x,labels=split_[,40], col=colors[split_[,40]]) }

cairo_pdf(file="tsne_50_rf39_split",width = 6, height = 5)
conv_iris = tsne(split_[,1:39], epoch_callback = ecb, perplexity=50)
dev.off()

# --
g <- ggplot(data=totalScores2, aes(x=PC1, y=PC2)) + stat_bin2d(binwidth = c(4,1)) + theme_light() 
#g <- g + geom_point(data=scores2[seq.int(1L,length(scores2$PC1),100L),],aes(color=metal))
print(g)

