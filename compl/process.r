rm(list=ls())

source("fix_bidentate_names.R")

readResults <- function(path){
  results <- read.csv(file = path, stringsAsFactors = F, sep = ",")
}

# name rows and cols (ie but them into rownames) and rm the actual col
readRacs <- function(path){
  racs <- read.csv(file = path, header = TRUE, stringsAsFactors = F)
  rownames(racs)<-racs$runs
  racs$runs<-NULL
  return(racs)
}

# NEUTRAL MONODENTATES
results1 <- readResults('unified_results_post.csv')
racs1 <- readRacs("consistent_descriptor_file.csv")
r1 <- cbind(results1, racs1)

resultsNan1 <- readResults('unified_results_post_nans.csv') # they happen to be only monodentates
racsNan1 <- readRacs('consistent_descriptor_file_nans.csv')
rN1 <- cbind(resultsNan1, racsNan1)
rN1 <- rN1[,colnames(r1)] # just use the colnames that are present in r1, s.t. the rbind doesnt fail

for (i in seq(1,nrow(rN1))){
  nan_row <- rN1[i,]
  r1 <- r1[!(as.character(r1$name) == as.character(nan_row$name)), ]
}
r1 <- rbind(r1,rN1)

results1 <- r1[,colnames(results1)]
racs1 <- r1[,colnames(racs1)]

print(dim(results1[results1$alpha == 20,])[1]/4)

# CHARGED MONODENTATES (from JP)
results2 <- readResults('unified_results_post_charged.csv')
racs2 <- readRacs("consistent_descriptor_file_charged.csv") 
r2 <- cbind(results2, racs2)

resultsNan2 <- readResults('unified_results_post_charged_nans.csv') # they happen to be only monodentates
racsNan2 <- readRacs('consistent_descriptor_file_charged_nans.csv') 
rN2 <- cbind(resultsNan2, racsNan2)
rN2 <- rN2[,colnames(r2)] # just use the colnames that are present in r1, s.t. the rbind doesnt fail

for (i in seq(1,nrow(rN2))){
  nan_row = rN2[i,]
  r2 <- r2[!(as.character(r2$name) == as.character(nan_row$name)), ]
}
r2 <- rbind(r2, rN2)

results2 <- r2[,colnames(results2)]
racs2 <- r2[,colnames(racs2)]

print(dim(results2[results2$alpha == 20,])[1]/4)

# NEUTRAL BIDENTATES
results3 <- readResults('unified_results_post_bi.csv')
racs3 <- readRacs("consistent_descriptor_file_bi.csv")
r3 <- cbind(results3, racs3)

r3 <- fix_bidentate_names(r3)

rm <- NULL
for (i in seq(1, dim(r3)[1])){
  rm <- c(rm, !(r3[i,]$eqlig == 'zzz'))
}
r3 <- r3[rm,]

results3 <- r3[,colnames(results3)]
racs3 <- r3[,colnames(racs3)]

print(dim(results3[results3$alpha == 20,])[1]/4)

# CHARGED BIDENTATES (from JP)
results4 <- readResults('unified_results_post_charged_bi.csv')
racs4 <- readRacs('consistent_descriptor_file_charged_bi.csv') 
r4 <- cbind(results4, racs4)

r4 <- fix_bidentate_names(r4)

rm <- NULL
for (i in seq(1, dim(r4)[1])){
  rm <- c(rm, !(r4[i,]$eqlig == 'zzz'))
}
r4 <- r4[rm,]

results4 <- r4[,colnames(results4)]
racs4 <- r4[,colnames(racs4)]

print(dim(results4[results4$alpha == 20,])[1]/4)

# because some have  more cols because they were produced with a slightly earlier version of mAD
results1 <- results1[,colnames(results4)]
results2 <- results2[,colnames(results4)]
results3 <- results3[,colnames(results4)]
racs2 <- racs2[,colnames(racs1)]

# concat lists according to denticity
racs_m <- rbind(racs1, racs2)
racs_b <- rbind(racs3, racs4) 
results_m <- rbind(results1, results2)
results_b <- rbind(results3, results4)

#for iecr generate a list with only neutral monos
racs_neut_m <- racs1
results_neut_m <- results1
results_neut_m$eqn <- as.numeric(as.character(gsub('smi', '', results_neut_m$eqlig)))

# transform eqlig from format 'smiXXX' to a numeric (and add 405 to the bidentates. not anymore bc we apply the fix function upon import)
results_m$eqn <- as.numeric(as.character(gsub('smi', '', results_m$eqlig)))
results_b$eqn <- as.numeric(as.character(gsub('smi', '', results_b$eqlig))) + 0

# combine all results
racs <- rbind(racs_m, racs_b)
results <- rbind(results_m, results_b)


#CLARA only monodentates:
racs_clara <- racs_m
results_clara <- results_m

# generate all columns that are interesting to us
colstokeep <- c("name","gene","alpha","metal","ox2RN","ox3RN","ox_2_split","ox_3_split", 'eqn',"ox_2_HS_charge","ox_3_HS_charge")

# we iterate over ox_{2,3} and {LS,HS} because the string of the descriptor name is always constructed this way
for (props in c("energy","time","flag_oct","num_coord_metal", 'status')){
  for (oxs in c("ox_2","ox_3")){
    for (spins in c("LS","HS")){  
      colstokeep <- c(colstokeep,paste(oxs,spins,props,sep="_"))
      results[,paste(oxs,spins,props,sep="_")]<- as.numeric(as.character(results[,paste(oxs,spins,props,sep="_")]))
    }
  }
}


results[is.na(results$ox_2_HS_flag_oct), ]$ox_2_HS_flag_oct <- 0
results[is.na(results$ox_2_LS_flag_oct), ]$ox_2_LS_flag_oct <- 0
results[is.na(results$ox_3_HS_flag_oct), ]$ox_3_HS_flag_oct <- 0
results[is.na(results$ox_3_LS_flag_oct), ]$ox_3_LS_flag_oct <- 0 


# generate the splitting energies and a convergence check
results$betterConvergence <- results$ox_2_HS_flag_oct + results$ox_2_LS_flag_oct + results$ox_3_LS_flag_oct + results$ox_3_HS_flag_oct
results$ox_2_split <- results$ox_2_HS_energy - results$ox_2_LS_energy
results$ox_3_split <- results$ox_3_HS_energy - results$ox_3_LS_energy


# set up another df for ox_2 complexes that have all the needed energies non-nan 
ox_2_split_results <- results[!(is.na(results$ox_2_split)) & results$ox_2_HS_flag_oct == 1 &
                                results$ox_2_LS_flag_oct == 1 ,colstokeep]
ox_2_split_results$split  <- ox_2_split_results$ox_2_split
ox_2_split_results$ox_2_split <- NULL
ox_2_split_results$ox_3_split <- NULL
# ox_2_split_results$lfc<-ox_2_split_results$ox_2_HS_charge - 2 if you want to add charge to this spin split thing.

# select only racs for ox_2_split_result
matching_racs <- racs[as.character(ox_2_split_results$name),]
ox_2_split_results <- cbind(ox_2_split_results,matching_racs)
ox_2_split_results$ox <- 2
ox_2_split_results$rRN <- ox_2_split_results$ox2RN

# same for ox_3
ox_3_split_results <- results[!(is.na(results$ox_3_split)) & results$ox_3_HS_flag_oct == 1 & results$ox_3_LS_flag_oct == 1 ,colstokeep]
ox_3_split_results$split  <- ox_3_split_results$ox_3_split
ox_3_split_results$ox_3_split<-NULL
ox_3_split_results$ox_2_split<-NULL
matching_racs<-racs[as.character(ox_3_split_results$name),]
ox_3_split_results<-cbind(ox_3_split_results,matching_racs)
ox_3_split_results$ox <-3
ox_3_split_results$rRN <-ox_3_split_results$ox3RN

# put together all valid splitting rows
split_results <- rbind(ox_2_split_results,ox_3_split_results)
require(plyr)
split_results$metal <- revalue(as.factor(split_results$metal),c("0"="Cr","1"="Mn","2"="Fe","3"="Co"))
split_results$ox <- as.factor(split_results$ox)


# require(ggplot2)
# gg<-ggplot(data=split_results,aes(x=V104,y=split,color=ox,group=ox)) +geom_point() + facet_grid("metal ~.") + theme_light()
# print(gg)
# 
# require(reshape2) # enables melt (wide to long format) and cast (vv). plyr and ggplot2 both use long format.
# results.m <- melt(results, measure.vars = c("ox_2_LS_time","ox_3_LS_time","ox_2_HS_time","ox_3_HS_time"), na.rm = TRUE)
# results.m$betterConvergence <- as.factor(results.m$betterConvergence)
# gg<-ggplot(data=results.m, aes(x = value / 3600, group = betterConvergence, fill = betterConvergence)) + 
#     geom_histogram(position = position_dodge(), bins=20) + theme_light()
# print(gg)
# 
# 
# results$metal<-as.factor(results$metal)
# gg<-ggplot(data=results) + geom_point(size=4,aes(y=ox_2_LS_flag_oct,x=ox_2_LS_time/3600,color=metal),shape=1) +
#                            geom_point(size=4,aes(y=ox_2_HS_flag_oct,x=ox_2_HS_time/3600,color=metal),shape=2) +
#                            geom_point(size=4,aes(y=ox_3_LS_flag_oct,x=ox_3_LS_time/3600,color=metal),shape=3) +
#                            geom_point(size=4,aes(y=ox_3_HS_flag_oct,x=ox_3_HS_time/3600,color=metal),shape=4) +
#                            geom_text(size=4,aes(y=ox_2_LS_flag_oct,x=ox_2_LS_time/3600,color=metal,label=ox2RN),
#                                       position = position_jitter(),data=results[results$ox_2_LS_time>36000, ]) +
#                            geom_text(size=4,aes(y=ox_3_LS_flag_oct,x=ox_3_LS_time/3600,color=metal,label=ox3RN),
#                                       position = position_jitter(),data=results[results$ox_3_LS_time>36000, ]) +
#                            theme_light()
# print(gg)
# #df <- results[,]
# # df$ss3 <- results$ox_3_HS_energy - results$ox_3_LS_energy
# # df$ipls <- results$ox_3_LS_energy - results$ox_2_LS_energy
# # df$iphs <- results$ox_3_HS_energy - results$xo_2_HS_energy
