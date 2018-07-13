fix_bidentate_names<-function(df){
  ## build map
  seq1 <- paste('smi',as.character(seq(0,148)),sep="") # current labs
  seq2 <- paste('smi',as.character(seq(405+1,405+128)),sep="") # target labs
  bad_list <- paste('smi',as.character(c(21-1, 22-1,  52-1,  53-1,  54-1,  
                                         55-1,  73-1,  74-1,  77-1,  81-1,
                                         85-1,  86-1, 116-1, 117-1, 118-1,
                                         119-1, 137-1, 138-1, 141-1, 145-1)),sep="")
  
  count <-0
  final_mapping <- list()
  for (i in seq1){
    if (i %in% bad_list)
    {
      final_mapping[i]<-'zzz'
    }
    else{
      count <- count + 1 
      final_mapping[i]<-seq2[count]
    }
  }
  
  
  ## make new df 
  newdf <- df
  ## NOTE we need to convert the type to char from factor, or 
  ## else add the additional levels to the factor variable
  newdf$eqlig <- factor(newdf$eqlig,levels = c(seq1,seq2,'zzz'))
  newdf$axlig1 <- factor(newdf$axlig1,levels = c(seq1,seq2,'zzz'))
  newdf$axlig2 <- factor(newdf$axlig2,levels = c(seq1,seq2,'zzz'))
  
  ## now fix ligands
  for (i in seq(1,nrow(newdf))){
    if (as.character(newdf[i,]$eqlig) %in% names(final_mapping) ){
      newdf[i,]$eqlig<-final_mapping[newdf[i,]$eqlig]
      newdf[i,]$axlig1<-final_mapping[newdf[i,]$axlig1]
      newdf[i,]$axlig2<-final_mapping[newdf[i,]$axlig2]}
    else{
      print('Error! key not in final mapping!')
      print(newdf[i,]$eqlig)
    }
  }
  
  ## now fix genes
  require(plyr)
  newdf$temp_metal <- factor(newdf$metal)
  newdf$temp_metal <- revalue(newdf$temp_metal,
                              c("0"="cr","1"="mn","2"="fe","3"="co"))
  
  newdf$newgene =  paste(newdf$temp_metal,
                         newdf$eqlig,newdf$axlig1,newdf$axlig2,sep="_")
  
  newdf$gene <- newdf$newgene
  
  ## clean up
  newdf$temp_metal<-NULL
  newdf$newgene<-NULL
  return(newdf)
}


# ## load df
# df<-read.csv(file = "unified_results_post_charged_bi.csv")
# 
# ## test function
# newdf<-fix_bidentate_names(df)
# 
# ## check results
# require(ggplot2)
# require(gtools)
# plotdf <- newdf
# plotdf$oldlig <- df$eqlig
# plotdf$eqlig<-factor(plotdf$eqlig,levels=mixedsort(levels(plotdf$eqlig)))
# plotdf$oldlig<-factor(plotdf$oldlig,levels=mixedsort(levels(plotdf$oldlig)))
# 
# g<- ggplot(data=plotdf,aes(x=eqlig,y=oldlig)) + geom_point() +
#   xlab("new") +ylab('old') +theme(axis.text.x = element_text(angle=90))
# print(g)


