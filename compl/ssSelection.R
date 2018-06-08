# ss[ss$Eq %in% [result(!ISudo/spin)]$Eq,),]

loadSet <-function(name){
  nameStr = deparse(substitute(name))
  pipeStr <- paste("awk 'BEGIN{i=0}{i++;if (i%1==0) print $1}' < RAC_", nameStr, '.csv',sep='')
  name <- read.csv(pipe(pipeStr), header = TRUE)
  name$set <- nameStr
  return(name)
}

ssadc_racs <- loadSet(ssadc)
ssadc_b_racs <- loadSet(ssadc_b)
ssac_b_racs <- loadSet(ssac_b)
ss_racs <- rbind(rbind(ssadc_racs, ssadc_b_racs), ssac_b_racs)

ss_racs$eq <- gsub('_ax_.*', '', gsub('.._eq_', '', ss_racs$runs))
split_results$eq <- gsub('_ax1_.*', '', gsub('.._eq_smi', '', split_results$name))

validChildren <- ss_racs[ss_racs$eq %in% split_results$eq,]

ss_racs$eqn <- as.numeric(as.character(ss_racs$eq))
h <- ggplot(data=ss_racs,aes(x=eqn)) + geom_histogram(binwidth=0.5,center=0.5)
print(h)