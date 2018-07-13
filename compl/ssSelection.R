# ss[ss$Eq %in% [result(!ISudo/spin)]$Eq,),]

loadSet <-function(name){
  nameStr = deparse(substitute(name))
  pipeStr <- paste("awk 'BEGIN{i=0}{i++;if (i%1==0) print $1}' < RAC_", nameStr, '.csv',sep='')
  name <- read.csv(pipe(pipeStr), header = TRUE)
  name$set <- nameStr
  return(name)
}
# load footprint for ligands
lig_footprint <- read.csv('../enum/sp.out', header = FALSE)

# loading ss sets
ssadc_racs <- loadSet(ssadc)
ssadc_b_racs <- loadSet(ssadc_b)
ssac_b_racs <- loadSet(ssac_b)
ss_racs <- rbind(rbind(ssadc_racs, ssadc_b_racs), ssac_b_racs)

# replace the string format with actual numbers of the ligands and plot them
ss_racs$eq <- gsub('_ax_.*', '', gsub('.._eq_', '', ss_racs$runs))
ss_racs$ax <- gsub('.*_ax_', '', ss_racs$runs)
ss_racs$eqn <- as.numeric(as.character(ss_racs$eq))
ss_racs$axn <- as.numeric(as.character(ss_racs$ax))

# for (i in 1:dim(ss_racs)[1]){
#   if (ss_racs[i, 'ax'] > 404){
#     if (ss_racs[i, 'eqn'] > 404){
#       ss_racs[i, 'charge'] <- lig_footprint[ss_racs[i, 'eqn']+1, 1] * 2 + lig_footprint[ss_racs[i, 'axn']+1, 1] * 1
#     } else {
#       ss_racs[i, 'charge'] <- lig_footprint[ss_racs[i, 'eqn']+1, 1] * 4 + lig_footprint[ss_racs[i, 'axn']+1, 1] * 1
#     }
#   } else {
#     if (ss_racs[i, 'eqn'] > 404){
#       ss_racs[i, 'charge'] <- lig_footprint[ss_racs[i, 'eqn']+1, 1] * 2 + lig_footprint[ss_racs[i, 'axn']+1, 1] * 2
#     } else {
#       ss_racs[i, 'charge'] <- lig_footprint[ss_racs[i, 'eqn']+1, 1] * 4 + lig_footprint[ss_racs[i, 'axn']+1, 1] * 2
#     }
#   }
# }
# # filter out everything not in [-4,0]
# save(ss_racs,file="ss_racs_w_charge.Rdata")
# ss_racs_validCharges <- ss_racs[ss_racs$charge > -5 & ss_racs$charge < 1, ]
# 
# h1 <- ggplot(data=ss_racs,aes(x=eqn)) + geom_histogram(binwidth=0.5,center=0.5) + theme(axis.text.x=element_blank())
# # select only ss complexes that have working ligands in the DFT
# validChildren <- ss_racs_validCharges[ss_racs_validCharges$eq %in% split_results$eq,]
# h2 <- ggplot(data=validChildren, aes(x=eqn)) + geom_histogram(binwidth=0.5,center=0.5)
# 
# grid.arrange(h1, h2)

ebi$eqn0 < ebi$eqn - 1

validChildren <- ss_racs[ss_racs$eqn %in% (ebi[ebi$gc_sum > 10,]$eqn-1) & ss_racs$axn %in% (ebi[ebi$gc_sum > 5,]$eqn-1),]

dim(validChildren)

ggplot(data=validChildren, aes(x=eqn)) + geom_histogram(binwidth=0.5,center=0.5)

# now after clara is done and we have pa1000 with all the medoids:
validChildren_medoids <- validChildren[validChildren$runs %in% rownames(pa1000$medoids),]

run_metal_ax_eq <- cbind(validChildren_medoids$runs, sqrt(validChildren_medoids$mc.Z.0.all))
run_metal_ax_eq <- cbind(run_metal_ax_eq, validChildren_medoids$axn)
run_metal_ax_eq <- cbind(run_metal_ax_eq, validChildren_medoids$eqn)
colnames(run_metal_ax_eq) <- c('run', 'metal', 'ax', 'eq')
write.csv(run_metal_ax_eq,file="run_metal_ax_eq.csv")


