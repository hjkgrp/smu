# use results that come from process.R and contains only these: results <- results_m for evaluating the viability for children ss complexes
require(reshape2)

res <- results
racs_iecr <- racs

racs_iecr <- racs_iecr[res$alpha == 20, ]
res <- res[res$alpha == "20", ]

for (i in seq(1, dim(res)[1])){
  if (as.character(res$ox_2_HS_flag_oct[i]) == 'undef'){
    res$ox_2_HS_flag_oct[i] = 0
  }
  if (as.character(res$ox_2_LS_flag_oct[i]) == 'undef'){
    res$ox_2_LS_flag_oct[i] = 0
  }
  if (as.character(res$ox_3_LS_flag_oct[i]) == 'undef'){
    res$ox_3_LS_flag_oct[i] = 0
  }
  if (as.character(res$ox_3_HS_flag_oct[i]) == 'undef'){
    res$ox_3_HS_flag_oct[i] = 0
  }
}

res$goodConvergence <- (as.numeric(as.character(res$ox_2_HS_flag_oct))) +
  (as.numeric(as.character(res$ox_2_LS_flag_oct))) +
  (as.numeric(as.character(res$ox_3_LS_flag_oct))) +
  (as.numeric(as.character(res$ox_3_HS_flag_oct)))

lig_footprint <- read.csv('../enum/sp.out', header = FALSE)

for (i in seq(1, dim(res)[1])){
  for (j in seq(1, dim(lig_footprint)[1])){
    if (res$eqn[i] == j){
      if (j > 405){
        res$charge[i] <- lig_footprint[j,1] * 3
      } else {
        res$charge[i] <- lig_footprint[j,1] * 6
      }
    }
  }
}

ebi <- dcast(res, res$eqlig ~ res$goodConvergence)
ebi$gc_sum <- ebi$`0`*0 + ebi$`1`*1 + ebi$`2`*2 + ebi$`3`*3 + ebi$`4`*4
ebi$eqn <- as.numeric(gsub('smi','',ebi$`res$eqlig`))

res$charge <- as.factor(res$charge)
# , fill=charge
gg<-ggplot(data=ebi, aes(x = gc_sum)) + geom_histogram(binwidth = 1, position=position_dodge()) + theme_light()
# cairo_pdf(file="parentLigandsByCharge.pdf",width = 6, height = 5)
gg<-ggplot(data=ebi, aes(x = gc_sum)) + stat_bin(aes(y=cumsum(..count..)),geom="step") + theme_light()

print(gg)
# dev.off()


