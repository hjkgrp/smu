# read in example
# results_neut_m <- read.csv(file = "../compl/unified_results_neut_mults_post.csv")

# use results_neut_mults that come from process.R and contains only these: results_neut_mults <- results_neut_mults_neut_m (neut + mono) for iecr
racs_neut_m <- racs_neut_m
results_neut_m <- results_neut_m

racs_neut_m <- racs_neut_m[results_neut_m$alpha == '20', ]
results_neut_m <- results_neut_m[results_neut_m$alpha == '20', ]

high_score_lig <- c(9,14,21,136,232,237,293,150,153,155,161,163,166,238,245,307)

for (i in seq(1, dim(results_neut_m)[1])){
  if (as.character(results_neut_m$ox_2_HS_flag_oct[i]) == 'undef'){
    results_neut_m$ox_2_HS_flag_oct[i] = 0
  }
  if (as.character(results_neut_m$ox_2_LS_flag_oct[i]) == 'undef'){
    results_neut_m$ox_2_LS_flag_oct[i] = 0
  }
  if (as.character(results_neut_m$ox_3_LS_flag_oct[i]) == 'undef'){
    results_neut_m$ox_3_LS_flag_oct[i] = 0
  }
  if (as.character(results_neut_m$ox_3_HS_flag_oct[i]) == 'undef'){
    results_neut_m$ox_3_HS_flag_oct[i] = 0
  }
}

results_neut_m$goodConvergence <- (as.numeric(as.character(results_neut_m$ox_2_HS_flag_oct))) +
                       (as.numeric(as.character(results_neut_m$ox_2_LS_flag_oct))) +
                       (as.numeric(as.character(results_neut_m$ox_3_LS_flag_oct))) +
                       (as.numeric(as.character(results_neut_m$ox_3_HS_flag_oct)))

gg<-ggplot(data=results_neut_m, aes(x = goodConvergence)) + geom_histogram(binwidth = 1) + theme_light()
print(gg)
# plug in charge in netc

results_neut_m_20 <- results_neut_m[results_neut_m$eqn %in% high_score_lig & results_neut_m$goodConvergence == 4,]
racs_neut_m_20 <- racs_neut_m[results_neut_m$eqn %in% high_score_lig & results_neut_m$goodConvergence == 4,]

write.csv(results_neut_m_20,file="results_smu30.csv")
write.csv(racs_neut_m_20,file="racs_smu30.csv")
