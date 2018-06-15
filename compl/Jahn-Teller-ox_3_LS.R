results$ox_3_LS_ax1_MLB<-as.numeric(as.character(results$ox_3_LS_ax1_MLB))
results$ox_3_LS_eq_MLB<-as.numeric(as.character(results$ox_3_LS_eq_MLB))
results$metal<-factor(results$metal)
p<-ggplot(data=results[!(is.na(results$ox_3_LS_eq_MLB)) &!(is.na(results$ox_3_LS_ax1_MLB))& results$alpha==20 &
                         results$ox_3_LS_flag_oct ==1,],
          aes(x=ox_3_LS_eq_MLB,y=ox_3_LS_ax1_MLB,color=metal)) +geom_abline(color='gray')+ 
  geom_label_repel(aes(label = eqlig),
                   data =results[!(is.na(results$ox_3_LS_eq_MLB)) &!(is.na(results$ox_3_LS_ax1_MLB))& results$alpha==20 &
                                   results$ox_3_LS_flag_oct ==1 & abs(results$ox_3_LS_eq_MLB - results$ox_3_LS_ax1_MLB)> 0.1,] )+
  geom_point(size=4)
print(p)