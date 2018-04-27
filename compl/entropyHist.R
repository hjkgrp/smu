dent <- c('mono', 'mono', 'mono','mono','mono','mono','mono','mono','bi','bi','bi','bi','bi','bi','bi','bi')
set <- c('ho','fo','ssac', 'ssadc','ft','ws','eaac','eaadc')
ent <- c(19.7,13.7,1,12.7,12.7,8.1,1,1,15.63,1,9.47,8.53,9.47,1,10.04,1)
df <- data.frame(ent, dent, set)
# df$set <- factor(df$set, levels = df$set)
df$set <- factor(df$set, levels = (unique(df$set)), ordered=TRUE)

# melt the data frame for plotting
# df <- melt(df, id.vars='set')

# plot everything
gg<-ggplot(df, aes(x = set, y=ent)) +   theme_light() +
  geom_bar(aes(fill = dent), position = "dodge", stat="identity")+
  ylab('Entropy')+ xlab('Set')

cairo_pdf(file="file.pdf",width = 6, height = 5)
print(gg)
dev.off()

# set <- c('homoleptics','5+1','4+2', 'strongly sym.','eq. asym.','weakly sym')
# ent <- c(1,2,3,4,5,6)
# df <- data.frame(ent, set)
# df$set <- factor(df$set, levels = (unique(df$set)), ordered=TRUE)
# 
# # melt the data frame for plotting
# # df <- melt(df, id.vars='set')
# 
# # plot everything
# gg<-ggplot(df, aes(x = set, y=ent)) +   theme_light() +
#   geom_bar(position = "dodge", stat="identity")+
#   ylab('Entropy')+ xlab('Set')
# 
# # cairo_pdf(file="file.pdf",width = 6, height = 5)
# print(gg)
# # dev.off()
