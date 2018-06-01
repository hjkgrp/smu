df <- as.data.frame(homolepticRF41KRRSpins)



g <- ggplot(data=df, aes(x=idu, y=e1, color='red'), size = 1) + geom_point()
g <- g + geom_point(data=df, aes(y=e2, color='blue'))
g <- g + geom_point(data=df, aes(y=e3, color='green'))
g <- g + geom_point(data=df, aes(y=e4, color='orange'))

print(g)
