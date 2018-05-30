df <- as.data.frame(homolepticRF41KRRspins$krrsplit[(0*405+1):(1*405)])
df$e1 <- as.data.frame(homolepticRF41KRRspins$krrsplit[(0*405+1):(1*405)])
df$e2 <- as.data.frame(homolepticRF41KRRspins$krrsplit[(1*405+1):(2*405)])
df$e3 <- as.data.frame(homolepticRF41KRRspins$krrsplit[(2*405+1):(3*405)])
df$e4 <- as.data.frame(homolepticRF41KRRspins$krrsplit[(3*405+1):(4*405)])
df$e5 <- as.data.frame(homolepticRF41KRRspins$krrsplit[(4*405+1):(5*405)])
df$e6 <- as.data.frame(homolepticRF41KRRspins$krrsplit[(5*405+1):(6*405)])
df$e7 <- as.data.frame(homolepticRF41KRRspins$krrsplit[(6*405+1):(7*405)])
df$e8 <- as.data.frame(homolepticRF41KRRspins$krrsplit[(7*405+1):(8*405)])

colnames(df)[1] <- 'e1_res'
colnames(df)[2] <- 'e1'
colnames(df)[3] <- 'e2'
colnames(df)[4] <- 'e3'
colnames(df)[5] <- 'e4'
colnames(df)[6] <- 'e5'
colnames(df)[7] <- 'e6'
colnames(df)[8] <- 'e7'
colnames(df)[9] <- 'e8'

df$idu <- as.numeric(row.names(df))

g <- ggplot(data=df, aes(x=idu, y=e1, color='red'), size = 1) + geom_point()
g <- g + geom_point(data=df, aes(y=e2, color='blue'))
g <- g + geom_point(data=df, aes(y=e3, color='green'))
g <- g + geom_point(data=df, aes(y=e4, color='orange'))

print(g)
