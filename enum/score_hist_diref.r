require(ggplot2)

dat_di <- read.csv('DHAL_data.csv', header = FALSE)
dat <- dat_di[, c(1, 28)]
dat$set <- 'di' 
colnames(dat) <- c('smiles','score','set')

dat_mo <- read.csv('MHAL_data.csv', header = FALSE)
dat_mo$V16 <- dat_mo$V16 / 10 * 17
dat_mo <- dat_mo[, c(1, 16)]
dat_mo$set <- 'mo' 
colnames(dat_mo) <- c('smiles','score', 'set')

dat_bi <- read.csv('THAL_data.csv', header = FALSE)
dat_bi$V28 <- dat_bi$V28 / 16 * 17
dat_bi <- dat_bi[, c(1, 28)]
dat_bi$set <- 'bi'
colnames(dat_bi) <- c('smiles','score','set')

dat <- rbind(dat, dat_mo)
dat <- rbind(dat, dat_bi)

diref <- read.csv('diref_occ', header = FALSE)
gdb4 <- read.csv('gdb-4', header = FALSE)
chembl <- read.csv('chembl_smiles', header = FALSE)

dat$diatomic <- as.factor(as.integer(dat[,1] %in% diref[,1]))
dat$gdb4 <- as.factor(as.integer(dat[,1] %in% gdb4[,1]))
dat$chembl <- as.factor(as.integer(dat[,1] %in% chembl[,1]))

j = 0
dat$diref <- NULL
dat$direfNumber <- NULL
for (i in seq(1, dim(dat)[1])){
  if (dat$diatomic[i] == 1){
    j = j + 1
    dat$diref[i] <- as.integer(as.logical(diref[j,2]))
    dat$direfNumber[i] <- diref[j,2]
  }
  else{
    dat$diref[i] <- 0 # could be 2 (or anoither number except 0/1 to denote non-diatoms)
    dat$direfNumber[i] <- 0
  }
}

# is it known to either 
for (i in seq(1, dim(dat)[1])){
  if (dat$diref[i] == 1 || dat$gdb4[i] == 1 || dat$chembl[i] == 1){  
    dat$known[i] <- 1
  }
  else{
    dat$known[i] <- 0
  }
}


dat$known <- factor(dat$known, levels=c('1','0'))
g <- ggplot(data=dat, aes(x=score, fill=known)) + geom_histogram(position=position_stack(), bins = 10, binwidth = 1) +
  xlab('score') + scale_fill_discrete(name = "known") + xlim(0,17)

h <- ggplot(data=dat[dat$direfNumber != 0,], aes(x=score, y=direfNumber)) + geom_point() +   geom_smooth(method=lm) +
  xlab('score') + xlim(0,17)

print(h)

