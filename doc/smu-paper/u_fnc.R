require(ggplot2)

# octet
fun.1 <- function(x) 10+2*(8 - x)
fun.2 <- function(x) 10-1*(8 - x)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + stat_function(fun = fun.1, xlim = c(8,16))  + stat_function(fun = fun.2, xlim = c(0,8)) + xlim(0,16) +
  ylab("utility") +
  xlab("VE")

# charge
fun.1 <- function(x) 0
fun.2 <- function(x) 3
fun.3 <- function(x) 1
fun.4 <- function(x) 0

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + stat_function(fun = fun.1, xlim = c(0,2))  +
  stat_function(fun = fun.2, xlim = c(-2,0)) +
  stat_function(fun = fun.3, xlim = c(-3,-2)) +
  stat_function(fun = fun.4, xlim = c(-4,-3)) +
  xlim(-4,2) +
  ylab("utility") +
  xlab("charge")




