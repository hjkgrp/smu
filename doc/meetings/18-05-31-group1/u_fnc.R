require(ggplot2)

fun.1 <- function(x) 10+2*(8 - x)
fun.2 <- function(x) 10-1*(8 - x)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

# octet
p + stat_function(fun = fun.1, xlim = c(8,16))  + stat_function(fun = fun.2, xlim = c(0,8)) + xlim(0,16) +
  ylab("utility") +
  xlab("VE")

# charge




