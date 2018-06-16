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


# entropy demo
p9 <- ggplot(data.frame(x = c(0, 4)), aes(x = x)) +
  stat_function(fun = dnorm, args = list(0.5, 0.01),
                aes(colour = "Low Entropy"), size = 1.5) +
  
  stat_function(fun = dnorm, args = list(1.5, 0.1),
                aes(colour = "Meh Entropy"), size = 1.5) +
  
  stat_function(fun = dnorm, args = list(2.5, 0.2),
                aes(colour = "Nice Entropy"), size = 1.5) +
  
  stat_function(fun = dnorm, args = list(3.5, 0.7),
                aes(colour = "High Entropy"), size = 1.5) +
  
  scale_x_continuous(name = "Entropy",
                     breaks = seq(0, 4, 0.2),
                     limits=c(0, 4)) +
  scale_y_continuous(name = "A.U.") +
  ggtitle("Entropy Comparison") +
  scale_colour_brewer(palette="Accent") +
  theme_bw() + theme(legend.position="none") +
  annotate("text", x = 0.9, y = 5, label = "Low Entropy") +
  annotate("text", x = 3.6, y = 1, label = "High Entropy")

p9
