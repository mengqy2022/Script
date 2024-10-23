lprob_poisson <- function(x) {
  n <- length(x)
  sum_x <- sum(x)
  c <- sum(lfactorial(x))
  
  function(lambda) {
    log(lambda) * sum_x - n * lambda - c
  }
}

# f1 <- lprob_poisson(x1)
# f1(10)
# 
# optimise(f1, c(0,100), maximum = T)
