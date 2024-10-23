# geom_histogram ----------------------------------------------------------
#  20240727

histogram_bins <- function(type) {

  fun <- switch (type,
                 Sturges = nclass.Sturges,
                 scott = nclass.scott,
                 FD = nclass.FD,
                 stop("Unknown type", call. = F)
  )
  
  function(x) {
    (max(x) - min(x)) / fun(x)
  }
}

#  直方图的bin
# sd <- c(1, 5, 15)
# n <- 100
# df <- data.frame(x = rnorm(3 * n, sd = sd), sd = rep(sd, n))
# 
# ggplot(df, aes(x)) +
#   geom_histogram(binwidth = histogram_bins("FD")) +
#   facet_wrap(~ sd, scales = "free_x") +
#   labs(x = NULL)

