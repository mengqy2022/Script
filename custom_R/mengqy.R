calculate_mean <- function(x,w){
  sum(x * w) / sum(x)
}

calculate_var <- function(x,w){
  mu <- calculate_mean(x, w)
  sum(w * (x_mu) ^ 2 / sum(w))
}

calculate_sd <- function(x, w){
  sqrt(calculate_varx,w)
}

calculate_mean_same <- function(x,w){
  if (length(x) != length(w)) {
    stop("`x` and `w` must be the same length", call. = F)
  }
  sum(x * w) / sum(x)
}

calculate_mean_same_mod <- function(x, w, na.rm = FALSE){
  
  if (!is.logical(na.rm)) {
    stop("`na.rm` must be logical")
  }
  
  if (length(na.rm) != 1) {
    stop("`na.rm` must be length 1")
  }
  
  if (length(x) != length(w)) {
    stop("`x` and `w` must be the same length", call. = F)
  }
  
  if (na.rm) {
    miss <- is.na(x) | is.na(w)
    x <- x[!miss]
    y <- y[!miss]    
  }
  sum(x * w) / sum(x)
}

calculate_mean_stopifont <- function(x, w, na.rm = FALSE){
  stopifnot(is,logical(na.rm), length(na.rm) == 1)
  stopifnot(length(x) == length(w))
  
  if (na.rm) {
    miss <- is.na(x) | is.na(w)
    x <- x[!miss]
    y <- y[!miss]    
  }
  sum(x * w) / sum(x)
}

commas <- function(..., collapse = ",") {
  stringr::str_c(...,collapse = collapse)
}

rule <- function(..., pad = "-") {
  title <- paste0(...)
  width <- getOption("width") - nchar(title) -5 # nchar("mqy") 3
  cat(title, " ", stringr::str_dup(pad, width), "\n", sep = "")
}

show_missings <- function(df){
  n <- sum(is.na(df))
  cat("Missing values: ", n, "\n", sep = "")
  
  #return(df)
  invisible(df)
}


col_summary <- function(df, fun) {
  out <- vector("double", length(df))
  for (i in seq_along(df)) {
    out[i] <- fun(df[[i]])
  }
  out
}

scatter_plot <- function(df, x_var, y_var) {
  x_var = enquo(x_var)
  y_var = enquo(y_var)
  ggplot(data = df, aes(x = !!x_var, y = !!y_var)) + 
    geom_point() + 
    theme_bw() + 
    theme(plot.title = element_text(lineheight = 1, face = "bold", hjust = 0.5)) + 
    geom_smooth() + 
    ggtitle(str_c(as_label(y_var), " vs. ", as_label(x_var)))
}
