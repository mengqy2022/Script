#  bootstrap生成器

boot_permute <- function(df, var) {
  n <- nrow(df)
  force(var)
  
  function() {
    col <- df[[var]]
    col[sample(n, replace = T)]
  }
}

# boot_permute1 <- boot_permute(mtcars, 'mpg')
# head(boot_permute1())

boot_model <- function(df, formula) {
  mod <- lm(formula, data = df)
  fitted <- unname(fitted(mod))
  resid <- unname(resid(mod))
  rm(mod)
  
  function() {
    fitted + sample(resid)
  }
}

# boot_mtcars2 <- boot_model(mtcars, mpg ~ wt)
# head(boot_mtcars2())
