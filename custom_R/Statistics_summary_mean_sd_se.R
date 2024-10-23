summary_mean_sd_se <- function(data, .summary_var, ...) {
  summary_var = enquo(.summary_var)
  data %>% 
    group_by(...) %>% 
    summarise(mean = mean(!!summary_var, na.rm = T),
              sd = sd(!!summary_var, na.rm = T),
              se = sd / sqrt(n()))
}
