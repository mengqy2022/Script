#  根据箱线图去除异常值
remove_outliers_IQR <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# 这么一来谁是异常值一看就清楚了，如果变量一般近似正态分布，可以使用拉依达准则（3σ准则）把超过偏离均数三倍标准差的值认定为异常值：
remove_outliers_1 <- function(x) {
  
  lower_limit <- mean(dt)-3*sd(dt) #计算下限
  Upper_limit <- mean(dt)+3*sd(dt) #计算上限
  
  y <- x
  y[x > Upper_limit ] <- NA
  y[x < lower_limit] <- NA
  y
  
}

#  盖帽法
remove_outliers_2 <- function(x) {
  
  u=quantile(dt,0.99)
  l=quantile(dt,0.01)
  
  y <- x
  
  y[x < l] <- l
  y[x > u] <- u
  y
  
}

#  Hampel滤波
remove_outliers_2 <- function(x) {
  
  u = median(dt)-3*mad(dt) #中位数减去3倍的绝对中位差
  l = median(dt)+3*mad(dt) #中位数加上3倍的绝对中位差

  y <- x
  
  y[x > l] <- NA
  y
  
}

#  归一化
normalize <-function(x){
  (x-min(x))/(max(x)-min(x))
}

# 识别异常值
Identify_outliers <- function(x, method = "boxplot", k = NULL,
                          coef = NULL, lp = NULL, up = NULL) {
  switch(method,
         "sd" = {
           if (is.null(k)) k = 3
           mu = mean(x, na.rm = T)
           sd = sd(x, na.rm = T)
           LL = mu - k * sd
           UL = mu + k * sd
         },
         "boxplot" = {
           if (is.null(coef)) coef = 1.5
           Q1 = quantile(x, 0.25, na.rm = T)
           Q3 = quantile(x, 0.75, na.rm = T)
           iqr = Q3 - Q1
           LL = Q1 - coef * iqr
           UL = Q3 + coef * iqr
         },
         "percentiles" = {
           if (is.null(lp)) lp = 0.025
           if (is.null(up)) up = 0.975
           LL = quantile(x, lp, na.rm = T)
           UL = quantile(x, up, na.rm = T)
         })
  idx = which(x < LL | x > UL)
  n = length(idx)
  list(outliers = x[idx], outliers_idx = idx, outliers_num = n)
}
