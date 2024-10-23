# 2024/01/01
# mqy

first_obse <- function(data,n = 1, sep=",") {
  list <- as.list(data)
  #  保留第一个观测值
  for (v in 2:length(colnames(data))) {
    for (o in 1:length(rownames(data))) {
      list[[v]][o] <- list[[v]][[o]] %>% strsplit(sep)
      list[[v]][o] <- list[[v]][[o]][n]
    }
  }

  #  转变字符串
  for (v in 2:length(colnames(data))) {
    list[[v]] <- as.character(list[[v]])
  }

  #  转换成矩阵
  data.frame <- as.data.frame(list)
  return(data.frame )
}