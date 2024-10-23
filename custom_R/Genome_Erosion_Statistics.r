GES <- function(data) {
  for (o in 1:length(rownames(data))) {
    if (is.na(data[o,2]) & is.na(data[o,3]) &
        is.na(data[o,4]) & is.na(data[o,5]) &
        is.na(data[o,6])) {
      data[o,7] <- 5
    } else if (!is.na(data[o,2]) & is.na(data[o,3]) &
               is.na(data[o,4]) & is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 4
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               is.na(data[o,4]) & is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 4
    } else if (is.na(data[o,2]) & is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 4
    } else if (is.na(data[o,2]) & is.na(data[o,3]) &
               is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 4
    } else if (is.na(data[o,2]) & is.na(data[o,3]) &
               is.na(data[o,4]) & is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 4
    } else if (is.na(data[o,2]) & is.na(data[o,3]) &
               is.na(data[o,4]) & !is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (is.na(data[o,2]) & is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 3
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
              !is.na(data[o,4]) & is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 3
    } else if (!is.na(data[o,2]) & is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 3
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 3
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               is.na(data[o,4]) & is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 3
    } else if (!is.na(data[o,2]) & is.na(data[o,3]) &
               is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 3
    } else if (!is.na(data[o,2]) & !is.na(data[o,3]) &
               is.na(data[o,4]) & is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 3
    } else if (!is.na(data[o,2]) & is.na(data[o,3]) &
               is.na(data[o,4]) & !is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (!is.na(data[o,2]) & is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (!is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (is.na(data[o,2]) & is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               is.na(data[o,4]) & !is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (!is.na(data[o,2]) & !is.na(data[o,3]) &
               is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (!is.na(data[o,2]) & !is.na(data[o,3]) &
               is.na(data[o,4]) & is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (!is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (!is.na(data[o,2]) & is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 2
    } else if (!is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 1
    } else if (!is.na(data[o,2]) & is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 1
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 1
    } else if (is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 1
    } else if (!is.na(data[o,2]) & !is.na(data[o,3]) &
               is.na(data[o,4]) & !is.na(data[o,5]) &
               !is.na(data[o,6])) {
      data[o,7] <- 1
    } else if (!is.na(data[o,2]) & !is.na(data[o,3]) &
               !is.na(data[o,4]) & !is.na(data[o,5]) &
               is.na(data[o,6])) {
      data[o,7] <- 1
    }
  }
  return(data)
}