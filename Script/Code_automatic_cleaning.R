

rm(list = ls())
path <- '/data/nas1/liky/project/YQ547-12/'  #修改成自己的路径

#复制脚本---------------
setwd(path)
if (!dir.exists('R')) {
  dir.create('R')
}
getwd()

file_list <- list.files(path = './', pattern = "\\.r$", full.names = TRUE)

for (file_name in file_list) {
  base_name <- basename(file_name)
  
  # 检查文件名长度是否足够长，第 3 和第 4 位是否为数字，以及是否不包含 "unused"
  if (nchar(base_name) >= 4 &&
      grepl("^r\\.\\d{2}_.*\\.r$", base_name) &&
      !grepl("unused", base_name)) {
    # 复制文件到 "R" 文件夹中
    file.copy(from = file_name, to = './R/')
    message("Copied file: ", base_name)
  }
}




#删除批注----------------
setwd(path)
if (!dir.exists('R')) {
  dir.create('R')
}
setwd(paste0(path,'R'))
getwd()

file_list <- list.files(path = './')

for (file_name in file_list) {
  script_content <- readLines(file_name, warn = FALSE)
  
  for (j in 1:length(script_content)) {
    line <- script_content[j]
    clean_line <- ""
    inside_quotes <- FALSE
    
    for (i in 1:nchar(line)) {
      char <- substr(line, i, i)
      
      if (char %in% c('"', "'")) {
        inside_quotes <- !inside_quotes
      }
      
      # 如果遇到#，且不在引号内，则删除#及后续内容
      if (char == "#" && !inside_quotes) {
        break
      }
      
      clean_line <- paste0(clean_line, char)
    }
    
    script_content[j] <- clean_line
  }
  
  # 删除两行以上的连续空行
  cleaned_content <- c()
  empty_line_count <- 0
  
  for (line in script_content) {
    if (line == "") {
      empty_line_count <- empty_line_count + 1
    } else {
      empty_line_count <- 0
    }
    
    if (empty_line_count <= 1) {
      cleaned_content <- c(cleaned_content, line)
    }
  }
  
  # 将修改后的内容写回原文件
  writeLines(cleaned_content, file_name)
  
  message("Processed file: ", file_name)
}


#删除路径--------------------
setwd(path)
if (!dir.exists('R')) {
  dir.create('R')
}
setwd(paste0(path, 'R'))
getwd()

file_list <- list.files(path = './', pattern = "\\.r$", full.names = TRUE)

replace_text <- function(line, path_to_remove) {
  line <- gsub(path, "", line)
  line <- gsub('/data/nas1/huxiaohua/project/10_YQ547_12/', "", line)
  line <- gsub("/data/nas1/luchunlin", "", line)
  line <- gsub("/data/nas1/tancj", "", line)
  line <- gsub("/data/nas3/xieyunpeng/", "", line)
  line <- gsub("/data/nas1/xieyunpeng/", "", line)
  line <- gsub("/data/nas1/huxiaohua", "", line) #修改成自己的路径
  return(line)
}

for (file_name in file_list) {
  script_content <- readLines(file_name, warn = FALSE)
  processed_content <- sapply(script_content, replace_text, path_to_remove = path, USE.NAMES = FALSE)
  
  writeLines(processed_content, file_name)
  
  message("Processed file: ", basename(file_name))
}
