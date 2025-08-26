#!/usr/bin/env Rscript

library(argparse)

# 主函数
main <- function() {
  # 创建参数解析器
  parser <- ArgumentParser(description = "R脚本批处理工具：复制文件、删除注释、清理路径")
  
  # 添加参数
  parser$add_argument("--path", type = "character", required = TRUE,
                      help = "【必需】原始R脚本所在目录的绝对路径（示例：/data/user/project）")
  parser$add_argument("--output", type = "character", required = TRUE,
                      help = "【必需】输出目录名称（示例：processed_scripts）")
  parser$add_argument("--steps", type = "character", nargs = "+",
                      choices = c("copy", "comments", "paths", "all"),
                      default = "all",
                      help = "【可选】处理步骤选择：copy=复制脚本, comments=删除注释, paths=清理路径, all=全部执行（默认）")
  
  # 解析参数
  args <- parser$parse_args()
  
  # 检查路径是否存在
  if (!dir.exists(args$path)) {
    stop("【错误】指定路径不存在：", args$path)
  }
  
  message("\n========== 脚本批处理开始 ==========")
  message("原始路径：", args$path)
  message("输出目录：", file.path(args$path, args$output))
  
  # 执行选定的步骤
  if ("all" %in% args$steps) {
    args$steps <- c("copy", "comments", "paths")
  }
  
  # 步骤1：复制脚本
  if ("copy" %in% args$steps) {
    message("\n>>> [步骤1/3] 正在复制脚本...")
    copy_scripts(args$path, args$output)
  }
  
  # 步骤2：删除注释
  if ("comments" %in% args$steps) {
    message("\n>>> [步骤2/3] 正在删除注释...")
    remove_comments(args$path, args$output)
  }
  
  # 步骤3：清理路径
  if ("paths" %in% args$steps) {
    message("\n>>> [步骤3/3] 正在清理路径...")
    remove_paths(args$path, args$output)
  }
  
  message("\n========== 处理完成 ==========")
  message("输出文件位于：", file.path(args$path, args$output))
}

# 函数1：复制符合条件的R脚本
copy_scripts <- function(base_path, output_dir) {
  setwd(base_path)
  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
    message("创建输出目录：", output_dir)
  }
  
  # 获取所有R脚本文件
  file_list <- list.files(path = './', pattern = "\\.(r|R)$", full.names = TRUE)
  copied_files <- 0
  
  for (file_name in file_list) {
    base_name <- basename(file_name)
    
    # 检查文件名格式：r.数字_xxx.R
    if (nchar(base_name) >= 4 &&
        grepl("^r\\.\\d{2}_.*\\.(r|R)$", base_name) &&
        !grepl("unused", base_name, ignore.case = TRUE)) {
      file.copy(from = file_name, to = file.path(output_dir, base_name))
      message("已复制：", base_name)
      copied_files <- copied_files + 1
    }
  }
  
  message("共复制 ", copied_files, " 个符合条件的脚本文件")
}

# 函数2：删除脚本注释
remove_comments <- function(base_path, output_dir) {
  setwd(base_path)
  setwd(output_dir)
  
  file_list <- list.files(path = './', pattern = "\\.(r|R)$")
  processed_files <- 0
  
  for (file_name in file_list) {
    script_content <- readLines(file_name, warn = FALSE)
    cleaned_content <- character(length(script_content))
    
    # 逐行处理注释
    for (j in seq_along(script_content)) {
      line <- script_content[j]
      clean_line <- ""
      inside_quotes <- FALSE
      quote_char <- ""
      
      # 逐个字符分析
      for (i in seq_len(nchar(line))) {
        char <- substr(line, i, i)
        
        # 处理引号内内容
        if (char %in% c('"', "'")) {
          if (!inside_quotes) {
            inside_quotes <- TRUE
            quote_char <- char
          } else if (char == quote_char) {
            inside_quotes <- FALSE
          }
        }
        
        # 遇到注释符号且不在引号内时终止
        if (char == "#" && !inside_quotes) {
          break
        }
        
        clean_line <- paste0(clean_line, char)
      }
      
      cleaned_content[j] <- clean_line
    }
    
    # 移除多余空行（保留单空行）
    final_content <- cleaned_content[nzchar(cleaned_content) | 
                                       c(TRUE, !(nzchar(cleaned_content)[-length(cleaned_content)] & 
                                                   nzchar(cleaned_content)[-1]))]
    
    writeLines(final_content, file_name)
    processed_files <- processed_files + 1
    }
  
  message("已处理 ", processed_files, " 个文件的注释删除")
  }

# 函数3：清理路径字符串
remove_paths <- function(base_path, output_dir) {
  setwd(base_path)
  setwd(output_dir)
  
  file_list <- list.files(path = './', pattern = "\\.(r|R)$")
  processed_files <- 0
  
  for (file_name in file_list) {
    script_content <- readLines(file_name, warn = FALSE)
    
    # 替换所有路径字符串
    processed_content <- gsub(base_path, "", script_content, fixed = TRUE)
    
    # 额外清理常见路径模式
    processed_content <- gsub("/+", "/", processed_content)  # 清理多余斜杠
    processed_content <- gsub('^["\']/', "", processed_content)  # 清理引号开头的路径
    
    writeLines(processed_content, file_name)
    processed_files <- processed_files + 1
  }
  
  message("已清理 ", processed_files, " 个文件的路径信息")
}

# 执行主函数（非交互模式时）
if (!interactive()) {
  main()
}
