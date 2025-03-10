#!/usr/bin/env Rscript

arg <- commandArgs(T)

# 脚本特点----
# - 仅支持长选项
# - 可以轻松设置默认值
# - 长选项后的参数只允许用`=` 链接，不允许接空格分隔的参数。
# - 参数`允许多次指定`,即可写`--input=a.txt  --input=b.txt`，如果对于有参数的选项，如果只允许一个值的话，注意核查返回结果长度。
# - 可以设置多类型的选项`character`（默认）、`integer`、`numeric`、`logical`等等。这里是使用`as()`函数进行实现，该函数在流程中比较少见，但是不难，一看就会。
# - 可以设置逻辑值选项，即不需要接参数的选项（我不知道别人怎么称呼这种不需要接参数的选项的，因为我都是用于转换成逻辑值，从而管理流程，因此我习惯将其称作逻辑值选项）。注意逻辑值不得有`=` 。
# - 帮助文档中脚本名称需要手动改，不会直接获取。
# - 不要尝试空格分割选项和参数的写法了，更麻烦！ 尝试直接通过设置（set_options函数）的选项列表`arg_list`解析帮助文档写法，也很麻烦，也不够直观。
# - 本次只分享了脚本参数解析与选项设置部分，包含模块如下：
# 
# 脚本测试参数（仅作测试时启用）、帮助文档、设置相应选项（长选项、类型、默认值）、检查选项、解析选项值。各模块都是包含函数设定与函数使用，并未将所有函数整理到一起，更方便查错。

if(F){
  arg[1] <- "--input=b"
  #arg[2] <- "hi"
  
  arg[2] <- "--output=result.txt" # 超出边界为NA
  arg[3] <- "--threads=5" # 超出边界为NA
  arg[4] <- "--input=a"
  
  arg[5] <- "--force"
  # 正常写法
  arg <- c("--input=b.txt","--output=result.txt","--input=a.txt","--threads=5", "--force")
  
  # 各种查错
  # arg <- c("--input=","--hi","--threads=5.2", "--force=") #报错：Error: Option: --hi cannot be recognized!
  
  ## 删 "--hi"
  # arg <- c("--input=","--threads=5.2", "--force=")  #报错：  Error: Option: --force= cannot be recognized!
  
  # 删"--force=" 
  # arg <- c("--input=","--threads=5.2") #报错：Error in get_optins(arg = arg, arg_list_option = arg_list[[i]]) : Option: --input=<INPUT> must have a effective value!
  
  # 改"--force=input.txt" 
  # arg <- c("--input=input.txt","--threads=5.2") # threads 会被强制转为数值5(使用integer())
  
  # 改"--threads=no" 非数值
  # arg <- c("--input=input.txt","--threads=no") # Warning message:   In asMethod(object) : NAs introduced by coercion # 即无法转为数值，返回NA。
  # 这个警告提示不太方便观察，改的话代码需要单独写，麻烦！没改！
}


# 帮助文档----
{ 
  USAGE <- function(){
    message("USAGE:    Rscript_name.R --input=<INPUT> --output=<OUTPUT>\n")
    message("Function: The function of this script!\n")
    message("     --input=<INPUT>  The input file! [required]\n") # 允许多次指定，后面代码是否允许多个值输入？注意配套。
    message("     --output=<OUTPUT>  The output file! [default: ./test.txt]\n")
    message("     --threads=<Threads>  An integer, the threads file![default:5]\n")
    message("     --force  Rerun this script even if the result exits!\n") # 添加该选项的时候执行，则为TRUE
    message("     --help  Display usage\n") 
    message("Author: Shimao Zheng\n")
    message("Email: zhengshimao007@163.com\n")
  }
  USAGE()
}

## 设置相应选项与参数 ----
### 逻辑值选项：默认值为default， 即无该选项时的逻辑值。
{
  set_options <- function(arg_list = arg_list, long, type = "character" ,default = FALSE ){ # 默认参数为字符串类型，
    if(missing(long)) USAGE() & stop("The long option must exits!")
    #arg_list[[gsub("=","",long)]] <- c(long = long, type = type, default = default)
    arg_list[[long]] <- c(long = long, type = type, default = default)
    return(arg_list)
  }
  arg_list <- list()
  # 设置选项
  arg_list <- set_options(arg_list = arg_list, long = "--input=" , type = "character", default = NA) # 可以用NA或这NULL定义默认值，但是NA值最好。
  
  arg_list <- set_options(arg_list = arg_list, long = "--output=", type = "character", default = "./test.txt")
  
  arg_list <- set_options(arg_list = arg_list, long = "--threads=", type = "integer", default = "5")
  
  arg_list <- set_options(arg_list = arg_list, long = "--force", type = "logical", default = FALSE) # 注意type与default的统一，后面没写这个检查！
  
  arg_list <- set_options(arg_list = arg_list, long = "--help", type = "", default = NA) # 只检查长选项，如果选项包含在`--help`,是无法运行后type和default检查部分的，所以可以随便写。
  
  # arg_list # 查看
}

# 检查选项----
{
  # 检查是否有arg_list中未设置的选项出现在arg中
  
  if(any(grepl("--help", arg)))  USAGE()&rm(list=ls()) & q(save = "no")
  
  for (i in arg) {
    op <- gsub("=.*","=",i)
    if(!any(grepl(op, names(arg_list)))){
      USAGE() & stop("Option: ",op," cannot be recognized!")
    }
  }
}  

# 解析选项值----
{
  #arg_list_option = arg_list$`--input=`
  get_optins <- function(arg = arg, # 所有输入选项
                         arg_list_option # 其中的一个arg_list中的一个元素
  ){
    
    if(!grepl("^--\\w+",arg_list_option["long"])) USAGE() & stop("The option you set must match '--\\w+'")
    
    long_option <- grep(pattern = arg_list_option["long"], x= arg, value = TRUE) # 选项之间写法不得太相近，尤其是逻辑值
    # long_option <- gsub(arg_list_option["long"],"",long_option)
    
    # 逻辑值选项
    if(arg_list_option["type"] == "logical"){ 
      if(length(long_option) == 0){
        #option_value <- as(arg_list_option["default"],arg_list_option["type"])
        #return(option_value)
        return(as(arg_list_option["default"],arg_list_option["type"]))
      }else{
        #option_value <- !as(arg_list_option["default"],arg_list_option["type"])
        #return(option_value)
        return(!as(arg_list_option["default"],arg_list_option["type"]))
      }
      #return(option_value) # "logical"
      # 必须加参数选项  
    }else{ 
      ## 当未设置该选项时
      if(length(long_option) == 0){ 
        if(is.na(arg_list_option["default"])){ # 无默认值
          stop("Option: ",arg_list_option["long"],"<",toupper(gsub("[-=]","",arg_list_option["long"])),"> must have a effective value!")
        }else{ # 有默认值
          option_value <- as(arg_list_option["default"], arg_list_option["type"])
          names(option_value) <- NULL
          return(option_value)
          #(as(arg_list_option["default"], arg_list_option["type"]))
        }
        ## 当有设置该选项时
      }else{ 
        option_value <- gsub(arg_list_option["long"],"",long_option)
        option_value <- option_value[!option_value %in% ""]
        if(length(option_value) == 0){ # 设置选项，但是没设置值
          stop("Option: ",arg_list_option["long"],"<",toupper(gsub("[-=]","",arg_list_option["long"])),"> must have a effective value!")
        }else{
          # option_value <- as(option_value, arg_list_option["type"])
          # names(option_value) <- NULL
          # return(option_value)
          return(as(option_value, arg_list_option["type"]))
        }
        
      }
      
    }
    
  }
  
  # 解析参数
  parse_arg_list <- function(arg = arg, arg_list = arg_list){
    args <- list()
    for (i in names(arg_list)) {
      #args[[gsub("-", "",arg_list[[i]][["long"]])]] <- as(object = arg_list[[i]]["default"], Class = "logical" )
      args[[gsub("[-=]", "",arg_list[[i]]["long"] )]] <- get_optins(arg = arg, arg_list_option = arg_list[[i]])
      
    }
    return(args) # 返回一个list类型。
  }
  
  all_args <- parse_arg_list(arg = arg, arg_list = arg_list)
  
  # 测试查看：使用时最后要注销掉！
  all_args
}


# 开始正式写代码----
input <- all_args$input
output <- all_args$output
threads <- all_args$threads
force <- all_args$force