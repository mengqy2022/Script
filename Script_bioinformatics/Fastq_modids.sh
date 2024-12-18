#!/bin/bash

# 初始化变量
R1=""
R2=""
USE_ZCAT=false

# 解析命令行参数
while getopts "1:2:zc" opt; do
  case ${opt} in
    1 )
      R1=$OPTARG
      ;;
    2 )
      R2=$OPTARG
      ;;
    z )
      USE_ZCAT=true
      ;;
    c )
      USE_ZCAT=false
      ;;
    \? )
      echo "Usage: cmd [-1 R1_file] [-2 R2_file] [-z (for zcat) | -c (for cat)]"
      exit 1
      ;;
  esac
done

# 检查是否提供了必要的参数
if [[ -z "$R1" || -z "$R2" ]]; then
  echo "Both R1 and R2 files are required."
  exit 1
fi

# 设置命令
if $USE_ZCAT; then
  CMD1="zcat $R1"
  CMD2="zcat $R2"
else
  CMD1="cat $R1"
  CMD2="cat $R2"
fi

# 处理 R1 文件
eval "$CMD1" | perl -ne 'if($_=~/^\@/){$_=~s/ 1:N:.*$/\#\/1/g;print "$_";} else {print "$_";}' > "${R1%.fq.gz}_ids.fq" &

# 处理 R2 文件
eval "$CMD2" | perl -ne 'if($_=~/^\@/){$_=~s/ 2:N:.*$/\#\/2/g;print "$_";} else {print "$_";}' > "${R2%.fq.gz}_ids.fq" &

# 等待所有后台进程完成
wait

echo "Processing completed."
