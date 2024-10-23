#!/bin/bash/
# **************************************************
# 下载SRA数据
# Date   : 2023-10-17
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo "descriptive: Download SRA data"
    echo "usage: bash Prefetch_down.sh [Sraids.txt] [Output flodle]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

#  获取当前路径
path=$(pwd)

#  SRA数据下载
if [ -a $path/$2 ] ;then
    rm -rf $path/$2
fi

mkdir $path/$2 && cd $path/$2

echo -e  " \n 开始下载 \n "

for ids in `cat ../$1` ;do prefetch -O ./ $ids ;done

echo "`date '+[哦了] %D %T'`"

echo -e " \n 下载完了 \n "