#!/bin/bash/
# **************************************************
# 文件处理
# Date   : 2023-12-14
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then  #  如果string为空，则为真
    echo -e " \n `date '+[啥也没输入] %D %T'` \n "
    echo "本脚本用于转录组fq文件拼接"
    echo "使用: veen_separate_mod.sh [Input file] [Output file]"
    exit 1
fi

#  获取当前路径
path=$(pwd)

sed 's/\r/,/g' $path/$1 | xargs > $path/$2

sed -i 's/\ /"/g' $path/$2

sed -i 's/^/"/' $path/$2

sed -i s'/.$//' $path/$2

echo -e " \n ♕臭小子，干完了！ \n "