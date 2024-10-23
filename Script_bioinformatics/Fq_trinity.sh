#!/bin/bash/
# **************************************************
# 转录组fq文件拼接
# Date   : 2023-11-1
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then  #  如果string为空，则为真
    echo -e " \n `date '+[啥也没输入] %D %T'` \n "
    echo "本脚本用于转录组fq文件拼接"
    echo "使用: bash Fq_trinity.sh [Input flodle]"
    echo -e " \n 如果数量太多不要使用本脚本 \n "
    echo -e "             >>>> mqy <<<< \n "
    exit 1
fi

#  获取当前路径
path=$(pwd)

if [ -a $path/Fq_trinity ] ;then
    rm -rf $path/Fq_trinity
fi

mkdir $path/Fq_trinity && cd $path/Fq_trinity

while read sample fq1 fq2 ;do
    echo -e " \n 开始运行 \n "
    source /home/mengqingyao/miniconda3/bin/activate trinity_2.3.2
    Trinity --seqType fq --max_memory 10G --left ../$fq1 --right ../$fq2 --CPU 6 --output $sample &
done
echo -e " \n 干完了！ \n "