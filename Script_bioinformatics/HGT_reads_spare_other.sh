#!/bin/bash/
# **************************************************
# 将序列与注释结果相对应，分离序列，并更改序列名称
# Date   : 2024-02-21
# Author : 孟庆瑶
# **************************************************

set -e

name=$0

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo "额外加入查询序列，提前将查询序列加入到[blastdbcmd]获得的数据库中，在运行。"
    echo -e "usage: bash $name [SPLIT IDS] [BLAST FOLDER] [FASTA FILE] [OUTPUT FILE]"
    echo -e "\t[The order of documents must not be incorrect.]"
    echo -e "\t[请输入相对路径]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

path=$(pwd)


if [ -a $path/$4 ] ;then
    rm -rf $path/$4
fi

mkdir $path/$4

cd $path/$2

if [ -a temp ] ;then
    rm -rf temp
fi

mkdir temp && cd temp

for i in $(sed '/^$/d' $path/$1);do
    cut -f 2 ../${i}* > ${i}.ids
    cut -f 1 ../${i}* | uniq >> ${i}.ids
    echo -e "\n  >>>> [${i} 完成] <<<< \n"
    seqkit faidx $path/$3 -l ${i}.ids > $path/$4/${i}.fa ;done

rm -rf $path/$2/temp