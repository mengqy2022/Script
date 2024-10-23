#!/bin/bash/
# **************************************************
# 去接头
# Date   : 2024-01-19
# Author : 孟庆瑶
# **************************************************

set -e

#  获取脚本名称并输出
name=$(basename $0) 

if [ -z $1 ] ;then  #  如果string为空，则为真
    echo -e " \n `date '+[啥也没输入] %D %T'` \n "
    echo -e "去除接头和低质量过滤[cutadapt] \n"
    echo -e "需要按照自己需求更改接头序列类型\n"
    echo -e ">>>> 只有与双端序列，单端请参考官方文件 <<<<\n"
    echo -e "使用: bash $name [DATABASE] [OUTPUT FPREFILE] [read1 adapter] [read2 adapter] "
    echo -e "\n             >>>> mqy <<<< \n "
    echo -e "  文件输入顺序不能乱 \n "
    exit 1
fi

#  获取当前路径
path=$(pwd)

if [ -a $path/$2 ] ;then
    rm -rf $path/$2
fi

echo -e "\n The parametes are $*."

echo -e " \n 开始去接头，过滤！ \n "

mkdir $2

cd $1

for i in `ls *gz | cut -d "_" -f1 | sort -u`
do 
    #  第一次运行不设置-M，
    #  根据fastqc的结果，设置-M的值
    #  1百万个序列平均长度27bp，最大147bp，设置40bp
    cutadapt -a $3 -A $4 -q 30,30 -o $path/$2/${i}_1_cuta.fq.gz -p $path/$2/${i}_2_cuta.fq.gz $path/$1/${i}_1*.gz $path/$1/${i}_2*.gz
    seqkit fx2tab -j 8 -l -n -i -H $path/$2/${i}_1_cuta.fq.gz | cut -f 2 > $path/$2/${i}_1_cuta.txt
    seqkit fx2tab -j 8 -l -n -i -H $path/$2/${i}_2_cuta.fq.gz | cut -f 2 > $path/$2/${i}_2_cuta.txt

done

echo -e " \n 干完了！ \n "