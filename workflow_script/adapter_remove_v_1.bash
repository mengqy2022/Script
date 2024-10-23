#!/bin/bash/
# **************************************************
# 去接头
# Date   : 2024-01-16
# Author : 孟庆瑶
# **************************************************

et -e

name=$(basename $0) 

if [ -z $1 ] ;then  #  如果string为空，则为真
    echo -e " \n `date '+[啥也没输入] %D %T'` \n "
    echo -e "去除接头和低质量过滤[trimmomatic] \n"
    echo -e "需要按照自己需求更改接头序列类型\n"
    echo -e "使用: bash $name [DATABASE] [OUTPUT FPREFILE]"
    echo -e "\n             >>>> mqy <<<< \n "
    echo -e "  文件输入顺序不能乱 \n "
    exit 1
fi

#  获取当前路径
path=$(pwd)

if [ -a $path/$2 ] ;then
    rm -rf $path/$2
fi

echo -e " \n 开始去接头，过滤！ \n "

mkdir $2

cd $1

for i in `ls *gz | cut -d "_" -f1 | sort -u`
do 
    java -jar /home/mengqingyao/biosoftware/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $path/$1/${i}_1*.gz $path/$1/${i}_2*.gz \
        -baseout $path/$2/${i}_trim.fq.gz ILLUMINACLIP:/home/mengqingyao/biosoftware/Trimmomatic-0.39/adapters/gene_denovo.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:20
done

echo -e " \n 干完了！ \n "