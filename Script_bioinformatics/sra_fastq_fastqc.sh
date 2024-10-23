#!/bin/bash/
# **************************************************
# 下载SRA数据
# 获得fastq序并进行质量检测
# Date   : 2023-10-17
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo "descriptive: Download SRA data, convert to fastq, quality assessment"
    echo "usage: bash General_no1.sh [Sraids.txt] [Output flodle]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

#  获取当前路径
path=$(pwd)

#  SRA数据下载
if [ -a $path/$2 ] ;then
    rm -rf $path/$2
fi

mkdir $path/$2

echo -e  " \n 开始下载 \n "

for ids in `cat $1` ;do prefetch -O $path/$2 $ids ;done

echo "`date '+[哦了] %D %T'`"

echo -e " \n 下载完了 \n "

#  -e filename 如果 filename存在，则为真 [ -e /var/log/syslog ]
if [ -e $path/$2 ] ;then
    cd $path/$2
else
    echo "文件夹不存在"
    exit 1
fi

echo  -e "\n  处理SRA文件 !!! \n "

for ids in `cat ../$1` ;do cd $ids && mv $ids* ../ && cd ../ && rm -rf ./$ids ;done

#  目前还在SRA_database中

for ids in `cat ../$1` ;do echo "获得fastq文件中" && fasterq-dump -3 -e 12 -O ./ $ids*  && gzip *.fastq ;done

echo "`date '+[哦了] %D %T'`"

echo -e " \n  开始质量检测  \n " 

mkdir ./fastqc && cd ./fastqc

ls ../*.gz | xargs fastqc -t 12 -o ./

multiqc ./

cd ../

for ids in `cat ../$1` ;do gzip *.sra ;done

echo "`date '+[哦了] %D %T'`"
#  -e：激活转义字符。使用-e选项时，若字符串中出现以下字符，则特别加以处理，而不会将它当成一般文字输出：
echo -e  " \n 都干完了 \n "