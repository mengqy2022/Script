#!/bin/bash/
# **************************************************
# 获得fastq序并进行质量检测
# Date   : 2023-10-17
# Author : 孟庆瑶
# **************************************************

#set -e


if [ -z $1 ] ; then
    echo "`date '+[error] %D %T'`"
    echo "descriptive: SRA convert to fastq and quality assessment"
    echo "usage: bash Sra2fq_qc1.sh [SRA database]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

#  获取当前路径
path=$(pwd)

if [ -e $path/$1 ] ; then
    cd $path/$1
else
    echo "文件夹不存在"
    exit 1
fi

ls -lh | cut -d " " -f 9 | sed '1d'>sra.ids.txt

echo  -e "\n  处理SRA文件 !!! \n  " 

for ids in `cat sra.ids.txt`; do cd $ids && mv $ids* ../ && cd ../ && rm -rf ./$ids ;done

#  目前还在SRA_database中

for ids in `cat sra.ids.txt` ;do echo -e " \n 获得fastq文件中 \n " && fasterq-dump -3 -e 12 -O ./ $ids*  && gzip *.fastq ;done

rm -rf sra.ids.txt

echo "`date '+[哦了] %D %T'`"

echo -e " \n  开始质量检测  \n " 

mkdir ./fastqc && cd ./fastqc

ls ../*.gz | xargs fastqc -t 12 -o ./

source /home/mengqingyao/miniconda3/bin/activate

multiqc ./

conda deactivate

cd ../

gzip *.sra &

echo "`date '+[哦了] %D %T'`"
#  -e：激活转义字符。使用-e选项时，若字符串中出现以下字符，则特别加以处理，而不会将它当成一般文字输出：
echo -e  " \n 都干完了 \n "
