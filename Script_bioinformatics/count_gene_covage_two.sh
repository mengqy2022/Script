#!/bin/bash/
# **************************************************
# 用来替换文件名和序列名称
# Date   : 2024-03-13
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ; then
    echo -e "`date '+[error] %D %T'`\n"
    filename=$(basename $0)
    echo -e "描述：${filename} 用来替换文件名和序列名称\n"
    echo -e "usage: bash ${filename} [Genomics FOLDER]\n"
    echo -e "Genomics FOLDER: the folder contains the genomics data, including the fasta file.\n"
    echo -e "example: bash ${filename} /data/database/foldername"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

pwd=$(pwd)
echo "Current directory: $pwd"

if [ ! -d $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo -e "[Genomics FOLDER: $1 does not exist.]\n"
    echo -e "Please check the folder name and try again. \n"
fi

echo -e " \n >>>> Start <<<< \n "    
echo -e "PATH: $pwd/$1"

cd $pwd/$1

head -n 1 * | sed 's/==> //g' | sed  's/<==//g' | sed 's/\.fna//g' | sed 's/\.fasta//g' | sed 's/\.fa//g' | sed 's/\.fnd//g' | sed 's/\.faa//g' | sed 's/\.ffn//g' | sed 's/>//g' | \
    cut -d " " -f 1 | sed 's/\./_/g' | paste -s -d " " | awk '{for (i = 1; i <= NF; i++) {if(i % 2 == 0) {print $i} else {printf("%s\t", $i)}}; if((i - 1) % 2 != 0) {printf("\n")}}' > $pwd/cds_coverage_rename.txt

cd $pwd

awk 'BEGIN{print "name rename"; while (getline < "cds_coverage_rename.txt") print}' > 123 && mv 123 cds_coverage_rename.txt

echo -e " \n >>>> Done <<<< \n "
echo -e "Results is $pwd/$2.txt \n"
echo -e "Thank you for using our software! \n"
echo -e "Please contact us if you have any questions. \n"
echo -e "Email: <15877464851@163.com> \n"
echo -e "Sincerely, \n"
echo -e "MQY \n"