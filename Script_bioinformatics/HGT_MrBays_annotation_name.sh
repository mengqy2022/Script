#!/bin/bash/
# **************************************************
# 根换进化树的名称
# Date   : 2024-03-13
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ; then
    echo -e "`date '+[error] %D %T'`\n"
    filename=$(basename $0)
    echo -e "描述：${filename} 脚本用于替换进化树的名称。[MrBays]\n"
    echo -e "usage: bash ${filename} [FILE NAME]\n"
    echo -e "example: bash ${filename} MrBays_annotation.fasta"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

pwd=$(pwd)

file=$1

if [! -f $file ] ; then
    echo -e "`date '+[error] %D %T'`\n"
    echo -e "File not found: $file \n"
    exit 1    
fi  

echo -e " \n >>>> Start <<<< \n "    
echo -e "File: $file"

cat ${file} | grep ">" | cut -d ":" -f 1 | sed 's/>//g' | sed 's/\./_/g' | sed 's/\-/_/g' | cut -d " " -f 1 > 1.ids

#  想要替换的名称
cat 1.ids | cut -d "_" -f 1,2,3 > 2.ids

if [ -f ann_name.txt ] ;then
    rm -rf ann_name.txt
fi

paste 1.ids 2.ids > ann_name.txt && rm 1.ids 2.ids

awk 'BEGIN {print "name\trename"; while (getline < "ann_name.txt") print}' > 123 && mv 123 ann_name.txt

echo -e " \n >>>> Done <<<< \n "
echo -e "Results are stored in $pwd \n"
echo -e "Thank you for using our software! \n"
echo -e "Please contact us if you have any questions. \n"
echo -e "Email: <15877464851@163.com> \n"
echo -e "Sincerely, \n"
echo -e "MQY \n"