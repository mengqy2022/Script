#!/bin/bash/
# **************************************************
# 将GTF转变成外显子的bed文件
# Date   : 2024-06-04
# Author : 孟庆瑶
# **************************************************

set -e

name=$0

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo -e "\n提取GTF中的启动子"
    echo -e "usage: bash $name [INPUT GTF FILE] [OUTPUT FILE]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

path=$(pwd)

if [ -a $path/$2 ] ;then
    rm -rf $path/$2
fi

cat $path/$1 |grep -v  '^#' | sed 's/; /\t/g' |awk 'BEGIN{FS=OFS="\t"} $3 == "exon" { \
    gsub(/gene_id "/,"",$9); gsub(/"/,"",$9) \
    gsub(/transcript_id "/,"",$11);gsub(/"/,"",$11) \
    gsub(/_number "/,"",$13);gsub(/"/,"",$13) \ 
    print $1, $4-1, $5,$9"_"$11"_"$13,".",$7 \
    }' > $path/$2

echo -e " \n >>>> Done <<<< \n "
echo -e "Results are stored in $path \n"
echo -e "Thank you for using our software! \n"
echo -e "Please contact us if you have any questions. \n"
echo -e "Email: <15877464851@163.com> \n"
echo -e "Sincerely, \n"
echo -e "MQY \n"