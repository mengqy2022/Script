#!/bin/bash/
# **************************************************
# 提取GTF中的启动子
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

sed 's/"/\t/g' $path/$1 | awk 'BEGIN{OFS=FS="\t"} {if($3=="gene") {ensn=$10; symbol = $16; \
    if($7 == "+") {start=$4-1; up = stat - 1000; if(up < 0) up = 0; dw = start + 500; print $1, up, dw, ensn, symbol, $7;}\
    else if($7 == "-"){start = $5 - 1; up = start + 1000; dw = $4 - 1; if(dw < 0) dw = 0; print $1, dw, up, ensn, symbol, $7}}}'> $path/$2

echo -e " \n >>>> Done <<<< \n "
echo -e "Results are stored in $path \n"
echo -e "Thank you for using our software! \n"
echo -e "Please contact us if you have any questions. \n"
echo -e "Email: <15877464851@163.com> \n"
echo -e "Sincerely, \n"
echo -e "MQY \n"