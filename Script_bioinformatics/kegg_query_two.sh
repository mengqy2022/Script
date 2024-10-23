#!/bin/bash/
# **************************************************
# 转换KEGG注释格式
# Date   : 2024-07-31
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    filename=$(basename $0)
    echo -e "\nConvert KEGG annotation result format\n" 
    echo -e "usage: bash ${filename} [INPUT FILE] [OUTPUT FILE]"
    echo -e "example: bash ${filename} query.ko query_1.ko"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

pwd=$(pwd)
echo "Current directory: $pwd"

if [ -a $pwd/$2 ] ;then
    rm -rf $pwd/$2
fi

awk '{$2=($2?$2:"-"); print}' $pwd/$1 | grep -v '-'> $pwd/$2

echo -e " \n >>>> Done <<<< \n "
echo -e "Results are saved in $2. \n"
echo -e "Thank you for using our software! \n"
echo -e "Please contact us if you have any questions. \n"
echo -e "Email: <15877464851@163.com> \n"
echo -e "Sincerely, \n"
echo -e "MQY \n"