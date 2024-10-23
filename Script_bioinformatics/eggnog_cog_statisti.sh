#!/bin/bash/
# **************************************************
# 处理eggnog注释文件
# Date   : 2024-01-04
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo "Modify the sequence name"
    echo -e "usage: bash eggnog_cog_statisti.sh [EGGNOG ANNOTATION]"
    echo -e "\t[The order of documents must not be incorrect.]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

cat $1 | grep -v "##" | grep -v "#" > temp

awk -F "\t" '{if ($7 != "S" && $7 != "-" && $7 != "R") {print $1"\t"$5"\t"$7}}' temp > out.emapper.annotations.COG.txt && rm -rf temp