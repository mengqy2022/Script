#!/bin/bash/
# **************************************************
# 将序列和开放阅读框整理，输出列表，转换列表
# Date   : 2024-02-12
# Author : 孟庆瑶
# **************************************************

set -e

name=$(basename $0) 

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo -e "将序列和开放阅读框整理，输出列表，转换列表。 \n"
    echo -e "usage: bash $name [BLAST] > [OUTFILE]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

python3 /mnt/d/Script_mqy/General_tools/shift_code_gene.py $1 > temp_bash

#  处理空行的命令 sed '/^\t*$/d'
cat temp_bash | paste -s | sed 's/>/\n/g' | sed '1d' | awk -F'\t' '{if($3!="") print $0}'

rm -rf temp_bash