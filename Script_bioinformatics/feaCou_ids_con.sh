#!/bin/bash/
# **************************************************
# 处理eggnog注释文件
# Date   : 2024-01-04
# Author : 孟庆瑶
# **************************************************

set -e

name=$(basename $0) 

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo -e "eggnog:获得ids之间的对应关系 \n"
    echo -e "usage: bash $name [GTF]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

echo -e "\n 输入的文件为：$* \n"

#提取gtf注释文件中gene_id等与gene_name的对应关系,便于下游id转换
gtf=$1

### gene_id to gene_name
#  -F 字段分隔符
grep 'gene_id' $gtf | awk -F 'gene_id \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_id_tmp
grep 'gene_id' $gtf | awk -F 'gene_name \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_name_tmp
paste gene_id_tmp gene_name_tmp >last_tmp
uniq last_tmp > gene_id_gene_name_gencode.txt
rm *_tmp

### transcript_id to gene_name
grep 'transcript_id' $gtf | awk -F 'transcript_id \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_id_tmp
grep 'transcript_id' $gtf | awk -F 'gene_name \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_name_tmp
paste gene_id_tmp gene_name_tmp >last_tmp
uniq last_tmp > transcript_gene_name_gencode.txt
rm *_tmp

echo -e " \n 老弟，干完了！ \n "