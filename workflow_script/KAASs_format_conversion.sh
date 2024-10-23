#!/bin/bash/
# **************************************************
# 处理KAAS预测文件
# Date   : 2023-09-26
# Author : 孟庆瑶
# **************************************************

#  发生错误停止脚本
set -e

#  获取脚本名称并输出
name=$(basename $0) 
echo 
echo "   脚本名称: $name"

#  设置参数以及处理
while getopts hi:d:o: opt
do
    case "$opt" in 
    h) 
        echo
        echo "脚本说明: 处理多个物种的KAAS预测文件，(*^▽^*)！"
        echo
        echo "使用说明: bash $name -i file_name.txt -d folder_prefix -o out_prefix"
        echo "------------------------------------------------------------------------------------"
        echo -e "\t-i: 输入文件名称前缀，文件为单列，每一列是多个物种的文件名称;"
        echo -e "\t-d: 存放物种文件夹的上一级文件夹;"
        echo "`date '+Date: %D %T'`"
        echo -e "\n             >>>> mqy <<<< \n "
        echo "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "     !!!!   保证所有文件都在一个工作目录下   !!!!"
        echo "     !!!!          可以处理多个物种          !!!!"
        echo "     !!!!每个物种有唯一的文件夹名称，名称与-i!!!!"
        echo "     !!!!文件种名称相同，全部存储在-d文件夹下!!!!"
        echo "     !!!!，输出结果在各自物种文件夹下。      !!!!"        
        echo "     !!!! 输出结果中包含:                    !!!!"
        echo "     !!!!               基因注释文件;        !!!!"
        echo "     !!!!               通路注释文件;        !!!!"
        echo "     !!!!               修改的K注释文件;     !!!!"
        echo "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo
        exit 0
        ;;
    i) 
        echo "输入-i文件为: $OPTARG"
        FILENAME=$OPTARG
        ;;
    d) 
        echo "输入-d文件为: $OPTARG"
        FOLDNAME=$OPTARG
        ;;
    *) 
        echo "未知参数: $opt"
        ;;
    esac
done

#####################-z 判断字符串是否为0

if [ -z $FILENAME ]
then
    echo "请输入-h，查看帮助文档！"
    exit 1
fi

for i in `cat ./$FILENAME`;do \
 awk '$1!="" && $2!="" {print $2"\t"$1}' ./$FOLDNAME/${i}/query.ko > ./$FOLDNAME/${i}/query_1.ko | \
 python3 /biosoftware/kegg_trans.py \
 ./$FOLDNAME/${i}/q00001.keg ./$FOLDNAME/${i}/query_1.ko ./$FOLDNAME/${i}/${i}_gene_anno.txt ./$FOLDNAME/${i}/${i}_pathway_anno.txt \
;done