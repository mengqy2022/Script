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
OUTPREFIX="基因注释信息"

while getopts hq:k:o:m: opt
do
    case "$opt" in 
    h) 
        echo
        echo "脚本说明: 处理单个物种的KAAS预测文件，(*^▽^*)！"
        echo
        echo "使用说明: bash $name -q query.ko -k q00001.keg -o out_prefix"
        echo "------------------------------------------------------------------------------------"
        echo -e "\t-q: 输入功能注释文件;"
        echo -e "\t-k: KEGG Orthology (KO) Download htext;"
        echo -e "\t-o: 结果文件前缀, 默认:基因注释信息;"
        echo "`date '+Date: %D %T'`"
        echo -e "\n             >>>> mqy <<<< \n "
        echo "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "     !!!!   保证所有文件都在一个工作目录下   !!!!"
        echo "     !!!!          可以处理单个物种          !!!!"
        echo "     !!!! 输出结果中包含:                    !!!!"
        echo "     !!!!               基因注释文件;        !!!!"
        echo "     !!!!               通路注释文件;        !!!!"
        echo "     !!!!               修改的K注释文件;     !!!!"
        echo "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo
        exit 0
        ;;
    q) 
        echo "输入-q文件为: $OPTARG"
        QUERYFILE=$OPTARG
        ;;
    k) 
        echo "输入-k文件为: $OPTARG"
        KOFILE=$OPTARG
        ;;
    o) 
        echo "设置输出文件名前缀为: $OPTARG"
        OUTPREFIX=$OPTARG
        ;;
    *) 
        echo "未知参数: $opt"
        ;;
    esac
done

#####################-z 判断字符串是否为0

if [ -z $QUERYFILE ]
then
    echo "请输入-h，查看帮助文档！"
    exit 1
fi

awk '$1!="" && $2!="" {print $2"\t"$1}' $QUERYFILE > $QUERYFILE\_mod.ko 

python3 /biosoftware/kegg_trans.py ./$KOFILE ./$QUERYFILE\_mod.ko ./$OUTPREFIX\_gene_anno.txt ./$OUTPREFIX\_pathway_anno.txt
