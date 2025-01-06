#!/bin/bash/
# **************************************************
# 将不同文件夹中，文件名称一样的文件进行合并(cat)，并构建ML树。
# Date   : 2024-02-21
# Author : 孟庆瑶
# **************************************************

set -e

name=$0

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo "将不同文件夹中，文件名称一样的文件进行合并(cat)，并构建ML树。"
    echo -e "usage: bash $name [FOLDER ONE] [FOLDER TOW] [OUTPUT FOLDER]"
    echo -e "\t[The order of documents must not be incorrect.]"
    echo -e "\t        >>>> [请全部输入绝对路径] <<<<"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

path=$(pwd)


if [ -a $path/$3 ] ;then
    rm -rf $path/$3
fi

mkdir $path/$3

ids=$(ls $path/$1) && echo ${ids} > file.ids && sed -i 's/ /\n/g' file.ids

for i in $(cat file.ids) ;do
    cat $path/$1/${i} $path/$2/${i} > $path/$3/${i} 
    mafft --auto --inputorder $path/$3/${i} > $path/$3/${i}.mafft 
    trimal -in $path/$3/${i}.mafft -automated1 -fasta -out $path/$3/${i}_trimal.fa 
    FastTree $path/$3/${i}_trimal.fa > $path/$3/${i}.tree ;done

rm -rf file.ids

echo -e " \n >>>> Done <<<< \n "
echo -e "Results are stored in $path \n"
echo -e "Thank you for using our software! \n"
echo -e "Please contact us if you have any questions. \n"
echo -e "Email: <15877464851@163.com> \n"
echo -e "Sincerely, \n"
echo -e "MQY \n"