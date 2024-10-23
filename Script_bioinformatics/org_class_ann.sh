#!/bin/bash/
# **************************************************
# 将注释的文件按照直系同源基因分开
# Date   : 2023-12-21
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ; then
    echo "`date '+[error] %D %T'`"
    echo "descriptive: SRA convert to fastq and quality assessment"
    echo "usage: bash org_class_ann.sh [FILE NAME LIST] [ANNOTATION] [DATABASE] [OUTPUT PREFIX]"
    echo -e "\t[The order of documents must not be incorrect.]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

if [ -a ids ] ;then
    rm -rf ids
fi

mkdir ids && cd ids
for i in `sed 's/\r//g' ../$1` ;do cat ../$3/$i* | grep ">" | sed 's/>//g' > $i\_ids.txt ;done
a=`ls | sed 's/\.txt//g'`
echo ${a} | xargs -n 1 > all_ids.txt

cd ../

if [ -a $4 ] ;then
    rm -rf $4
fi

mkdir $4 && cd $4
cat ../$2 | grep -v "##" > 1 && mv 1 $2
for i in `cat ../ids/all_ids.txt` ;do \
    for a in `cat ../ids/$i.*` ;do awk -v VAR=${a} '{if($1==VAR) {print$0}}' $2 >> $i\_anno.txt ;done \
    ;done

rm -rf $2
rm -rf ../ids

echo  -e "\n  处理完成 !!! \n  " 