#!/bin/bash/
# **************************************************
# 为了进行下游分析，统一序列的名称 fastq
# Date   : 2023-12-20
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo "Modify the sequence name"
    echo -e "usage: bash reads_name_mod_FQ.sh [MOD FILE NAME] [R1_fastq] [R2_fastq]"
    echo -e "\t[The order of documents must not be incorrect.]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

awk -v VAR '$@/$@/{pritn "@VAR" ++i; next} {print}' ${1}_R1.fastq > ${1}_R1mode.fastq &
awk -v VAR '$@/$@/{pritn "@VAR" ++i; next} {print}' ${1}_R2.fastq > ${1}_R2mode.fastq &