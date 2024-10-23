#!/bin/bash/
# **************************************************
# 找到共同的KO号
# Date   : 2024-01-04
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo "Modify the sequence name"
    echo -e "usage: bash eggnog_comm_kegg.sh [EGGNOG KEGG ONE] [EGGNOG KEGG TWO]"
    echo -e "\t[The order of documents must not be incorrect.]"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

cut -f 2 $1 > 1 && cut -f 2 $2 > 2 

comm -12 <(sort 1) <(sort 2) | uniq > comm_kegg.ids && rm -rf 1 2 

