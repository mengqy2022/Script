#!/bin/bash/
# **************************************************
# 共生菌blat
# Date   : 2024-01-12
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then  #  如果string为空，则为真
    echo -e " \n `date '+[啥也没输入] %D %T'` \n "
    echo "本脚本用于将fastq转变fasta文件的blat比对"
    echo "使用: bash sym_blat.sh [IDlist] [SPECIES FILE] [database] [OUTPUT FPREFILE]"
    echo -e "\n             >>>> mqy <<<< \n "
    echo -e "  文件输入顺序不能乱 \n "
    exit 1
fi

#  获取当前路径
path=$(pwd)

if [ -a $path$4 ] ;then
    rm -rf $path$4
fi

echo -e " \n 开始比对！ \n "

source /home/mengqingyao/miniconda3/bin/activate blat

mkdir $4

for ids in `cat $1`;do blat  $path/$2 $path/$3/${ids}.* -minIdentity=99 -out=blast8 $path/$4/${ids}\_blat.txt && \
    python3 /home/mengqingyao/biosoftware/Script_mqy/blast_tophit.py $path/$4/${ids}\_blat.txt $path/$4/${ids}\_blat_tophit.txt && \
    awk '{if ($4>=50) print $0}' $path/$4/${ids}\_blat_tophit.txt | less > $path/$4/1 && \
    awk '{if ($5<=2) print $0}' $path/$4/1 > $path/$4/2 && mv $path/$4/2 $path/$4/${ids}\_blat_tophit_filter.txt && rm -rf $path/$4/1 && \
    cat $path/$4/${ids}\_blat_tophit_filter.txt | cut -f 1 > $path/$4/${ids}\_blat_tophit_filter.ids ;done

echo -e " \n 干完了！ \n "