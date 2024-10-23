#!/bin/bash/
# **************************************************
# RNA-seq,mapping and file exchange
# Date   : 2024-01-19
# Author : 孟庆瑶
# **************************************************

set -e

#  获取脚本名称并输出
name=$(basename $0) 
echo 
echo -e "   脚本名称: $name\n"

#  设置参数以及处理
OUTPREFIX="STAR_result"
INDEX="STAR_index"

#  设置封装
while getopts hg:i:d:G:o: opt
do
    case "$opt" in 
    h)
        echo
        echo -e "   脚本说明: [STAR + samtools] \n "
        echo 
        echo -e "使用说明: bash $name -g genome.fasta -G annotation.gft -i index prefix -d database -o output"
        echo "===================================================================================="
        echo -e "\t-g: 输入物种基因组 [.fasta]"
        echo -e "\t-i: index prefix 默认:[STAR_index]"
        echo -e "\t-d: 数据库文件夹"
        echo -e "\t-G: GTF注释文件"
        echo -e "\t-o: 输出文件名称 默认:[STAR_result]"
        echo -e "\n >>>> [输入文件文压缩为.gz，其他格式请参考官方说明文件] <<<< \n"
        echo "`date '+Date: %D %T'`"
        echo -e "\n             >>>> mqy <<<< \n "
        echo -e "             #        \"  " 
        echo -e "             # mm   mmm   "
        echo -e "             #\"  #    #   "
        echo -e "             #   #    #   "
        echo -e "             #   #  mm\#mm  "
        exit 0
        ;;
    g) 
        echo "输入-g文件为: $OPTARG"
        GEN=$OPTARG
        ;;
    i) 
        echo "输入-i文件为: $OPTARG"
        INDEX=$OPTARG
        ;;
    d) 
        echo "输入-d文件为: $OPTARG"
        DATABASE=$OPTARG
        ;;
    G) 
        echo "输入-G文件为: $OPTARG"
        GTF=$OPTARG
        ;;
    o) 
        echo "设置结果文件名前缀为: $OPTARG"
        OUTPREFIX=$OPTARG
        ;;
    *) 
        echo "未知参数: $opt"
        ;;
    esac
done

#####################-z 判断字符串是否为0
if [ -z $GEN ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    exit 1
fi

if [ -z $INDEX ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    exit 1
fi

if [ -z $DATABASE ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    exit 1
fi

if [ -z $OUTPREFIX ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    exit 1
fi

#  获取当前路径
path=$(pwd)

if [ -a $path/$OUTPREFIX ] ;then
    rm -rf $path/$OUTPREFIX
fi

mkdir $path/$OUTPREFIX

if [ -a $path/$INDEX ] ;then
    rm -rf $path/$INDEX
fi

mkdir $path/$INDEX

echo -e " \n **** 建立索引！[STAR] **** \n "

STAR --runThreadN 10 \
    --runMode genomeGenerate \
    --genomeDir $path/$INDEX/$INDEX \
    --genomeFastaFiles $path/$GEN \
    --sjdbGTFfile $path/$GTF \
    --sjdbOverhang 149 \
    --genomeSAindexNbases 9 # 根据基因组的大小，1711280db。


cd $path/$DATABASE

echo -e " \n 开始比对！[STAR] \n "

for i in `ls *gz | cut -d "_" -f1 | sort -u`
do
    echo -e " \n ****开始比对：$path/$DATABASE/${i}_1\t$path/$DATABASE/${i}_2 **** \n "
    STAR --runMode alignReads --runThreadN 20 \
        --genomeDir $path/$INDEX/$INDEX \
        --outTmpDir $path/$OUTPREFIX/${i}_tmp \
        --readFilesIn $path/$DATABASE/${i}*_1*.gz $path/$DATABASE/${i}*_2*.gz \
        --readFilesCommand zcat \
        --sjdbGTFfile $path/$GTF \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $path/$OUTPREFIX/${i} \
        --quantMode GeneCounts \
        --twopassMode Basic  \
        --outFilterMismatchNoverReadLmax 0
        # 一条read最多比对到10个不同位置，否则认为是unmapped
        #--outFilterMultimapNmax 500 \
        # 每对read最大的mismatch数目
        #--outFilterMismatchNmax 999 \
        # 每对read最大的mismatch占read长度的比例。例如PE100数据，则100*2*0.04=8bp
        #--outFilterMismatchNoverReadLmax 0.04 \
        # minimum intron size: genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion
        #--alignIntronMin 21 \
        #--alignIntronMax 1000000 \
        # PE数据之间最长gap为1000000bp
        #--alignMatesGapMax 1000000 \
        
    samtools index $path/$OUTPREFIX/${i}Aligned.sortedByCoord.out.bam
    samtools flagstat $path/$OUTPREFIX/${i}Aligned.sortedByCoord.out.bam > $path/$OUTPREFIX/${i}Aligned.sortedByCoord.out.bam.flagstat
done

echo -e " \n 干完了！ \n "