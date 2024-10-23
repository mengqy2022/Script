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
OUTPREFIX="result"
INDEX="index"

#  设置封装
while getopts hg:i:d:c:G:o: opt
do
    case "$opt" in 
    h)
        echo
        echo -e "   脚本说明: [Global_alignment + File_conversion] \n "
        echo 
        echo -e "使用说明: bash $name -g genome.fasta -c hisat2,STAR -G annotation.gft -i index prefix -d database -o output"
        echo "===================================================================================="
        echo -e "\t-g: 输入物种基因组 [.fasta]"
        echo -e "\t-i: index prefix 默认:[bacteria]"
        echo -e "\t-d: 数据库文件夹"
        echo -e "\t-c: 选择进行mapping的软件"
        echo -e "\t-G: GTF注释文件 [TBtools转换来的GTF文件可以正常运行]"
        echo -e "\t-o: 输出文件名称 默认:[result]"
        echo -e "\n >>>> [当-c为STAR时，才需要gft文件] <<<< \n"
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
    c) 
        echo "输入-c类型为: $OPTARG"
        CLASS=$OPTARG
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

if [ -z $CLASS ]
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

if [$CLASS=="hisat2"]
then
    echo -e " \n 建立索引！[hisat2] \n "
    #  建立索引
    hisat2-build $path/$GEN $path/$INDEX/$INDEX

    echo -e " \n 开始比对！[hisat2] \n "

    #  比对
    cd $path/$DATABASE

    for i in `ls *gz | cut -d "_" -f1 | sort -u`
    do
        cd $path/$DATABASE
        echo -e " \n ****开始比对：$path/$DATABASE/${i}_1\t$path/$DATABASE/${i}_2 \n**** "
        hisat2 -p 20 -x $path/$INDEX/$INDEX -1 $path/$DATABASE/${i}*_1P*.gz -2 $path/$DATABASE/${i}*_2P*.gz -S $path/$OUTPREFIX/${i}_hisat.sam #&

        cd $path/$OUTPREFIX
        echo -e " \n ****运行中：$path/$OUTPREFIX/${i}_hisat.bam 转换为 $path/$OUTPREFIX/${i}_hisat.sam \n**** "
        samtools view -bS $path/$OUTPREFIX/${i}_hisat.sam -o $path/$OUTPREFIX/${i}_hisat.bam && rm -rf $path/$OUTPREFIX/${i}_hisat.sam #&
        samtools sort $path/$OUTPREFIX/${i}_hisat.bam -O BAM > $path/$OUTPREFIX/${i}_hisat_sort.bam && rm -rf $path/$OUTPREFIX/${i}_hisat.bam
        samtools index $path/$OUTPREFIX/${i}_hisat_sort.bam
        samtools flagstat $path/$OUTPREFIX/${i}_hisat_sort.bam > $path/$OUTPREFIX/${i}_hisat_sort.flagstat
    done
else

    echo -e " \n 建立索引！[STAR] \n "
    STAR --runThreadN 10 --runMode genomeGenerate  \
    --genomeDir $path/$INDEX/$INDEX \
    --genomeFastaFiles $path/$GEN \
    --sjdbGTFfile $path/$GTF \
    --sjdbOverhang 149 \
    --genomeSAindexNbases 9   #  根据基因组的大小，1711280db。
    
    cd $path/$DATABASE

    echo -e " \n 开始比对！[STAR] \n "

    for i in `ls *gz | cut -d "_" -f1 | sort -u`
    do
        echo -e " \n ****开始比对：$path/$DATABASE/${i}_1\t$path/$DATABASE/${i}_2 \n**** "

        STAR --runMode alignReads --runThreadN 20 \
        --genomeDir $path/$INDEX/$INDEX \
        --readFilesIn $path/$DATABASE/${i}*_1P*.gz $path/$DATABASE/${i}*_2P*.gz \
        --readFilesCommand zcat \
        --sjdbGTFfile $path/$GTF \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $path/$OUTPREFIX/${i} \
        --quantMode GeneCounts \
        --twopassMode Basic --twopass1readsN -1 --seedSearchStartLmax 15 --outSJfilterOverhangMin 15 8 8 8 --outFilterMismatchNoverReadLmax 0.1

        samtools index $path/$OUTPREFIX/${i}Aligned.sortedByCoord.out.bam
        samtools flagstat $path/$OUTPREFIX/${i}Aligned.sortedByCoord.out.bam > $path/$OUTPREFIX/${i}Aligned.sortedByCoord.out.bam.flagstat
    done
fi

echo -e " \n 干完了！ \n "