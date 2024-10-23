#!/bin/bash/
# **************************************************
# 测序数据处理，去接头和低质量过滤，以及mapping
# Date   : 2024-01-19
# Author : 孟庆瑶
# **************************************************

set -e

#  第一层输出程序

if [ -z $1 ]
then  #  如果string为空，则为真
    #  获取脚本名称并输出
    name=$(basename $0) 
    echo 
    echo "   脚本名称: $name"
    echo -e " \n `date '+[啥也没输入] %D %T'` \n "
    echo -e "   脚本说明: [Remove_adapter + Global_alignment] \n "
    echo -e "使用: bash $name [INPUT NAME] \n"
    echo -e "\t [Remove_adapter] \n"
    echo -e "\t [Global_alignment] "
    echo -e "\n             >>>> mqy <<<< \n "
    echo -e "             #        \"  " 
    echo -e "             # mm   mmm   "
    echo -e "             #\"  #    #   "
    echo -e "             #   #    #   "
    echo -e "             #   #  mm\#mm  "
    exit 1

    elif [$1=="Global_alignment"]
    then
        name="sym_mapping_v_2.bash Global_alignment"
        echo 
        echo "   脚本名称: $name"
        echo
        #  设置参数以及处理
        OUTPREFIX="result"
        INDEX="bacteria"

        #  设置封装
        while getopts hg:i:d:o: opt
        do
            case "$opt" in 
            h)
                echo
                echo -e "   脚本说明: [全局比对(hista2) + 文件转换(samtools)] \n "
                echo 
                echo -e "使用说明: bash sym_mapping_v_2.bash Remove adapter -g genome.fasta -i index prefix -d database -o output"
                echo "===================================================================================="
                echo -e "\t-g: 输入物种基因组 [.fasta]"
                echo -e "\t-i: index prefix 默认:[bacteria]"
                echo -e "\t-d: 数据库文件夹"
                echo -e "\t-o: 输出文件名称 默认:[result]"
                echo
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
                echo "输入-r文件为: $OPTARG"
                DATABASE=$OPTARG
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

        #  建立索引
        hisat2-build $path/$GEN $path/$INDEX/$INDEX

        echo -e " \n 开始比对！ \n "

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
            # ((count++))
            # #每10个样本等待任务完成再继续
            # if [ $count -eq 2 ];then 
            #     wait 
            #     count=0 
            # fi
        done
        echo -e " \n 干完了！ \n "


fi

















