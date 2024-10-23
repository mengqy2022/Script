#!/bin/bash/
# **************************************************
# 共生菌基因预测以及去除假基因
# Date   : 2024-04-01
# Author : 孟庆瑶
# **************************************************

#  发生错误停止脚本
set -e

#  获取脚本名称并输出
name=$(basename $0) 
echo 
echo "   脚本名称: $name >>>> 未测试 <<<<"

#  设置参数以及处理
OUTPREFIX="Bacterial_results"

while getopts hi:g:o:v opt
do
    case "$opt" in 
    h)
        echo
        echo "脚本说明: 基因预测、假基因预测和获得去除假基因的CDs序列，(*￣︶￣)，东西不多，一字一字看！"
        echo
        echo "使用说明: bash $name -i genome_file_name -g genome_fasta -o out_prefix"
        echo "------------------------------------------------------------------------------------"
        echo -e "\t-i: 输入文本文件、包含基因文件名称、不带后缀、单列文件，一列为一个基因组；"
        echo -e "\t-g: 基因文件夹名称，文件夹中可以包含一个或者多个基因组；"
        echo -e "\t-o: 文件前缀, 默认:Bacterial_results；"
        echo -e "\t-v: 显示版本信息；"
        echo
        echo "`date '+Date: %D %T'`"
        echo -e "\n              >>>> mqy <<<< \n "
        echo -e "   如有疑问请联系: <15877464851@163.com>\n"
        echo "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "     !!!!    请将基因组放入同一个文件夹中。  !!!!"
        echo "     !!!!                                    !!!!"
        echo "     !!!!      -i文件中为基因组文件名称      !!!!"
        echo "     !!!!       -g为存放基因组的文件夹       !!!!"
        echo "     !!!!   -i内容和-g中的基因组名称要相同   !!!!"
        echo "     !!!!                                    !!!!"
        echo "     !!!!  输出结果文件包含:                 !!!!"
        echo "     !!!!               基因预测文件;        !!!!"
        echo "     !!!!               假基因预测文件;      !!!!"
        echo "     !!!!               去除假基因序列文件。 !!!!"
        echo "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo
        exit 0
        ;;
    i)
        echo "输入-i文件为: $OPTARG"
        FILENAME=$OPTARG
        ;;
    g)
        echo "输入-g文件为: $OPTARG"
        FOLDERNAME=$OPTARG
        ;;
    o)
        echo "设置输出文件名前缀为: $OPTARG"
        OUTPREFIX=$OPTARG
        ;;
    v)
        echo -e "\n版本信息: v1.0\n"
        echo -e "增加了并行运算的功能，在处理多个基因组时可以提高运行速度。\n"
        echo -e "增加了路径识别能力，可以自动识别输入文件和输出文件夹的路径。\n"
        exit 0
        ;;
    *)
        echo "未知参数: $opt"
        ;;
    esac
done

#####################-z 判断字符串是否为0
if [ -z $FILENAME ]
then
    echo -e "\n使用说明: bash $name -i genome_file_name -g genome_fasta -o out_prefix"
    echo -e "\n   请输入-h，查看帮助文档！\n"
    exit 1
fi

##################### 输出结果文件

path_run=$(pwd)

if [ -d $OUTPREFIX ]; then
    rm -rf $OUTPREFIX 
fi

mkdir $OUTPREFIX && cd $OUTPREFIX

# 识别输出路径
path_out=$(pwd)

############设置单个物种文件夹

cd $path_run

for i in $(cat $FIlENAME) ;do mkdir $FOLDERNAME/${i}  &&  mv $FOLDERNAME/${i}.* $FOLDERNAME/${i} ;done

############Prokka
source /home/mengqingyao/miniconda3/bin/activate prokka

cd $path_out

if [ -d "prokka"]; then
    rm -rf prokka
fi

mkdir prokka && cd prokka

## 最大并行数量
MAX_JOBS=8

run_prokka() {
    prokka $FOLDERNAME/${1}/${1}.* --cpus 0 --force --outdir ${1}_prokka --prefix ${1} --addmrna --compliant 200
}

for i in $(cat ../../$FILENAME) ; do 
    #-eq: 测试两个整数是否相等；比如 $A -eq $B
    #-ne: 测试两个整数是否不等；不等，为真；相等，为假；
    #-gt: 测试一个数是否大于另一个数；大于，为真；否则，为假；
    #-lt: 测试一个数是否小于另一个数；小于，为真；否则，为假；
    #-ge: 大于或等于
    #-le：小于或等于 
    
    while [ $(jobs | wc -l) -le ${MAX_JOBS} ]; do
        sleep 1;
    done

    echo "Prokka $i"
    run_prokka $i & 
done

wait
echo "Prokka 完成！"

############pseudofinder

conda activate pseudofinder

cd $path_out

if [ -d "pseudogene"]; then
    rm -rf pseudogene
fi

mkdir pseudogene && cd pseudogene

## 最大并行数量
MAX_JOBS=8

run_pseudofinder() {
    pseudofinder.py annotate --lenome ../prokka/${1}_prokka/${1}.gb* --outprefix ${1} \
    -di --database /data/mengqy/database/diamond_nr/diamond_makedb.dmnd --threads 16 -e 1e-15
}

for i in $(cat ../../$FILENAME) ; do 
    while [ $(jobs | wc -l) -le ${MAX_JOBS} ]; do
        sleep 1;
    done

    echo "pseudofinder $i"
    run_pseudofinder $i & 
done

wait
echo "pseudofinder 完成！"

#############获得去除假基因的序

cd $path_out

if [ -d "remove_pseudogene"]; then
    rm -rf remove_pseudogene
fi

mkdir remove_pseudogene && cd remove_pseudogene

## 最大并行数量
MAX_JOBS=8

run_remove_pseudogene() {
    cat ../prokka/${1}_prokka/${1}.ffn | grep ">" | cut -c 2- | cut -d " " -f 1 > ${1}_cds.ids && \
    awk '{if(match($0,"old_locus_tag=")) {print substr($0,RSTART+RLENGTH) }}' ../pseudogene/${1}_pseudos.gff \
    | sed 's/,/\t/g' | awk '{print $1}' | sed '/^$/d' >> ${1}_pseudofene.ids && comm -23 <(sort ${1}_cds.ids) <(sort ${1}_pseudofene.ids) > ${1}_remove.ids && \
    python3 /home/mengqingyao/Script_mqy/get_ids_fasta.py ../prokka/${1}_prokka/${1}.faa ${1}_remove.ids ${1}_remove_ids.faa && \
    python3 /home/mengqingyao/Script_mqy/get_ids_fasta.py ../prokka/${1}_prokka/${1}.ffn ${1}_remove.ids ${1}_remove_ids.ffn
}

for i in $(cat ../../$FILENAME) ; do 
    while [ $(jobs | wc -l) -le ${MAX_JOBS} ]; do
        sleep 1;
    done

    echo "获得去除假基因的序 $i"
    run_remove_pseudogene $i & 
done

wait
echo "获得去除假基因的序 完成！"

#############输出结果

echo "输出结果文件为: $path_out"
echo "输出结果中包含:"
echo "               基因预测文件;"
echo "               假基因预测文件;"
echo "               去除假基因序列文件。"
echo ""
echo "   运行结束！"
echo "   感谢使用本脚本！"
echo "   版本信息: v1.0"
echo "   日期: `date '+%Y-%m-%d %H:%M:%S'`"
echo "   作者: 孟庆瑶"
echo "   如有疑问请联系: <15877464851@163.com>"